#!/usr/bin/env bash
# alt-depth-ab.sh — measure the depth a reference with ALT contigs costs.
#
# Usage: scripts/ci/alt-depth-ab.sh <work_dir>
#
# Maps the fixture's CYP2D6-region and MHC-region read pairs to two whole
# references with step 02's aligner and preset (minimap2 -x sr, not
# ALT-aware), and prints mosdepth's mean depth over CYP2D6, CYP2D7, their
# flanks, HLA-A and HLA-B, for all reads and for MAPQ >= 1 (callers ignore
# MAPQ 0 reads):
#   no-ALT    NCBI's GRCh38 no-ALT analysis set, today's default (195 sequences)
#   with-ALT  the Broad hg38 FASTA the pipeline used before (3,366 sequences,
#             ALT, HLA and decoy contigs)
# It also counts the primary alignments each reference places on ALT or HLA
# contigs, which a caller or T1K reading chr6 and chr22 never sees.
#
# The whole-genome minimap2 index is built in parts of IDX_PART bases (-I),
# so it fits a 16 GB runner, and mapped with --split-prefix, which merges the
# parts into the MAPQ one index would give. Needs Docker, gh (GH_TOKEN), curl,
# md5sum and about 25 GB of disk. Writes ab.tsv and ab.md into <work_dir>,
# and ab.md to the job summary.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"

OUT_ARG=${1:?Usage: $0 <work_dir>}
mkdir -p "$OUT_ARG"
W="$(cd "$OUT_ARG" && pwd)"
THREADS=${THREADS:-4}
IDX_PART=${IDX_PART:-1G}
GH_REPO=${GITHUB_REPOSITORY:-GeiserX/Personal-Genome-Pipeline}
TAG=${FIXTURE_TAG:-$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")}
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}

NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
NOALT_NAME=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
BROAD_URL=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta

# The read sets: the fixture's slices of the two loci (build-fixture.sh).
LOCI=("cyp2d chr22:42000000-42300000" "mhc chr6:29900000-33100000")
# Regions, GRCh38, 0-based BED. CYP2D6 and CYP2D7 as Cyrius 1.1.1 (step 21)
# defines them in data/CYP2D6_region_38.bed (CYP2D6 with REP6); the flanks are
# two 50 kb stretches of the CYP2D slice outside the CYP2D6-CYP2D8 cluster;
# HLA-A and HLA-B the gene windows build-fixture.sh sends to VEP.
REGIONS_BED=$'chr22\t42050000\t42100000\tCYP2D flanks\nchr22\t42123192\t42132032\tCYP2D6\nchr22\t42139676\t42145745\tCYP2D7\nchr22\t42200000\t42250000\tCYP2D flanks\nchr6\t29941259\t29945884\tHLA-A\nchr6\t31353871\t31357188\tHLA-B'
REGION_ORDER=("CYP2D6" "CYP2D7" "CYP2D flanks" "HLA-A" "HLA-B")

for tool in docker gh curl md5sum awk sort; do
  command -v "$tool" >/dev/null || { echo "ERROR: ${tool} is required" >&2; exit 1; }
done

in_image() {
  local image=$1; shift
  docker run --rm -i -u "$(id -u):$(id -g)" -e HOME=/tmp --memory 14g -v "${W}:/w" -w /w "$image" "$@"
}
sam() { in_image "$SAMTOOLS_IMAGE" samtools "$@"; }

# --- 1. Reads -------------------------------------------------------------------
echo "=== Fixture ${TAG}: read pairs of ${LOCI[*]} ==="
deadline=$(( $(date +%s) + 90 * 60 ))
until gh release download "$TAG" -R "$GH_REPO" -D "${W}/fixture" --clobber \
        -p HG002_slice.bam -p HG002_slice.bam.bai -p SHA256SUMS 2>/dev/null; do
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "ERROR: release ${TAG} not downloadable after 90 minutes (the E2E workflow's build-fixture job publishes it)." >&2
    exit 1
  fi
  echo "Release ${TAG} is not published yet; waiting."
  sleep 60
done
(cd "${W}/fixture" && grep -E ' HG002_slice\.bam(\.bai)?$' SHA256SUMS | sha256sum -c -)
for l in "${LOCI[@]}"; do
  read -r name region <<< "$l"
  sam view -b -o "/w/${name}.bam" /w/fixture/HG002_slice.bam "$region"
  sam collate -u -O "/w/${name}.bam" "/w/${name}.collate" \
    | sam fastq -n -1 "/w/${name}_R1.fastq.gz" -2 "/w/${name}_R2.fastq.gz" -0 /dev/null -s /dev/null -
  rm -f "${W}/${name}.bam"
  echo "  ${name}: $(( $(gzip -dc "${W}/${name}_R1.fastq.gz" | wc -l) / 4 )) read pairs from ${region}"
done
printf '%s\n' "$REGIONS_BED" > "${W}/regions.bed"

# --- 2. Each reference: download, check, index, map, measure --------------------
# get_ref KIND: write ${W}/ref.fa and print where it came from and its check.
get_ref() {
  local want got url
  case "$1" in
    no-ALT)
      url="${NCBI}/${NOALT_NAME}"
      curl -fsSL --retry 5 --retry-delay 10 -o "${W}/ref.fa.gz" "$url"
      want=$(curl -fsSL --retry 5 "${NCBI}/md5checksums.txt" | awk -v f="./${NOALT_NAME}" '$2 == f {print $1}')
      got=$(md5sum "${W}/ref.fa.gz" | awk '{print $1}')
      [ -n "$want" ] && [ "$want" = "$got" ] || { echo "ERROR: ${NOALT_NAME} md5 ${got}, NCBI lists '${want}'" >&2; return 1; }
      gzip -dc "${W}/ref.fa.gz" > "${W}/ref.fa"
      rm -f "${W}/ref.fa.gz"
      echo "${url} (md5 ${got}, as NCBI's md5checksums.txt lists it)" ;;
    with-ALT)
      url="$BROAD_URL"
      curl -fsSL --retry 5 --retry-delay 10 -o "${W}/ref.fa" "$url"
      # The bucket reports the object's md5 as base64 in x-goog-hash.
      want=$(curl -fsSI "$url" | tr -d '\r' | awk -F'md5=' 'tolower($0) ~ /^x-goog-hash:/ && NF > 1 {print $2}' \
        | base64 -d | od -An -v -tx1 | tr -d ' \n')
      got=$(md5sum "${W}/ref.fa" | awk '{print $1}')
      [ -n "$want" ] && [ "$want" = "$got" ] || { echo "ERROR: Broad hg38 md5 ${got}, the bucket reports '${want}'" >&2; return 1; }
      echo "${url} (md5 ${got}, as the bucket's x-goog-hash reports it)" ;;
  esac
}

: > "${W}/ab.tsv"
: > "${W}/placement.tsv"
: > "${W}/sources.txt"
for l in "${LOCI[@]}"; do read -r name _ <<< "$l"; : > "${W}/ab_${name}.tsv"; done
for kind in no-ALT with-ALT; do
  echo "=== ${kind} reference ==="
  src=$(get_ref "$kind")
  echo "  ${src}"
  printf '%s: %s; %s sequences\n' "$kind" "$src" "$(grep -c '^>' "${W}/ref.fa")" >> "${W}/sources.txt"
  echo "  minimap2 -x sr index in parts of ${IDX_PART} bases"
  start=$(date +%s)
  in_image "$MINIMAP2_IMAGE" minimap2 -x sr -I "$IDX_PART" -t "$THREADS" -d /w/ref.mmi /w/ref.fa
  echo "  index: $(( $(date +%s) - start )) s, $(du -h "${W}/ref.mmi" | cut -f1)"
  rm -f "${W}/ref.fa"
  for l in "${LOCI[@]}"; do
    read -r name _ <<< "$l"
    bam="${kind}_${name}.bam"
    in_image "$MINIMAP2_IMAGE" minimap2 -t "$THREADS" -a -x sr --split-prefix "/w/split_${kind}_${name}" \
        /w/ref.mmi "/w/${name}_R1.fastq.gz" "/w/${name}_R2.fastq.gz" \
      | sam sort -@ "$THREADS" -o "/w/${bam}" -
    sam index "/w/${bam}"
    # Primary alignments, and how many of them sit on an ALT or HLA contig.
    sam view -F 0x904 "/w/${bam}" | awk -F'\t' -v k="$kind" -v n="$name" '
      { total++ } $3 ~ /_alt$/ || $3 ~ /^HLA-/ { alt++ }
      END { printf "%s\t%s\t%d\t%d\n", k, n, total, alt }' >> "${W}/placement.tsv"
    for q in 0 1; do
      in_image "$MOSDEPTH_IMAGE" mosdepth -n -t 2 -Q "$q" -b /w/regions.bed "/w/${kind}_${name}_q${q}" "/w/${bam}"
      # chrom start end name mean -> a length-weighted mean per region name.
      gzip -dc "${W}/${kind}_${name}_q${q}.regions.bed.gz" | awk -F'\t' -v k="$kind" -v q="$q" '
        { l = $3 - $2; s[$4] += $5 * l; n[$4] += l }
        END { for (r in s) printf "%s\t%s\t%d\t%.2f\n", r, k, q, s[r] / n[r] }' >> "${W}/ab_${name}.tsv"
    done
  done
  rm -f "${W}/ref.mmi"
done

# A region's depth comes from the read set of its own locus.
for l in "${LOCI[@]}"; do
  read -r name _ <<< "$l"
  case "$name" in
    cyp2d) awk -F'\t' '$1 ~ /^CYP2D/' "${W}/ab_${name}.tsv" ;;
    mhc) awk -F'\t' '$1 ~ /^HLA-/' "${W}/ab_${name}.tsv" ;;
  esac
done >> "${W}/ab.tsv"

# Every region has a depth for both references and both MAPQ floors, and
# every all-reads depth is above 0: a short or empty table is a broken run,
# not a measurement, so it fails before any table is written. (A MAPQ >= 1
# depth of 0 is a possible result.)
rows=$(grep -c . "${W}/ab.tsv" || true)
zero=$(awk -F'\t' '$3 == 0 && $4 <= 0' "${W}/ab.tsv" | wc -l | tr -d ' ')
if [ "$rows" -ne $(( ${#REGION_ORDER[@]} * 4 )) ] || [ "$zero" -ne 0 ]; then
  echo "ERROR: ab.tsv has ${rows} depths (want $(( ${#REGION_ORDER[@]} * 4 ))) and ${zero} regions with no reads at all" >&2
  exit 1
fi

# --- 3. Table -------------------------------------------------------------------
{
  echo "### Depth with and without ALT contigs (fixture ${TAG}, HG002 at about 30x)"
  echo
  echo "Mean depth from mosdepth. Callers ignore reads with MAPQ 0, so the MAPQ >= 1 column is the depth they see."
  echo
  echo "| Region | Reference | Depth, all reads | Depth, MAPQ >= 1 | MAPQ >= 1 / all |"
  echo "|---|---|---|---|---|"
  for r in "${REGION_ORDER[@]}"; do
    for kind in no-ALT with-ALT; do
      awk -F'\t' -v r="$r" -v k="$kind" '
        $1 == r && $2 == k { d[$3] = $4 }
        END { printf "| %s | %s | %.1f | %.1f | %s |\n", r, k, d[0], d[1], (d[0] > 0 ? sprintf("%.2f", d[1] / d[0]) : "n/a") }' "${W}/ab.tsv"
    done
  done
  echo
  echo "| Reference | Read set | Primary alignments | On ALT or HLA contigs |"
  echo "|---|---|---|---|"
  awk -F'\t' '{printf "| %s | %s | %d | %d |\n", $1, $2, $3, $4}' "${W}/placement.tsv"
  echo
  echo "References:"
  echo
  sed 's/^/- /' "${W}/sources.txt"
  echo
  echo "Images: ${MINIMAP2_IMAGE}, ${SAMTOOLS_IMAGE}, ${MOSDEPTH_IMAGE}. Index built with -I ${IDX_PART}, mapped with --split-prefix."
} > "${W}/ab.md"
cat "${W}/ab.md"
cat "${W}/ab.md" >> "$SUMMARY"
