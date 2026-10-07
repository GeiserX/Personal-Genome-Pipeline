#!/usr/bin/env bash
# write_manifest.sh: record what produced a sample's outputs.
#
# Usage: GENOME_DIR=... bin/write_manifest.sh <sample> <written_by> [declared_sex]
#
# Writes ${GENOME_DIR}/<sample>/run_manifest.tsv, three tab-separated columns
# (section, key, value):
#   run      written_utc, written_by, git_commit, declared_sex
#   version  every NAME=VALUE line of versions.env
#   image    each *_IMAGE reference and the digest Docker resolved it to
#            ("not pulled" when the image is not on this machine)
#   data     clinvar_file_date (##fileDate of the ClinVar file), vep_cache,
#            pcgr_bundle, pypgx_bundle, hla_database (newest IPD-IMGT/HLA
#            release named in t1k_idx/hlaidx/hla.dat), pgs:<ID> (the header
#            lines of each PGS scoring file)
#
# run-all.sh calls it when a run starts and again just before the reports, so
# an image a step pulled during the run gets its digest; the two report scripts
# call it when the sample has no manifest yet. The HTML report prints it in its footer.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample> <written_by> [declared_sex]}
WRITTEN_BY=${2:?Usage: $0 <sample> <written_by> [declared_sex]}
DECLARED_SEX=${3:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=../scripts/lib/common.sh
. "$(dirname "$0")/../scripts/lib/common.sh"
validate_sample "$SAMPLE"

OUT="${GENOME_DIR}/${SAMPLE}/run_manifest.tsv"
mkdir -p "${GENOME_DIR}/${SAMPLE}"
TMP="${OUT}.tmp"
: > "$TMP"
row() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "$TMP"; }

row run written_utc "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
row run written_by "$WRITTEN_BY"
commit="unknown (not a git checkout)"
if git -C "$PGP_ROOT" rev-parse --verify -q HEAD >/dev/null 2>&1; then
  commit=$(git -C "$PGP_ROOT" rev-parse HEAD)
  if [ -n "$(git -C "$PGP_ROOT" status --porcelain --untracked-files=no 2>/dev/null)" ]; then
    commit="${commit} (with uncommitted changes)"
  fi
fi
row run git_commit "$commit"
[ -z "$DECLARED_SEX" ] || row run declared_sex "$DECLARED_SEX"

# versions.env as written (comments and blank lines dropped)
while IFS= read -r line || [ -n "$line" ]; do
  [[ "$line" =~ ^([A-Za-z0-9_]+)=(.*)$ ]] || continue
  value=${BASH_REMATCH[2]}
  value=${value%%#*}
  value=${value%"${value##*[![:space:]]}"}
  value=${value#\"}
  value=${value%\"}
  row version "${BASH_REMATCH[1]}" "$value"
done < "${PGP_ROOT}/versions.env"

# Every image, with the digest this machine holds for it
while IFS= read -r line || [ -n "$line" ]; do
  [[ "$line" =~ ^([A-Z0-9_]+_IMAGE)= ]] || continue
  name=${BASH_REMATCH[1]}
  ref=${!name:-}
  [ -n "$ref" ] || continue
  digest=$("$CONTAINER_ENGINE" image inspect --format '{{join .RepoDigests " "}}' "$ref" 2>/dev/null || true)
  if [ -z "$digest" ]; then
    id=$("$CONTAINER_ENGINE" image inspect --format '{{.Id}}' "$ref" 2>/dev/null || true)
    digest=${id:+local image ${id}}
  fi
  row image "$ref" "${digest:-not pulled}"
done < "${PGP_ROOT}/versions.env"

# header_lines FILE [MAX]: the leading '#' lines of a (gzipped) text file.
# Its consumer must read to the end (no early exit), or the pipe fails under pipefail.
header_lines() {
  { gzip -dcf "$1" 2>/dev/null || true; } | awk -v max="${2:-200}" '/^#/ {print; if (++n >= max) exit; next} {exit}'
}

clinvar_date="not found"
for f in clinvar/clinvar_pathogenic_chr.vcf.gz clinvar/clinvar_chr.vcf.gz clinvar/clinvar.vcf.gz; do
  [ -f "${GENOME_DIR}/${f}" ] || continue
  # The reader goes to the end of the header: an awk that exits at the first
  # match would end the pipe early and stop this script with SIGPIPE (141).
  d=$(header_lines "${GENOME_DIR}/${f}" | awk -F= '/^##fileDate=/ && d == "" {d = $2} END {print d}')
  if [ -n "$d" ]; then clinvar_date="${d} (${f})"; break; fi
done
row data clinvar_file_date "$clinvar_date"

caches=""
for d in "${GENOME_DIR}"/vep_cache/homo_sapiens/*_GRCh38; do
  [ -f "${d}/info.txt" ] && caches="${caches:+${caches}, }$(basename "$d")"
done
row data vep_cache "${caches:-none installed}"
if [ -d "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data" ]; then
  row data pcgr_bundle "${PCGR_DATA_BUNDLE} (installed)"
else
  row data pcgr_bundle "${PCGR_DATA_BUNDLE} (not installed)"
fi
bundle="${GENOME_DIR}/reference/pypgx-bundle"
if [ -d "$bundle" ]; then
  ref=$(git -C "$bundle" describe --tags --always 2>/dev/null || echo "present, not a git checkout")
  row data pypgx_bundle "$ref"
else
  row data pypgx_bundle "not installed"
fi
hla="${GENOME_DIR}/t1k_idx/hlaidx/hla.dat"
if [ -f "$hla" ]; then
  rel=$(awk '/^DT/ && match($0, /Rel\. [0-9]+\.[0-9]+\.[0-9]+/) {print substr($0, RSTART + 5, RLENGTH - 5)}' "$hla" \
    | sort -t. -k1,1n -k2,2n -k3,3n | tail -1)
  row data hla_database "IPD-IMGT/HLA ${rel:-release not found in hla.dat}"
else
  row data hla_database "not installed (step 8 builds it)"
fi
for f in "${GENOME_DIR}"/prs_scores/*.txt.gz; do
  [ -f "$f" ] || continue
  id=$(basename "$f" .txt.gz)
  head=$(header_lines "$f" 60 | grep -E '^#(pgs_id|pgs_name|trait_reported|format_version|HmPOS_build|HmPOS_date|variants_number)=' \
    | sed 's/^#//' | paste -sd ';' - || true)
  row data "pgs:${id}" "${head:-no header lines}"
done

mv "$TMP" "$OUT"
echo "Run manifest: ${OUT} ($(grep -c '^image' "$OUT" || true) images)"
