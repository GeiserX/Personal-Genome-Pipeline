#!/usr/bin/env bash
# 26b-ancestry-fine.sh — Fine-grained ancestry: admixture fractions and a chromosome painting
# Usage: ./scripts/26b-ancestry-fine.sh <sample_name>
#
# Opt-in. Step 26 names the most similar of five continental groups. This
# step goes further against an open reference panel of 80 populations
# (gnomAD's HGDP + 1000 Genomes, phased), installed by
# scripts/setup.sh --ancestry-fine-panel:
#   1. genotypes the panel's SNVs from the gVCF (so the sites where you match
#      the reference count), or from a VCF that holds reference calls;
#   2. ADMIXTURE projects you on the allele frequencies of the panel's groups
#      at two levels (continental and regional groups; 5 and 17 on the
#      default panel, which leaves Oceania out), and a block bootstrap of the
#      SNVs gives the spread of each fraction;
#   3. Beagle phases your genotypes against the panel;
#   4. RFMix paints each segment of each haplotype with the continental group
#      it is most similar to;
#   5. bin/ancestry_fine_report.py writes the fractions table, the painting
#      (BED and PNG) and a plain-text summary.
# A group name is ours: a set of reference populations, not a nationality.
# docs/26b-ancestry-fine.md says where the resolution ends. ADMIXTURE is free
# for academic use; RFMix for academic research use only, and any other use,
# a personal one included, needs a licence from Stanford.
#
# Without the panel the step says so in one line and exits 0.
# Env: ANCESTRY_FINE_DIR (another panel folder inside GENOME_DIR, or none to
# skip), THREADS (default 4), ANCESTRY_FINE_BOOTSTRAP (bootstrap replicates,
# default 100; 0 for none).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image BCFTOOLS_IMAGE PYTHON_IMAGE ADMIXTURE_IMAGE BEAGLE_IMAGE RFMIX_IMAGE

PANEL=${ANCESTRY_FINE_DIR:-${GENOME_DIR}/reference/ancestry_fine/${ANCESTRY_FINE_PANEL}}
S="${GENOME_DIR}/${SAMPLE}"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"
GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
OUT="${S}/ancestry_fine"
WORK="${OUT}/work"
BOOT=${ANCESTRY_FINE_BOOTSTRAP:-100}

if [ "${ANCESTRY_FINE_DIR:-}" = none ] || [ ! -s "${PANEL}/model/DONE" ]; then
  echo "Step 26b skipped: no fine-grained ancestry panel at ${PANEL}; install it with scripts/setup.sh --ancestry-fine-panel ${GENOME_DIR}"
  exit 0
fi
P=$(cpath "$PANEL")
for f in "$VCF" "${VCF}.tbi" "$REF_FASTA"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
case "$BOOT" in ''|*[!0-9]*) echo "ERROR: ANCESTRY_FINE_BOOTSTRAP must be a whole number, got '${BOOT}'" >&2; exit 2 ;; esac

echo "============================================"
echo "  Step 26b: Fine-grained ancestry (${ANCESTRY_FINE_PANEL})"
echo "  Tools: ADMIXTURE, Beagle, RFMix"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUT}"
echo "============================================"
echo ""

rm -rf "$WORK" "${OUT}/${SAMPLE}_ancestry_fine".* "${OUT}/${SAMPLE}_painting".*
mkdir -p "${WORK}/admixture" "${WORK}/phase" "${WORK}/paint"
W=$(cpath "$WORK")

# --- 1. The panel's SNVs in this sample ---------------------------------------------
# From the gVCF as step 25 does: gvcf2vcf expands each reference block over
# a panel site into a 0/0 record, -T keeps the panel's positions,
# --trim-alt-alleles drops the <*> allele and a no-call is dropped. The awk
# keeps a record only when its REF is the panel's and its ALT is the panel's
# (a 0/0 has none: it gets the panel's), so every record is a genotype of the
# panel's two alleles, and counts them per chromosome with the first and the
# last position.
GENO="${WORK}/genotypes.vcf"
# The awk below reads the first sample column only, so a file with more than
# one sample would give that person's ancestry under this sample's name: the
# inner bash stops first unless the file holds exactly one sample.
if [ -f "$GVCF" ] && [ -f "${GVCF}.tbi" ]; then
  INPUT_KIND=gvcf
  echo "[1/5] Genotyping the panel's SNVs from the gVCF (${GVCF})..."
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  CMD=(bash -euo pipefail -c 'n=$(bcftools query -l "$3" | wc -l | tr -d " "); [ "$n" -eq 1 ] || { echo "ERROR: $3 holds $n samples; step 26b reads a one-sample file (keep one with bcftools view -s <sample>)." >&2; exit 1; }; bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Ov'
       _ "$REF_FASTA_C" "${P}/panel/positions.tsv" "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz")
else
  INPUT_KIND=vcf
  echo "[1/5] No gVCF beside the VCF: genotyping the panel's SNVs from the VCF (${VCF})..."
  # shellcheck disable=SC2016  # $1 and $2 belong to the inner bash
  CMD=(bash -euo pipefail -c 'n=$(bcftools query -l "$2" | wc -l | tr -d " "); [ "$n" -eq 1 ] || { echo "ERROR: $2 holds $n samples; step 26b reads a one-sample file (keep one with bcftools view -s <sample>)." >&2; exit 1; }; bcftools view -T "$1" -Ov "$2"' _ "${P}/panel/positions.tsv" "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz")
fi
run_in --cpus 2 --memory 4g "$BCFTOOLS_IMAGE" "${CMD[@]}" \
  | awk -F'\t' -v OFS='\t' -v sites="${PANEL}/panel/sites.tsv" -v counts="${WORK}/counts.tsv" '
      BEGIN { while ((getline l < sites) > 0) { split(l, f, "\t"); ref[f[1] ":" f[2]] = f[3]; alt[f[1] ":" f[2]] = f[4] } }
      /^##/ { if ($0 !~ /^##(INFO|FILTER|FORMAT)=/) print; next }
      /^#CHROM/ { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; NF = 10; print; next }
      {
        k = $1 ":" $2
        if (!(k in ref) || (k in seen) || $4 != ref[k]) next
        if ($5 == ".") $5 = alt[k]
        else if ($5 != alt[k]) next
        n = split($9, fmt, ":"); split($10, val, ":"); gt = "."
        for (i = 1; i <= n; i++) if (fmt[i] == "GT") gt = val[i]
        if (gt !~ /^[01][\/|][01]$/) next
        seen[k] = 1; if (!($1 in per)) first[$1] = $2
        per[$1]++; last[$1] = $2; if (gt ~ /^0[\/|]0$/) homref++
        print $1, $2, ".", $4, $5, ".", ".", ".", "GT", gt
      }
      END { for (c in per) print c "\t" per[c] "\t" first[c] "\t" last[c] > counts; print "homref\t" homref + 0 > counts }' > "$GENO"
N_GENO=$(awk -F'\t' '$1 != "homref" {n += $2} END {print n + 0}' "${WORK}/counts.tsv")
N_HOMREF=$(awk -F'\t' '$1 == "homref" {print $2}' "${WORK}/counts.tsv")
N_SITES=$(wc -l < "${PANEL}/panel/sites.tsv" | tr -d ' ')
echo "  ${N_GENO} of the panel's ${N_SITES} SNVs genotyped, ${N_HOMREF} of them as reference (0/0)"
# A VCF without reference calls has only the sites where you differ from the
# reference: fractions from them alone lean toward the groups where those
# alleles are common. Stop rather than report that.
if [ "$INPUT_KIND" = vcf ] && [ "$(( N_HOMREF * 5 ))" -lt "$N_GENO" ]; then
  echo "ERROR: ${VCF} has reference calls (0/0) at only ${N_HOMREF} of the ${N_GENO} panel SNVs it genotypes." >&2
  echo "  Step 26b needs the gVCF beside it (${GVCF}, step 03 writes it) or a VCF that holds reference calls" >&2
  echo "  (chip-to-vcf.sh output does): from variant sites alone the fractions would be biased." >&2
  exit 1
fi
# shellcheck disable=SC2016  # $1 belongs to the inner bash
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" bash -euo pipefail -c \
  'bcftools view -Oz -o "$1.gz" "$1" && bcftools index -t "$1.gz" && rm -f "$1"' _ "${W}/genotypes.vcf"
# The chromosomes to phase and paint: 100 genotyped SNVs or more spanning
# 5 cM or more (bin/ancestry_fine_report.py chroms says why).
CHROMS=$(run_in --cpus 1 --memory 1g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/ancestry_fine_report.py chroms --counts "${W}/counts.tsv" --maps "${P}/maps")

# --- 2. ADMIXTURE: projection on the groups' frequencies ----------------------------
# sample-bed writes the inputs and run_admixture.sh, which projects the sample
# at each level and then runs the bootstrap replicates (ADMIXTURE's own -B
# crashes with -P).
echo "[2/5] ADMIXTURE projection, ${BOOT} bootstrap replicates..."
LEVELS=()
while IFS=$'\t' read -r level k _; do LEVELS+=("--level" "${level}:${k}:${P}/model/${level}.P"); done < "${PANEL}/model/levels.tsv"
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/ancestry_fine_report.py sample-bed --bim "${P}/model/ref.bim" --vcf "${W}/genotypes.vcf.gz" \
    --sample "$SAMPLE" "${LEVELS[@]}" --bootstrap "$BOOT" --out-dir "${W}/admixture"
# shellcheck disable=SC2016  # $1 and $2 belong to the inner sh
run_in --cpus "$THREADS" --memory 4g "$ADMIXTURE_IMAGE" sh -c 'cd "$1" && sh run_admixture.sh "$2"' _ "${W}/admixture" "$THREADS" \
  || { echo "ERROR: ADMIXTURE failed (log above)" >&2; exit 1; }

if [ -z "$CHROMS" ]; then
  echo "[3/5] No chromosome has 100 genotyped SNVs spanning 5 cM: no phasing and no painting."
  echo "[4/5] Skipped."
else
  # --- 3. Beagle: phasing against the panel ------------------------------------------
  echo "[3/5] Phasing with Beagle against the panel (${CHROMS})..."
  # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
  run_in --cpus "$THREADS" --memory 8g "$BEAGLE_IMAGE" bash -euo pipefail -c '
      for c in $4; do
        beagle -Xmx6g gt="$2/genotypes.vcf.gz" chrom="$c" ref="$1/panel/$c.vcf.gz" map="$1/maps/$c.map" \
          out="$2/phase/$c" impute=false seed=1 nthreads="$3" > "$2/phase/$c.out" 2>&1 \
          || { tail -n 20 "$2/phase/$c.out"; exit 1; }
      done' _ "$P" "$W" "$THREADS" "$CHROMS" \
    || { echo "ERROR: Beagle failed (log above)" >&2; exit 1; }

  # --- 4. RFMix: the painting --------------------------------------------------------
  # A chromosome RFMix fails on is named in paint/failed.txt, which the
  # report lists as not painted, and the next one is still painted: one
  # failure no longer costs the painting of every later chromosome.
  echo "[4/5] Painting with RFMix ($(wc -l < "${PANEL}/model/paint_map.tsv" | tr -d ' ') reference individuals)..."
  # shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
  run_in --cpus "$THREADS" --memory 8g "$RFMIX_IMAGE" bash -euo pipefail -c '
      for c in $4; do
        rfmix -f "$2/phase/$c.vcf.gz" -r "$1/panel/$c.vcf.gz" -m "$1/model/paint_map.tsv" -g "$1/maps/$c.rfmix.map" \
          -o "$2/paint/$c" --chromosome="$c" --n-threads="$3" --random-seed=1 > "$2/paint/$c.log" 2>&1 \
          || { tail -n 20 "$2/paint/$c.log"; echo "$c" >> "$2/paint/failed.txt"; rm -f "$2/paint/$c.msp.tsv"; }
      done' _ "$P" "$W" "$THREADS" "$CHROMS" \
    || { echo "ERROR: RFMix failed (log above)" >&2; exit 1; }
  if [ -s "${WORK}/paint/failed.txt" ]; then
    echo "WARNING: RFMix failed on $(paste -sd' ' "${WORK}/paint/failed.txt") (last lines above, logs in ${WORK}/paint); the report lists them as not painted" >&2
  fi
fi

# --- 5. The report ---------------------------------------------------------------------
echo "[5/5] Writing the fractions, the painting and the summary..."
echo ""
run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYTHON_IMAGE" \
  python3 /pgp-bin/ancestry_fine_report.py report --sample "$SAMPLE" --model "${P}/model" \
    --admixture "${W}/admixture" --paint "${W}/paint" --panel "$ANCESTRY_FINE_PANEL" --out "$(cpath "$OUT")"
for f in "${OUT}/${SAMPLE}_ancestry_fine.tsv" "${OUT}/${SAMPLE}_painting.png"; do
  [ -s "$f" ] || { echo "ERROR: step 26b wrote no ${f}" >&2; exit 1; }
done
echo ""
echo "  Input: the ${INPUT_KIND}"
echo "  Fractions: ${OUT}/${SAMPLE}_ancestry_fine.tsv"
echo "  Painting:  ${OUT}/${SAMPLE}_painting.bed and ${SAMPLE}_painting.png"
echo "  Intermediate files (genotypes, phased haplotypes, ADMIXTURE and RFMix output): ${WORK}"
echo "  See docs/26b-ancestry-fine.md for what the names mean and where the resolution ends."
