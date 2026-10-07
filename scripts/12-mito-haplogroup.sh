#!/usr/bin/env bash
# Mitochondrial Haplogroup — Determine maternal lineage from mtDNA variants,
# and check the mtDNA for a second person's reads (contamination)
# Input: step 20's Mutect2 chrM calls (mito/<sample>_chrM_filtered.vcf.gz)
#        when they exist, else the chrM records of the step 03 VCF
# Output: mito/<sample>_haplogroup.txt (haplogrep3), and with the Mutect2
#         calls mito/<sample>_haplocheck.txt (haplocheck)
#
# Mutect2 in mitochondrial mode is the caller made for chrM: it reports each
# allele's fraction, so haplogrep3 sees heteroplasmies and haplocheck can look
# for two haplogroups mixed in one sample. DeepVariant calls chrM as a diploid
# nuclear contig, so its calls are the fallback, without the contamination check.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image HAPLOGREP3_IMAGE HAPLOCHECK_IMAGE BCFTOOLS_IMAGE
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
MUTECT2="${GENOME_DIR}/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/mito"
CHRM="${OUTPUT_DIR}/${SAMPLE}_chrM_for_haplogroup.vcf.gz"
CHECK="${OUTPUT_DIR}/${SAMPLE}_haplocheck.txt"

echo "=== Mitochondrial Haplogroup: ${SAMPLE} ==="
mkdir -p "$OUTPUT_DIR"
rm -f "$CHECK" "${CHECK%.txt}.html"

# Step 1: the chrM calls haplogrep3 reads
if [ -f "$MUTECT2" ]; then
  SOURCE=Mutect2
  echo "Input: Mutect2 chrM calls of step 20 (${MUTECT2})"
  # PASS records, one allele per record.
  run_in --cpus 1 --memory 1g "${BCFTOOLS_IMAGE}" bash -euo pipefail -c \
    'bcftools view -f PASS "$1" | bcftools norm -m-any -Oz -o "$2" && bcftools index -f -t "$2"' \
    _ "$(cpath "$MUTECT2")" "$(cpath "$CHRM")"
elif [ -f "$VCF" ]; then
  SOURCE=DeepVariant
  echo "Input: chrM of the step 03 VCF (${VCF}); run scripts/20-mtoolbox.sh first for Mutect2's chrM calls and the contamination check"
  run_in --cpus 1 --memory 1g "${BCFTOOLS_IMAGE}" bash -euo pipefail -c \
    'bcftools view -r chrM "$1" -Oz -o "$2" && bcftools index -f -t "$2"' \
    _ "$(cpath "$VCF")" "$(cpath "$CHRM")"
else
  echo "ERROR: no chrM calls: neither ${MUTECT2} (step 20) nor ${VCF} (step 03)" >&2
  exit 1
fi

# Step 2: Run haplogrep3
# The image has no entrypoint; haplogrep3 is on PATH.
echo "Classifying haplogroup..."
run_in \
  --cpus 2 --memory 2g \
  "${HAPLOGREP3_IMAGE}" \
  haplogrep3 classify \
    --tree phylotree-fu-rcrs@1.2 \
    --input "$(cpath "$CHRM")" \
    --output "/genome/${SAMPLE}/mito/${SAMPLE}_haplogroup.txt" \
    --extend-report

# Step 3: contamination, from the allele fractions Mutect2 reports
if [ "$SOURCE" = Mutect2 ]; then
  echo "Checking for contamination (haplocheck)..."
  run_in --cpus 1 --memory 2g "${HAPLOCHECK_IMAGE}" \
    haplocheck --out "$(cpath "$CHECK")" "$(cpath "$CHRM")"
  [ -s "$CHECK" ] || { echo "ERROR: haplocheck wrote no report (${CHECK})" >&2; exit 1; }
fi

echo "=== Haplogrep3 complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt (from ${SOURCE} calls)"
if [ -f "${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt" ]; then
  echo ""
  echo "Haplogroup:"
  head -5 "${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt"
fi
if [ -s "$CHECK" ]; then
  # Columns by header name; haplocheck quotes its values.
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/"/, "", h); c[h] = i }; next }
    NR == 2 { s = $c["Contamination Status"]; l = $c["Contamination Level"]; gsub(/"/, "", s); gsub(/"/, "", l)
              print "Contamination (haplocheck): " s " (level " l ")" }' "$CHECK"
  echo "  YES: two mtDNA haplogroups are mixed in the reads, so another person's DNA is likely in the sample."
  echo "  Report: ${CHECK}"
else
  echo "Contamination (haplocheck): not checked (needs step 20's Mutect2 calls)"
fi
echo ""
echo "Quality > 0.9 = high confidence. Common European: H, U, J, T, K, V, W, X"
