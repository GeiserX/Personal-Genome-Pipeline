#!/usr/bin/env bash
# Mitochondrial Haplogroup — Determine maternal lineage from mtDNA variants
# Input: VCF.gz (full genome — chrM will be extracted)
# Output: haplogroup classification with quality score
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/mito"

echo "=== Mitochondrial Haplogroup: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Step 1: Extract chrM variants
echo "Extracting chrM variants..."
docker run --rm \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/genome/${SAMPLE}/vcf" \
  staphb/bcftools:1.21 \
  bcftools view -r chrM "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz"

docker run --rm \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/genome/${SAMPLE}/vcf" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz"

# Step 2: Run haplogrep3
echo "Classifying haplogroup..."
source "$(dirname "$0")/../versions.env" 2>/dev/null || HAPLOGREP3_IMAGE="jtb114/haplogrep3:latest"
docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}/${SAMPLE}:/genome/${SAMPLE}" \
  ${HAPLOGREP3_IMAGE} \
  classify \
    --tree phylotree-fu-rcrs@1.2 \
    --input "/genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz" \
    --output "/genome/${SAMPLE}/mito/${SAMPLE}_haplogroup.txt" \
    --extend-report

echo "=== Haplogrep3 complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt"
if [ -f "${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt" ]; then
  echo ""
  echo "Haplogroup:"
  head -5 "${OUTPUT_DIR}/${SAMPLE}_haplogroup.txt"
fi
echo ""
echo "Quality > 0.9 = high confidence. Common European: H, U, J, T, K, V, W, X"
