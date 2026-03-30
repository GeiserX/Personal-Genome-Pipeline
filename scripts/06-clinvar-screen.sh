#!/usr/bin/env bash
# ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
# Finds known disease-causing variants the person carries
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
CLINVAR="${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/clinvar"

echo "=== ClinVar Pathogenic Screen: ${SAMPLE} ==="

for f in "$VCF" "${VCF}.tbi" "$CLINVAR" "${CLINVAR}.tbi"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Step 1: Extract PASS variants only
docker run --rm --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -f PASS "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz"

docker run --rm --cpus 1 --memory 1g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz"

# Step 2: Intersect with ClinVar pathogenic
docker run --rm --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools isec -p "/genome/${SAMPLE}/clinvar/isec" \
    "/genome/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz" \
    /genome/clinvar/clinvar_pathogenic_chr.vcf.gz

echo "=== ClinVar screen complete ==="
echo "Shared variants: ${OUTPUT_DIR}/isec/0002.vcf (in both sample AND ClinVar)"
echo "Count: $(grep -c -v '^#' "${OUTPUT_DIR}/isec/0002.vcf" 2>/dev/null || echo 0) pathogenic hits"
