#!/usr/bin/env bash
# ClinVar Pathogenic Screen — intersect sample VCF with ClinVar pathogenic/LP variants
# Finds known disease-causing variants the person carries
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}
VCF="${GENOMA_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
CLINVAR="${GENOMA_DIR}/reference/clinvar_pathogenic_chr.vcf.gz"
OUTPUT_DIR="${GENOMA_DIR}/${SAMPLE}/clinvar"

echo "=== ClinVar Pathogenic Screen: ${SAMPLE} ==="
mkdir -p "$OUTPUT_DIR"

# Step 1: Extract PASS variants only
docker run --rm --cpus 2 --memory 2g \
  -v "${GENOMA_DIR}:/genoma" \
  staphb/bcftools:1.21 \
  bcftools view -f PASS "/genoma/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    -Oz -o "/genoma/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz"

docker run --rm --cpus 1 --memory 1g \
  -v "${GENOMA_DIR}:/genoma" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genoma/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz"

# Step 2: Intersect with ClinVar pathogenic
docker run --rm --cpus 2 --memory 2g \
  -v "${GENOMA_DIR}:/genoma" \
  staphb/bcftools:1.21 \
  bcftools isec -p "/genoma/${SAMPLE}/clinvar/isec" \
    "/genoma/${SAMPLE}/clinvar/${SAMPLE}_pass.vcf.gz" \
    /genoma/reference/clinvar_pathogenic_chr.vcf.gz

echo "=== ClinVar screen complete ==="
echo "Shared variants: ${OUTPUT_DIR}/isec/0002.vcf (in both sample AND ClinVar)"
echo "Count: $(grep -c -v '^#' "${OUTPUT_DIR}/isec/0002.vcf" 2>/dev/null || echo 0) pathogenic hits"
