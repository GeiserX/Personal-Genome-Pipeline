#!/usr/bin/env bash
# Imputation Prep — Split VCF by chromosome for Michigan Imputation Server upload
# Creates PASS-only, bgzipped, tabix-indexed VCFs per chromosome
# NOTE: MIS requires 20+ samples per job. Single WGS = mainly useful for phasing.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/imputation"
MIS_DIR="${OUTPUT_DIR}/mis_ready"

echo "=== Imputation Prep: ${SAMPLE} ==="
mkdir -p "$MIS_DIR"

# Step 1: Split by chromosome
for chr in $(seq 1 22); do
  echo "Splitting chr${chr}..."
  docker run --rm --cpus 2 --memory 2g \
    -v "${GENOME_DIR}/${SAMPLE}:/data" \
    staphb/bcftools:1.21 \
    bcftools view -r "chr${chr}" "/data/vcf/${SAMPLE}.vcf.gz" \
      -Oz -o "/data/imputation/${SAMPLE}_chr${chr}.vcf.gz"
  docker run --rm --cpus 1 --memory 1g \
    -v "${GENOME_DIR}/${SAMPLE}:/data" \
    staphb/bcftools:1.21 \
    bcftools index "/data/imputation/${SAMPLE}_chr${chr}.vcf.gz"
done

# Step 2: Create MIS-ready copies (PASS-only + tabix)
# IMPORTANT: Use bcftools -Oz (not bgzip pipe) — bgzip is NOT in staphb/bcftools PATH
for chr in $(seq 1 22); do
  echo "MIS-ready chr${chr}..."
  docker run --rm --cpus 2 --memory 2g \
    -v "${GENOME_DIR}/${SAMPLE}:/data" \
    staphb/bcftools:1.21 bash -c "
      bcftools view -f PASS -Oz -o /data/imputation/mis_ready/${SAMPLE}_chr${chr}.vcf.gz /data/imputation/${SAMPLE}_chr${chr}.vcf.gz
      bcftools index -t /data/imputation/mis_ready/${SAMPLE}_chr${chr}.vcf.gz
    "
done

echo "=== Imputation prep complete ==="
echo "MIS-ready VCFs: ${MIS_DIR}/${SAMPLE}_chr{1-22}.vcf.gz"
echo ""
echo "Next steps:"
echo "  1. Register at https://imputationserver.sph.umich.edu"
echo "  2. Upload all 22 chr VCFs"
echo "  3. Select: refpanel=TOPMed r2, population=eur, build=hg38"
echo "  4. Results expire after 7 days — download immediately"
