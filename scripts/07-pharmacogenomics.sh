#!/usr/bin/env bash
# PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
# Input: VCF.gz + GRCh38 reference
# Output: HTML + JSON reports with metabolizer status for 23 pharmacogenes
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/vcf"

echo "=== PharmCAT: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "Outputs: ${OUTPUT_DIR}/${SAMPLE}.report.html and ${OUTPUT_DIR}/${SAMPLE}.report.json"

# Validate inputs (PharmCAT requires both VCF and its tabix index)
for f in "$VCF" "${VCF}.tbi" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    if [ "$f" = "${VCF}.tbi" ]; then
      echo "  PharmCAT requires a tabix index. Generate it with:" >&2
      echo "  docker run --rm -v \"\${GENOME_DIR}:/genome\" staphb/bcftools:1.21 bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" >&2
    fi
    exit 1
  fi
done

# Step 1: Preprocess VCF (normalize, filter to PGx positions)
docker run --rm \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "${GENOME_DIR}/reference:/ref" \
  "${PHARMCAT_IMAGE}" \
  python3 /pharmcat/pharmcat_vcf_preprocessor \
    -vcf "/data/${SAMPLE}.vcf.gz" \
    -refFna /ref/Homo_sapiens_assembly38.fasta \
    -o /data/ \
    -bf "$SAMPLE"

# Step 2: Run PharmCAT on preprocessed VCF
docker run --rm \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${PHARMCAT_IMAGE}" \
  java -jar /pharmcat/pharmcat.jar \
    -vcf "/data/${SAMPLE}.preprocessed.vcf.bgz" \
    -o /data/ \
    -bf "$SAMPLE" \
    -reporterJson \
    -reporterHtml

echo "=== PharmCAT complete ==="
echo "Reports: ${OUTPUT_DIR}/${SAMPLE}.report.html and ${OUTPUT_DIR}/${SAMPLE}.report.json"
echo ""
echo "Key genes covered: CYP2C19, CYP2D6, CYP2B6, CYP3A5, UGT1A1, DPYD, NAT2, TPMT"
echo "NOTE: CYP2D6 may return 'Not called' — use Cyrius (BAM-based) for CYP2D6."
