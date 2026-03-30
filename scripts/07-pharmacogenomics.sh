#!/usr/bin/env bash
# PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
# Input: VCF.gz + GRCh38 reference
# Output: HTML report with metabolizer status for 23 pharmacogenes
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/vcf"

echo "=== PharmCAT: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.report.html"

# Validate inputs
for f in "$VCF" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

docker run --rm \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "${GENOME_DIR}/reference:/ref" \
  pgkb/pharmcat:2.15.5 \
  java -jar /pharmcat/pharmcat.jar \
    -vcf "/data/${SAMPLE}.vcf.gz" \
    -refFasta /ref/Homo_sapiens_assembly38.fasta \
    -o /data/ \
    -bf "$SAMPLE"

echo "=== PharmCAT complete ==="
echo "Report: ${OUTPUT_DIR}/${SAMPLE}.report.html"
echo ""
echo "Key genes covered: CYP2C19, CYP2D6, CYP2B6, CYP3A5, UGT1A1, DPYD, NAT2, TPMT"
echo "NOTE: CYP2D6 may return 'Not called' — use Cyrius (BAM-based) for CYP2D6."
