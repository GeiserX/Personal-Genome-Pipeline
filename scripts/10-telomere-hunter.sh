#!/usr/bin/env bash
# TelomereHunter — Estimate telomere length from WGS BAM
# Output: tel_content metric (GC-corrected telomeric reads per million)
# Higher values = longer telomeres. Provides biological age baseline.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}
SAMPLE_DIR="${GENOMA_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/telomere/${SAMPLE}"

echo "=== TelomereHunter: ${SAMPLE} ==="
echo "WARNING: This reads the entire BAM (~30-40GB). Takes 30-60 minutes."
mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 4 --memory 4g \
  --user root \
  -v "${SAMPLE_DIR}/aligned:/bam" \
  -v "${OUTPUT_DIR}:/output" \
  lgalarno/telomerehunter:latest \
  telomerehunter -ibt "/bam/${SAMPLE}_sorted.bam" -o /output -p "$SAMPLE"

echo "=== TelomereHunter complete ==="
SUMMARY="${OUTPUT_DIR}/tumor_TelomerCnt_${SAMPLE}/${SAMPLE}_tumor_summary.tsv"
if [ -f "$SUMMARY" ]; then
  echo "Summary: $SUMMARY"
  TEL_CONTENT=$(awk -F'\t' 'NR==2 {print $NF}' "$SUMMARY")
  echo "Telomere content: ${TEL_CONTENT}"
fi
