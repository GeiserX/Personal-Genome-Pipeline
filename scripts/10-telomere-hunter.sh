#!/usr/bin/env bash
# TelomereHunter — Estimate telomere length from WGS BAM
# Output: tel_content metric (GC-corrected telomeric reads per million)
# Higher values = longer telomeres. Provides biological age baseline.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/telomere/${SAMPLE}"

echo "=== TelomereHunter: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"
echo "WARNING: This reads the entire BAM (~30-40GB). Takes 30-60 minutes."

for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 4 --memory 4g \
  --user root \
  -v "${GENOME_DIR}:/genome" \
  lgalarno/telomerehunter:latest \  # No versioned tags from publisher
  telomerehunter \
    -ibt "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    -o "/genome/${SAMPLE}/telomere/${SAMPLE}" \
    -p "$SAMPLE"

echo "=== TelomereHunter complete ==="
SUMMARY="${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_summary.tsv"
if [ -f "$SUMMARY" ]; then
  echo "Summary: $SUMMARY"
  TEL_CONTENT=$(awk -F'\t' 'NR==2 {print $11}' "$SUMMARY")
  echo "Telomere content: ${TEL_CONTENT}"
fi
