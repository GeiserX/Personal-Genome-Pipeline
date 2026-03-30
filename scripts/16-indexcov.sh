#!/usr/bin/env bash
# indexcov (goleft) — Rapid coverage QC and sex chromosome check from BAM index
# Input: sorted BAM with .bai index
# Output: HTML report with per-chromosome coverage + sex check
# Instant (~5 seconds, reads only the BAM index)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/indexcov"

echo "=== indexcov: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/goleft:0.2.4--h9ee0642_1 \
  goleft indexcov \
    --directory "/genome/${SAMPLE}/indexcov" \
    "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

echo "=== indexcov complete ==="
echo "Results: ${OUTPUT_DIR}/"
echo ""
echo "Key outputs:"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.html  — interactive coverage plot"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.ped   — sex chromosome inference"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.roc   — coverage uniformity"
echo ""
if [ -f "${OUTPUT_DIR}/indexcov-indexcov.ped" ]; then
  echo "Sex check (from X/Y coverage ratio):"
  tail -1 "${OUTPUT_DIR}/indexcov-indexcov.ped" | awk '{print "  Predicted sex: " ($6==1 ? "male" : "female")}'
fi
