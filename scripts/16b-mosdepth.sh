#!/usr/bin/env bash
# mosdepth — fast per-base and per-region coverage statistics from BAM
# Input: sorted BAM + BAI from $GENOME_DIR/<sample>/aligned/
# Output: coverage distributions, thresholds, summary in $GENOME_DIR/<sample>/mosdepth/
#
# Complements indexcov (step 16): indexcov reads only the BAM index for a quick
# estimate; mosdepth reads actual alignments for precise per-base depth.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/mosdepth"

echo "=== mosdepth Coverage: ${SAMPLE} ==="
echo "BAM: ${BAM}"

# Validate inputs
if [ ! -f "$BAM" ]; then
  echo "ERROR: BAM not found: ${BAM}" >&2
  exit 1
fi
if [ ! -f "${BAM}.bai" ]; then
  echo "ERROR: BAM index not found: ${BAM}.bai" >&2
  echo "Create it with: samtools index ${BAM}" >&2
  exit 1
fi

# Skip if output already exists
if [ -f "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt" ]; then
  echo "mosdepth output already exists in ${OUTPUT_DIR}/, skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

mkdir -p "$OUTPUT_DIR"

# mosdepth flags:
#   --by 500        Window size for per-region coverage (500bp bins)
#   --fast-mode     Skip per-base output (saves disk + time), keep distributions
#   --threads       Decompression threads (mosdepth uses 1 main + N decompression)
#   --thresholds    Report fraction of bases at these coverage thresholds
#   --no-per-base   Omit the large per-base BED (redundant with --fast-mode)
echo "Computing coverage statistics..."
docker run --rm --user root \
  --cpus "${THREADS}" --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0 \
  mosdepth \
    --by 500 \
    --fast-mode \
    --threads "${THREADS}" \
    --thresholds 1,5,10,15,20,30,50 \
    "/genome/${SAMPLE}/mosdepth/${SAMPLE}" \
    "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"

echo "=== mosdepth complete ==="
echo "Summary:        ${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
echo "Distribution:   ${OUTPUT_DIR}/${SAMPLE}.mosdepth.global.dist.txt"
echo "Region coverage: ${OUTPUT_DIR}/${SAMPLE}.regions.bed.gz"
echo "Thresholds:     ${OUTPUT_DIR}/${SAMPLE}.thresholds.bed.gz"

# Print summary
if [ -f "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt" ]; then
  echo ""
  echo "Coverage summary:"
  head -5 "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
fi
