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
#   --by 500        Window size for per-region coverage (500bp bins; WGS default)
#   --by <BED>      Capture BED for on-target coverage (WES; set CAPTURE_BED env var)
#   --fast-mode     Use faster algorithm (does NOT skip per-base output by itself)
#   --no-per-base   Omit the large per-base BED (saves ~3-5 GB per WGS sample)
#   --threads       Decompression threads (mosdepth uses 1 main + N decompression)
#   --thresholds    Report fraction of bases at these coverage thresholds

# Use capture BED for WES on-target coverage, or 500bp bins for WGS
if [ -n "${CAPTURE_BED:-}" ]; then
  if [ ! -f "${CAPTURE_BED}" ]; then
    echo "ERROR: CAPTURE_BED not found: ${CAPTURE_BED}" >&2
    exit 1
  fi
  # Compute container-relative path for the BED file
  CAPTURE_BED_REL="${CAPTURE_BED#"${GENOME_DIR}/"}"
  BY_FLAG="/genome/${CAPTURE_BED_REL}"
  echo "Using capture BED for WES on-target coverage: ${CAPTURE_BED}"
else
  BY_FLAG="500"
fi

echo "Computing coverage statistics..."
docker run --rm --user root \
  --cpus "${THREADS}" --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0 \
  mosdepth \
    --by "${BY_FLAG}" \
    --fast-mode \
    --no-per-base \
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
