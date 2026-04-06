#!/usr/bin/env bash
# MultiQC — aggregate QC reports into a single HTML dashboard
# Input: all QC outputs from previous steps (fastp, mosdepth, samtools, etc.)
# Output: single HTML report in $GENOME_DIR/<sample>/multiqc/
#
# MultiQC auto-discovers supported tool outputs by scanning the sample directory.
# Supported tools in this pipeline: fastp (JSON), mosdepth, samtools flagstat/stats.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT_DIR="${SAMPLE_DIR}/multiqc"

echo "=== MultiQC: ${SAMPLE} ==="
echo "Scanning: ${SAMPLE_DIR}/"

if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

# Skip if output already exists
if [ -f "${OUTPUT_DIR}/multiqc_report.html" ]; then
  echo "MultiQC report already exists, skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

mkdir -p "$OUTPUT_DIR"

# Generate samtools flagstat if BAM exists and flagstat doesn't
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
FLAGSTAT="${SAMPLE_DIR}/aligned/${SAMPLE}_flagstat.txt"
if [ -f "$BAM" ] && [ ! -f "$FLAGSTAT" ]; then
  echo "Generating samtools flagstat for MultiQC..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${SAMTOOLS_IMAGE}" \
    samtools flagstat "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    > "$FLAGSTAT" 2>/dev/null || true
fi

# Run MultiQC
# Flags:
#   -f            Force overwrite existing reports
#   -o            Output directory
#   -n            Report filename
#   --title       Report title shown in HTML
#   --no-data-dir Skip creating multiqc_data/ directory (just the HTML)
echo "Running MultiQC..."
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${MULTIQC_IMAGE}" \
  multiqc \
    "/genome/${SAMPLE}" \
    -f \
    -o "/genome/${SAMPLE}/multiqc" \
    -n "multiqc_report.html" \
    --title "${SAMPLE} — Personal Genome Pipeline QC"

echo "=== MultiQC complete ==="
echo "Report: ${OUTPUT_DIR}/multiqc_report.html"
echo "Open in browser to view aggregated QC dashboard."
