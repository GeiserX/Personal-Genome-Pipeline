#!/usr/bin/env bash
# MultiQC — aggregate QC reports into a single HTML dashboard
# Input: all QC outputs from previous steps (fastp, mosdepth, samtools, etc.)
# Output: single HTML report in $GENOME_DIR/<sample>/multiqc/
#
# MultiQC auto-discovers supported tool outputs by scanning the sample directory.
# Supported tools in this pipeline: fastp (JSON), mosdepth, samtools flagstat/stats.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT_DIR="${SAMPLE_DIR}/multiqc"

echo "=== MultiQC: ${SAMPLE} ==="
echo "Scanning: ${SAMPLE_DIR}/"

if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

# No skip: the report takes seconds and reads every QC file, so it is rebuilt
# on every run and never shows an older state of the sample.
mkdir -p "$OUTPUT_DIR"

# Generate samtools flagstat if BAM exists and flagstat doesn't
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
FLAGSTAT="${SAMPLE_DIR}/aligned/${SAMPLE}_flagstat.txt"
# Written through a temporary name: a failed flagstat leaves no empty file
# that every later run would take as done, and does not stop the report.
if [ -f "$BAM" ] && [ ! -s "$FLAGSTAT" ]; then
  echo "Generating samtools flagstat for MultiQC..."
  if ! atomic_out "$FLAGSTAT" run_in --cpus 2 --memory 2g \
      "${SAMTOOLS_IMAGE}" \
      samtools flagstat "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"; then
    echo "WARNING: samtools flagstat failed; the report has no flagstat section."
  fi
fi

# Run MultiQC
# Flags:
#   -f            Force overwrite existing reports
#   -o            Output directory
#   -n            Report filename
#   --title       Report title shown in HTML
#   --no-version-check  Do not ask the MultiQC server for a newer release
#   --no-ai       No AI summary (it would send report data to an outside service)
echo "Running MultiQC..."
run_in --cpus 2 --memory 2g \
  "${MULTIQC_IMAGE}" \
  multiqc \
    "/genome/${SAMPLE}" \
    -f \
    --no-version-check \
    --no-ai \
    -o "/genome/${SAMPLE}/multiqc" \
    -n "multiqc_report.html" \
    --title "${SAMPLE} — Personal Genome Pipeline QC"

echo "=== MultiQC complete ==="
echo "Report: ${OUTPUT_DIR}/multiqc_report.html"
echo "Open in browser to view aggregated QC dashboard."
