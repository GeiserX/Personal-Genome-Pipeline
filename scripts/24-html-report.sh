#!/usr/bin/env bash
# 24-html-report.sh — Generate a self-contained HTML report of all pipeline results
# Usage: ./scripts/24-html-report.sh <sample_name>
#
# Reads every step's output once into ${SAMPLE}/summary.json
# (bin/collect_summary.py) and renders the HTML report from it
# (bin/render_report.py). generate-report.sh renders the text report from the
# same summary, so the two reports show the same values. A result from an
# earlier run whose step was skipped or failed in the latest run-all.sh run is
# marked stale with its date. No internet connection is needed to view it.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT="${SAMPLE_DIR}/${SAMPLE}_report.html"

if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

echo "============================================"
echo "  Step 24: HTML Report Generator"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""

# What produced these outputs: run-all.sh writes the manifest when a run
# starts; a sample without one gets it now.
if [ ! -f "${SAMPLE_DIR}/run_manifest.tsv" ]; then
  GENOME_DIR="$GENOME_DIR" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" 24-html-report.sh
fi

echo "Reading pipeline outputs and writing the report..."
rm -f "$OUTPUT"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/render_report.py \
    --sample "$SAMPLE" \
    --sample-dir "/genome/${SAMPLE}" \
    --json "/genome/${SAMPLE}/summary.json" \
    -o "/genome/${SAMPLE}/${SAMPLE}_report.html"

if [ ! -f "$OUTPUT" ]; then
  echo "ERROR: the report was not written: ${OUTPUT}" >&2
  exit 1
fi
REPORT_KB=$(( $(wc -c < "$OUTPUT") / 1024 ))

echo ""
echo "============================================"
echo "  HTML report generated: ${OUTPUT}"
echo "  Summary: ${SAMPLE_DIR}/summary.json"
echo "  Size: ${REPORT_KB} KB"
echo "============================================"
echo ""
echo "Open in your browser:"
echo "  open ${OUTPUT}          # macOS"
echo "  xdg-open ${OUTPUT}     # Linux"
echo "  start ${OUTPUT}         # Windows (WSL)"
