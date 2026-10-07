#!/usr/bin/env bash
# generate-report.sh — Create a plain-text summary report of all pipeline results
# Usage: ./scripts/generate-report.sh <sample_name>
#
# Reads every step's output once into ${SAMPLE}/summary.json
# (bin/collect_summary.py) and renders the text report from it
# (bin/render_report.py); step 24 renders the HTML report from the same
# summary. A result from an earlier run whose step was skipped or failed in the
# latest run-all.sh run is marked stale with its date. Run after completing all
# (or some) pipeline steps.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
REPORT="${SAMPLE_DIR}/${SAMPLE}_report.txt"

# Check sample directory exists
if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

if [ ! -f "${SAMPLE_DIR}/run_manifest.tsv" ]; then
  GENOME_DIR="$GENOME_DIR" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" generate-report.sh
fi

rm -f "$REPORT"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/render_report.py \
    --sample "$SAMPLE" \
    --sample-dir "/genome/${SAMPLE}" \
    --json "/genome/${SAMPLE}/summary.json" \
    -o "/genome/${SAMPLE}/${SAMPLE}_report.txt"

if [ ! -f "$REPORT" ]; then
  echo "ERROR: the report was not written: ${REPORT}" >&2
  exit 1
fi
cat "$REPORT"
echo ""
echo "Report saved to: ${REPORT}"
