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
#
# Then it writes what changed since the previous run of the sample
# (bin/run_delta.py, docs/delta-report.md): ${SAMPLE}_delta.txt, .html and
# .json, against the newest copy under history/, and keeps this run's
# summary.json and run_manifest.tsv as history/<UTC time>/ for the next run.
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

# What changed since the previous run, then keep this run as the next baseline.
HISTORY="${SAMPLE_DIR}/history"
PREV=$(find "$HISTORY" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | LC_ALL=C sort | tail -n 1 || true)
PREV_ARG=()
[ -z "$PREV" ] || PREV_ARG=(--previous "/genome/${SAMPLE}/history/$(basename "$PREV")")
echo ""
if ! run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/run_delta.py \
    --sample "$SAMPLE" \
    --summary "/genome/${SAMPLE}/summary.json" \
    --manifest "/genome/${SAMPLE}/run_manifest.tsv" \
    ${PREV_ARG[@]+"${PREV_ARG[@]}"} \
    --out "/genome/${SAMPLE}/${SAMPLE}_delta"; then
  echo "ERROR: the delta against the previous run was not written (bin/run_delta.py failed)" >&2
  exit 1
fi
# history/<UTC time>[-NN]/: sorts by name in run order.
if [ -f "${SAMPLE_DIR}/summary.json" ]; then
  stamp=$(date -u '+%Y-%m-%dT%H%M%SZ') dest="${HISTORY}/${stamp}" n=1
  while [ -e "$dest" ]; do n=$((n + 1)); dest=$(printf '%s/%s-%02d' "$HISTORY" "$stamp" "$n"); done
  mkdir -p "$dest"
  cp -p "${SAMPLE_DIR}/summary.json" "$dest/"
  [ ! -f "${SAMPLE_DIR}/run_manifest.tsv" ] || cp -p "${SAMPLE_DIR}/run_manifest.tsv" "$dest/"
  echo "Delta saved to: ${SAMPLE_DIR}/${SAMPLE}_delta.txt (this run kept as ${dest#"${SAMPLE_DIR}"/})"
else
  echo "WARNING: no ${SAMPLE_DIR}/summary.json, so this run is not kept in history/" >&2
fi
