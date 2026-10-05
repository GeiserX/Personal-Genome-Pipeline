#!/usr/bin/env bash
# 21-cyrius.sh — [OPT-IN] CYP2D6 star allele calling using Cyrius
# Usage: ./scripts/21-cyrius.sh <sample_name>
#
# CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
# and complex structural variants (gene deletions, duplications, hybrids).
# Cyrius uses depth-based analysis specifically designed for CYP2D6.
#
# OPT-IN: a default run leaves this step out (run-all.sh runs it with
# TOOLS=...,cyrius). Cyrius 1.1.1 (no release since 2021) is under the
# PolyForm Strict licence 1.0.0: non-commercial use only. `setup.sh --cyrius`
# installs it once, from the hash-locked scripts/cyrius-constraints.txt; this
# step runs it with no network. It is the second CYP2D6 caller step 36 needs
# before a CYP2D6 call reaches PharmCAT.
#
# Before Cyrius runs, mosdepth measures the depth over CYP2D6 and its flanks
# (bin/cyp2d6_depth_check.py). When the reads there are multi-mapped (a BAM
# aligned to a reference with ALT contigs), the call is marked: its Filter
# becomes CYP2D6_depth_unreliable and step 36 does not pass it on.
#
# Requires: Sorted BAM with index
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
BAI="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai"
OUTDIR="${GENOME_DIR}/${SAMPLE}/cyrius"
CYRIUS_DIR="${GENOME_DIR}/tools/cyrius-${CYRIUS_VERSION}"
mkdir -p "$OUTDIR"
# A run that stops anywhere below must not leave the last run's call behind.
RESULT_FILE="${OUTDIR}/${SAMPLE}_cyp2d6.tsv"
rm -f "$RESULT_FILE"

# This step feeds step 36's outside calls for PharmCAT. Remove the ones made
# from an earlier result, so step 07 never reads a call this run has not
# confirmed; step 36 writes them again.
rm -f "${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_outside_calls.tsv" \
  "${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv"

# Validate inputs
for FILE in "$BAM" "$BAI"; do
  if [ ! -f "$FILE" ]; then
    echo "ERROR: Required file not found: ${FILE}"
    exit 1
  fi
done
# The install must be the one setup.sh makes from today's lock file and image.
STAMP="python=${PYTHON_IMAGE} lock=$(_digest sha256 "${PGP_ROOT}/scripts/cyrius-constraints.txt")"
if [ "$(cat "${CYRIUS_DIR}/INSTALLED" 2>/dev/null)" != "$STAMP" ]; then
  echo "ERROR: Cyrius ${CYRIUS_VERSION} is not installed for this version of the pipeline (${CYRIUS_DIR})." >&2
  echo "  Install it once (it downloads from PyPI; Cyrius is for non-commercial use only):" >&2
  echo "    ./scripts/setup.sh --cyrius ${GENOME_DIR}" >&2
  exit 1
fi

echo "============================================"
echo "  Step 21: CYP2D6 Star Allele Calling"
echo "  Tool: Cyrius ${CYRIUS_VERSION} (Illumina), opt-in"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${BAM}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# [1/3] Depth over CYP2D6 and its flanks, all reads and MAPQ >= 1. The regions
# come from bin/cyp2d6_depth_check.py, which also judges them.
echo "[1/3] CYP2D6 depth check..."
DEPTH_DIR="${OUTDIR}/cyp2d6_depth"
CHECK="${OUTDIR}/${SAMPLE}_cyp2d6_depth_check.tsv"
mkdir -p "$DEPTH_DIR"
rm -f "$CHECK"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py bed > "${DEPTH_DIR}/regions.bed"
for Q in 0 1; do
  run_in --cpus 2 --memory 2g "${MOSDEPTH_IMAGE}" \
    mosdepth -n -c chr22 -t 2 -Q "$Q" -b "$(cpath "${DEPTH_DIR}/regions.bed")" \
      "$(cpath "${DEPTH_DIR}/q${Q}")" "$(cpath "$BAM")"
done
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py check \
    --all "$(cpath "${DEPTH_DIR}/q0.regions.bed.gz")" \
    --mapq1 "$(cpath "${DEPTH_DIR}/q1.regions.bed.gz")" \
    --out "$(cpath "$CHECK")"
# A check that wrote nothing counts as failed: the call is then marked.
DEPTH_STATUS=$(awk -F'\t' '$1 == "status" {print $2}' "$CHECK" 2>/dev/null || true)

# [2/3] Cyrius, from the install setup.sh made, with no network. The manifest
# (the BAM path) is created inside the container.
echo "[2/3] Running Cyrius CYP2D6 caller..."
# shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
run_in --cpus 4 --memory 8g -w /tmp \
  "${PYTHON_IMAGE}" \
  bash -euo pipefail -c '
    echo "$2" > /tmp/manifest.txt
    PYTHONPATH="$1" python3 -m cyrius \
      --manifest /tmp/manifest.txt \
      --genome 38 \
      --prefix "$3" \
      --outDir "$4" \
      --threads 4' \
  _ "$(cpath "$CYRIUS_DIR")" "$(cpath "$BAM")" "${SAMPLE}_cyp2d6" "$(cpath "$OUTDIR")/"

echo ""
echo "[3/3] Parsing results..."

if [ ! -f "$RESULT_FILE" ]; then
  echo "ERROR: Cyrius finished but wrote no ${RESULT_FILE}. Check the messages above." >&2
  exit 1
fi
if [ "$DEPTH_STATUS" != ok ]; then
  # Keep Cyrius's genotype for the record; the Filter says it cannot be used.
  awk -F'\t' -v OFS='\t' 'NR > 1 {$3 = "CYP2D6_depth_unreliable"} {print}' "$RESULT_FILE" > "${RESULT_FILE}.tmp"
  mv "${RESULT_FILE}.tmp" "$RESULT_FILE"
  MSG=$(awk -F'\t' '$1 == "message" {print $2}' "$CHECK" 2>/dev/null || true)
  echo "WARNING: ${MSG:-the CYP2D6 depth check wrote no result}"
  echo "  The call below is marked CYP2D6_depth_unreliable and step 36 does not pass it to PharmCAT."
fi
echo ""
echo "  CYP2D6 Results:"
echo "  ─────────────────"
column -t "$RESULT_FILE" 2>/dev/null || cat "$RESULT_FILE"
echo ""
DIPLOTYPE=$(awk -F'\t' 'NR==2 {print $2}' "$RESULT_FILE" 2>/dev/null || echo "N/A")
echo "  Diplotype: ${DIPLOTYPE}"
echo ""
echo "  Step 36 compares it with pypgx (step 32); only an agreed call reaches PharmCAT."
echo "  Interpret this with PharmGKB: https://www.pharmgkb.org/gene/PA128"

echo ""
echo "============================================"
echo "  CYP2D6 calling complete: ${SAMPLE}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
