#!/usr/bin/env bash
# 26-ancestry.sh — Genetic ancestry: your place among a reference panel
# Usage: ./scripts/26-ancestry.sh <sample_name>
#
# One sample cannot be placed by a PCA of its own. This step projects it onto
# the principal components of a reference panel of known samples instead,
# with pgsc_calc's ancestry projection (the same run that gives step 25 its
# percentiles), and writes the principal components and the panel population
# the sample is most similar to: ancestry/<sample>_ancestry.tsv.
#
# The panel is pgsc_calc's 1000 Genomes database (${PGSC_PANEL} in
# versions.env, about 7 GB), installed by scripts/setup.sh --ancestry-panel.
# Without it the step says so in one line and exits 0.
#
# Requires: the VCF of step 3 (its gVCF too, for a projection that counts the
# sites where you match the reference), Java 17+ and Nextflow. Runs step 25,
# which uses the panel when it is installed. Env: ANCESTRY_PANEL as in step 25.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

PANEL=${ANCESTRY_PANEL:-${GENOME_DIR}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst}
OUT="${GENOME_DIR}/${SAMPLE}/ancestry/${SAMPLE}_ancestry.tsv"

if [ "${ANCESTRY_PANEL:-}" = none ] || [ ! -f "$PANEL" ]; then
  echo "Step 26 skipped: no ancestry reference panel at ${PANEL}; install it with scripts/setup.sh --ancestry-panel ${GENOME_DIR}"
  exit 0
fi

echo "============================================"
echo "  Step 26: Ancestry (projection onto ${PANEL##*/})"
echo "  Tool: pgsc_calc ${PGSC_CALC_VERSION}, through step 25"
echo "  Sample: ${SAMPLE}"
echo "  Output: ${OUT}"
echo "============================================"
echo ""

rm -f "$OUT"
"$(dirname "$0")/25-prs.sh" "$SAMPLE"

if [ ! -s "$OUT" ]; then
  echo "ERROR: step 25 ran with the panel but wrote no ${OUT}; see ${GENOME_DIR}/${SAMPLE}/prs/pgsc_calc/pgsc_calc.log" >&2
  exit 1
fi
echo ""
echo "  Population most similar to ${SAMPLE}: $(awk -F'\t' '$1 == "population" {print $2}' "$OUT")"
echo "  Principal components and the probability of each panel population: ${OUT}"
echo "  See docs/26-ancestry.md for what the label means and what it does not."
