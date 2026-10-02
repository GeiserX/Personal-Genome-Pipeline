#!/usr/bin/env bash
# 27-cpic-lookup.sh — Look up CPIC drug-gene recommendations from PharmCAT results
# Usage: ./scripts/27-cpic-lookup.sh <sample_name>
#
# Reads PharmCAT's report.json (step 7) with bin/pgx_parse.py and writes the
# gene phenotypes and the medications PharmCAT's own report matches to them,
# so a gene PharmCAT calls is never dropped for missing from a hand-kept
# table. With pypgx output (step 32) it also writes the PharmCAT/pypgx
# comparison and warns about genes only pypgx could call.
#
# Requires: PharmCAT output from step 7. No internet connection needed.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
# Find PharmCAT JSON output
PHARMCAT_JSON=""
for DIR in "${GENOME_DIR}/${SAMPLE}/pharmcat" "${GENOME_DIR}/${SAMPLE}/vcf"; do
  for FILE in "${DIR}"/*.report.json "${DIR}"/*_pharmcat.json; do
    if [ -f "$FILE" ] 2>/dev/null; then
      PHARMCAT_JSON="$FILE"
      break 2
    fi
  done
done

if [ -z "$PHARMCAT_JSON" ]; then
  echo "ERROR: PharmCAT output not found. Run step 7 first."
  echo "  Expected in: ${GENOME_DIR}/${SAMPLE}/pharmcat/ or ${GENOME_DIR}/${SAMPLE}/vcf/"
  exit 1
fi

OUTDIR="${GENOME_DIR}/${SAMPLE}/cpic"
mkdir -p "$OUTDIR"
OUTPUT="${OUTDIR}/${SAMPLE}_cpic_recommendations.txt"

echo "============================================"
echo "  Step 27: CPIC Drug-Gene Recommendations"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${PHARMCAT_JSON}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""

# pypgx (step 32) calls CYP2D6 from read depth, which PharmCAT cannot. Its
# summary is compared with PharmCAT here rather than in step 32: run-all.sh
# starts steps 07 and 32 side by side, so step 32 could read a missing or
# previous-run report, while this step runs after both.
PYPGX_SUMMARY="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv"
COMPARISON="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pharmcat_comparison.tsv"
PYPGX_ARGS=()
rm -f "$COMPARISON"
if [ -f "$PYPGX_SUMMARY" ]; then
  echo "pypgx summary: ${PYPGX_SUMMARY} (comparison: ${COMPARISON})"
  PYPGX_ARGS=(--pypgx "$(cpath "$PYPGX_SUMMARY")" --comparison "$(cpath "$COMPARISON")")
else
  echo "pypgx summary not found (step 32 not run): no PharmCAT/pypgx comparison."
fi

# The parser, the drug lookup and the comparison live in bin/pgx_parse.py, the
# same code the CPIC_LOOKUP Nextflow module and the unit tests run. It exits 1
# when the report yields no gene, after writing a report that says so.
echo "Parsing PharmCAT results..."
rm -f "$OUTPUT" "${OUTDIR}/${SAMPLE}_phenotypes.tsv"
RC=0
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/pgx_parse.py cpic-report \
    --sample "$SAMPLE" \
    --report "$(cpath "$PHARMCAT_JSON")" \
    --outdir "$(cpath "$OUTDIR")" \
    ${PYPGX_ARGS[@]+"${PYPGX_ARGS[@]}"} || RC=$?

echo ""
echo "============================================"
echo "  CPIC recommendations: ${SAMPLE}"
echo "  Output: ${OUTPUT}"
echo "============================================"
echo ""
if [ -f "$OUTPUT" ]; then
  cat "$OUTPUT"
fi
if [ "$RC" -ne 0 ]; then
  echo "ERROR: the PharmCAT report ${PHARMCAT_JSON} could not be parsed into gene results (exit ${RC})." >&2
  exit "$RC"
fi
