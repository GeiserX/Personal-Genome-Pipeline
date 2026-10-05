#!/usr/bin/env bash
# 36-pgx-consensus.sh — What the BAM-based callers tell PharmCAT
# Usage: ./scripts/36-pgx-consensus.sh <sample_name>
#
# PharmCAT (step 07) reads a VCF, and from a VCF it types neither HLA nor
# CYP2D6. This step writes the outside-call file PharmCAT takes with -po:
#   HLA-A, HLA-B  T1K's types (step 08), cut to two fields (*57:01)
#   CYP2D6        only when pypgx (step 32) and Cyrius (step 21, opt-in) give
#                 the same diplotype and the CYP2D6 depth check passed; one
#                 caller alone, a disagreement or multi-mapped depth leaves it
#                 'indeterminate', and nothing reaches PharmCAT
# and a consensus table that says, for each gene, what every caller said and
# what was passed on (bin/pgx_outside_calls.py, the same code as the
# PGX_CONSENSUS Nextflow module).
#
# Run it after steps 08, 32 and (opt-in) 21, then steps 07 and 27: step 07
# gives PharmCAT the file when it exists, and step 27 lists the outside calls
# in the CPIC recommendations.
#
# Output: ${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_outside_calls.tsv
#         ${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

S="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${S}/pgx_consensus"
CALLS="${OUTDIR}/${SAMPLE}_outside_calls.tsv"
CONSENSUS="${OUTDIR}/${SAMPLE}_pgx_consensus.tsv"
mkdir -p "$OUTDIR"

echo "=== PGx consensus (outside calls for PharmCAT): ${SAMPLE} ==="

# input LABEL OPTION FILE: pass FILE to the script when it exists.
ARGS=()
input() {
  if [ -f "$3" ]; then
    echo "  ${1}: ${3}"
    ARGS+=("$2" "$(cpath "$3")")
  else
    echo "  ${1}: not found (${3})"
  fi
}
input "HLA types (step 08)" --hla "${S}/hla_t1k/${SAMPLE}_hla_genotype.tsv"
input "pypgx (step 32)" --pypgx "${S}/pypgx/${SAMPLE}_pypgx_summary.tsv"
input "Cyrius (step 21)" --cyrius "${S}/cyrius/${SAMPLE}_cyp2d6.tsv"
# Both CYP2D6 steps write the same depth check; step 32's, else step 21's.
DEPTH="${S}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv"
[ -f "$DEPTH" ] || DEPTH="${S}/cyrius/${SAMPLE}_cyp2d6_depth_check.tsv"
input "CYP2D6 depth check" --depth-check "$DEPTH"

rm -f "$CALLS" "$CONSENSUS"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/pgx_outside_calls.py \
    --calls "$(cpath "$CALLS")" \
    --consensus "$(cpath "$CONSENSUS")" \
    ${ARGS[@]+"${ARGS[@]}"}

echo ""
echo "Outside calls for PharmCAT (${CALLS}):"
if [ -s "$CALLS" ]; then
  sed 's/^/  /' "$CALLS"
else
  echo "  none"
fi
echo ""
echo "=== PGx consensus complete ==="
echo "Table: ${CONSENSUS}"
echo "Next: steps 07 (PharmCAT reads the outside calls) and 27 (CPIC lookup)."
