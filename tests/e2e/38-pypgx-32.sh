#!/usr/bin/env bash
# Step 32 (pypgx) returns a CYP2D6 result; the fixture has CYP2D6 and VDR (its
# control gene). The bundle branch follows the pinned pypgx image.
. "$(dirname "$0")/lib.sh"

BUNDLE="${GENOME_DIR}/reference/pypgx-bundle"
PYPGX_VERSION="${PYPGX_IMAGE##*:}"
PYPGX_VERSION="${PYPGX_VERSION%%--*}"
if [ ! -d "$BUNDLE" ]; then
  git clone -q --branch "$PYPGX_VERSION" --depth 1 https://github.com/sbslee/pypgx-bundle.git "$BUNDLE"
fi

run_step 32-pypgx.sh "$SAMPLE"
check_step_exit 32-pypgx.sh

SUMMARY_TSV="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv"
CYP2D6=$(awk -F'\t' '$1 == "CYP2D6" {print $2; exit}' "$SUMMARY_TSV" 2>/dev/null)
echo "CYP2D6: ${CYP2D6:-none}"
check "pypgx returns a CYP2D6 row" test -n "${CYP2D6:-}"
check "the CYP2D6 result is not FAILED or N/A" lacks '^(FAILED|N/A)$' "${CYP2D6:-FAILED}"

finish
