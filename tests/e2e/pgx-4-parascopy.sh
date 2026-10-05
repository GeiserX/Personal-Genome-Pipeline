#!/usr/bin/env bash
# Step 35 (Parascopy, opt-in) gives HG002's SMN1/SMN2 locus a copy number
# with its quality, on the fixture's chr5 slice. Its data comes from
# `setup.sh --parascopy-data`. The fixture BAM holds reads only in the
# slices, so the background depth comes from 100 bp windows over the chr20
# slice (PARASCOPY_DEPTH_BED) instead of Parascopy's genome-wide windows.
. "$(dirname "$0")/lib.sh"

echo "+ scripts/setup.sh --parascopy-data"
"${REPO}/scripts/setup.sh" --parascopy-data "$GENOME_DIR"
check "setup.sh --parascopy-data installed the homology table" \
  test -s "${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}/homology_table/GRCh38.bed.gz"

BED="${GENOME_DIR}/${SAMPLE}/parascopy_windows.bed"
awk 'BEGIN { for (s = 10050000; s < 10450000; s += 100) printf "chr20\t%d\t%d\n", s, s + 100 }' > "$BED"
PARASCOPY_DEPTH_BED="$BED" run_step 35-paralogs.sh "$SAMPLE"
check_step_exit 35-paralogs.sh

TSV="${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_smn_copy_number.tsv"
cat "$TSV" 2>/dev/null
check_ge "SMN1/SMN2 rows on chr5" "$(awk -F'\t' 'NR > 1 && $1 == "chr5"' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 1
# shellcheck disable=SC2016  # awk fields
check "a row with an aggregate copy number and its quality" \
  awk -F'\t' 'NR > 1 && $1 == "chr5" && $6 != "" && $6 != "*" && $7 ~ /^[0-9.]+$/ {found = 1} END {exit !found}' "$TSV"
echo "- Parascopy on the fixture: $(awk -F'\t' 'NR > 1 && $1 == "chr5" {printf "%s:%s-%s agCN %s (Q%s) psCN %s; ", $1, $2, $3, $6, $7, $9}' "$TSV" 2>/dev/null)" >> "$E2E_NOTES"

finish
