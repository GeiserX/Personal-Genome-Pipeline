#!/usr/bin/env bash
# Step 12 (haplogrep3) assigns a haplogroup from the chrM calls.
. "$(dirname "$0")/lib.sh"

run_step 12-mito-haplogroup.sh "$SAMPLE"
check_step_exit 12-mito-haplogroup.sh

OUTF="${GENOME_DIR}/${SAMPLE}/mito/${SAMPLE}_haplogroup.txt"
check "haplogroup file exists" test -s "$OUTF"
HG=$(awk -F'\t' 'NR == 2 {gsub(/"/, "", $2); print $2}' "$OUTF" 2>/dev/null)
check "haplogroup column is filled (${HG:-empty})" test -n "${HG:-}"

finish
