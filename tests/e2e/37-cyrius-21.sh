#!/usr/bin/env bash
# Step 21 (Cyrius) runs and writes its TSV, whatever the CYP2D6 call is.
. "$(dirname "$0")/lib.sh"

run_step 21-cyrius.sh "$SAMPLE"
check_step_exit 21-cyrius.sh

TSV="${GENOME_DIR}/${SAMPLE}/cyrius/${SAMPLE}_cyp2d6.tsv"
check_ge "lines in the Cyrius TSV (header and sample)" "$(grep -c . "$TSV" 2>/dev/null || true)" 2
[ -s "$TSV" ] && cat "$TSV"

finish
