#!/usr/bin/env bash
# Step 16b (mosdepth) reports the depth the fixture was sampled to.
. "$(dirname "$0")/lib.sh"

run_step 16b-mosdepth.sh "$SAMPLE"
check_step_exit 16b-mosdepth.sh

D="${GENOME_DIR}/${SAMPLE}/mosdepth"
check "summary exists" test -s "${D}/${SAMPLE}.mosdepth.summary.txt"
BINS=$(gzip -dc "${D}/${SAMPLE}.regions.bed.gz" 2>/dev/null \
  | awk '$1 == "chr20" && $2 >= 10000000 && $3 <= 10500000 && $4 >= 10' | wc -l | tr -d ' ')
check_ge "500 bp bins on the chr20 slice at 10x or more" "$BINS" 500

finish
