#!/usr/bin/env bash
# Step 11 (bcftools roh) writes per-site states for the sample.
. "$(dirname "$0")/lib.sh"

run_step 11-roh-analysis.sh "$SAMPLE"
check_step_exit 11-roh-analysis.sh

ROH="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt"
check_ge "per-site (ST) lines in the ROH output" "$(grep -c '^ST' "$ROH" 2>/dev/null || true)" 100

finish
