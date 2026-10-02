#!/usr/bin/env bash
# Step 06 finds the planted pathogenic record and reports it with its gene.
. "$(dirname "$0")/lib.sh"

run_step 06-clinvar-screen.sh "$SAMPLE"
check_step_exit 06-clinvar-screen.sh

GENE=$(planted gene)
DIR="${GENOME_DIR}/${SAMPLE}/clinvar"
check_ge "planted $(planted chrom):$(planted pos) is among the sample's ClinVar hits" \
  "$(sample_side_hits "$DIR" "")" 1
check_ge "the hit names its gene (${GENE}) in a file that holds the sample's call" \
  "$(sample_side_hits "$DIR" "$GENE")" 1

finish
