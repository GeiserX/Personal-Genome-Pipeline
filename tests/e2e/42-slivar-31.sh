#!/usr/bin/env bash
# Step 31 prioritises the VEP-annotated variants and runs the compound-het search.
. "$(dirname "$0")/lib.sh"

run_step 31-slivar.sh "$SAMPLE"
check_step_exit 31-slivar.sh

OUT=$(cat "$STEP_LOG")
check "no 'unbound variable' error" lacks 'unbound variable' "$OUT"
check "slivar compound-hets did not fail" lacks 'compound-hets failed' "$OUT"
check_ge "records in the prioritised VCF" "$(vcf_count "${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz")" 1
check "compound-het output is a readable VCF" vcf_ok "${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz"

finish
