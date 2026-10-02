#!/usr/bin/env bash
# Step 23 reads CSQ from step 30's output and keeps the impactful variants.
. "$(dirname "$0")/lib.sh"

run_step 23-clinical-filter.sh "$SAMPLE"
check_step_exit 23-clinical-filter.sh

check "the step found CSQ annotations" lacks 'No CSQ' "$(cat "$STEP_LOG")"
VCF="${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz"
check "clinical VCF is readable" vcf_ok "$VCF"
check_ge "records in the clinical VCF" "$(vcf_count "$VCF")" 1

finish
