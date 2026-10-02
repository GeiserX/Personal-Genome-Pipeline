#!/usr/bin/env bash
# Step 20 (GATK Mutect2, mitochondrial mode) calls chrM variants from the step-02 BAM.
. "$(dirname "$0")/lib.sh"

run_step 20-mtoolbox.sh "$SAMPLE"
check_step_exit 20-mtoolbox.sh

VCF="${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
check "filtered chrM VCF is readable" vcf_ok "$VCF"
check_ge "chrM records" "$(vcf_count "$VCF")" 5

finish
