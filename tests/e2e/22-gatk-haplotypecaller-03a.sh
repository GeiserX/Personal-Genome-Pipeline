#!/usr/bin/env bash
# Step 03a (GATK HaplotypeCaller) on the chr20 slice: GATK takes the sample from @RG.
. "$(dirname "$0")/lib.sh"

INTERVALS=chr20:10000000-10500000 run_step 03a-gatk-haplotypecaller.sh "$SAMPLE"
check_step_exit 03a-gatk-haplotypecaller.sh

VCF="${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
check "VCF is readable" vcf_ok "$VCF"
check_ge "records on the chr20 slice" "$(vcf_count "$VCF")" 200
check_eq "VCF sample column" "$(bcf query -l "$VCF" 2>/dev/null)" "$SAMPLE"

finish
