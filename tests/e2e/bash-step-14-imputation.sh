#!/usr/bin/env bash
# Step 14 writes one PASS-only, indexed VCF per chromosome, chr1-22 and chrX,
# from one container.
. "$(dirname "$0")/lib.sh"

run_step 14-imputation-prep.sh "$SAMPLE"
check_step_exit 14-imputation-prep.sh
D="${SAMPLE}/imputation/mis_ready"
check_eq "chromosome VCFs" "$(find "${GENOME_DIR}/${D}" -maxdepth 1 -name "${SAMPLE}_chr*.vcf.gz" ! -name '*.part.*' | wc -l | tr -d ' ')" 23
check_eq "indexes" "$(find "${GENOME_DIR}/${D}" -maxdepth 1 -name "${SAMPLE}_chr*.vcf.gz.tbi" | wc -l | tr -d ' ')" 23
check_eq "leftover .part files" "$(find "${GENOME_DIR}/${D}" -maxdepth 1 -name '*.part*' | wc -l | tr -d ' ')" 0
check_ge "chr20 records" "$(vcf_count "${D}/${SAMPLE}_chr20.vcf.gz")" 100
check_ge "chrX records" "$(vcf_count "${D}/${SAMPLE}_chrX.vcf.gz")" 1
check_eq "chr20 records that did not PASS" "$(vcf_count -i 'FILTER!="PASS" && FILTER!="."' "${D}/${SAMPLE}_chr20.vcf.gz")" 0
check_eq "chr20 file holds only chr20" "$(bcf query -f '%CHROM\n' "${D}/${SAMPLE}_chr20.vcf.gz" 2>/dev/null | sort -u | paste -sd, -)" chr20

finish
