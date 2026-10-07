#!/usr/bin/env bash
# Step 03 with `male`: no heterozygous call on the non-PAR chrX and chrY
# slices (HG002 is male), while the PAR1 slice stays diploid; a gVCF is
# written next to the VCF and the intermediate files are gone. Calls the BAM
# case 20 aligned, under another sample name so case 21's VCF stays as it is.
. "$(dirname "$0")/lib.sh"

X=HG002X
NONPAR="chrX:73700001-74000000,chrY:2700001-3000000"
PAR1="chrX:1000001-1200000"
mkdir -p "${GENOME_DIR}/${X}/aligned"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${X}/aligned/${X}_sorted.bam"
ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${X}/aligned/${X}_sorted.bam.bai"

INTERVALS="${NONPAR//,/ } ${PAR1}" run_step 03-deepvariant.sh "$X" male
check_step_exit 03-deepvariant.sh

VCF="${X}/vcf/${X}.vcf.gz"
GVCF="${X}/vcf/${X}.g.vcf.gz"
check "VCF is readable" vcf_ok "$VCF"
HETS=$(vcf_count -g het -r "$NONPAR" "$VCF")
ALTS=$(vcf_count -i 'GT="alt"' -r "$NONPAR" "$VCF")
PAR_HETS=$(vcf_count -g het -r "$PAR1" "$VCF")
echo "non-PAR chrX/chrY: ${ALTS} records with an ALT allele, ${HETS} heterozygous; PAR1: ${PAR_HETS} heterozygous"
bcf view -H -g het -r "$NONPAR" "$VCF" 2>/dev/null | head -n 5
check_ge "calls with an ALT allele on the non-PAR slices (the check below has calls to judge)" "$ALTS" 5
check_eq "heterozygous calls on non-PAR chrX and chrY of a male sample" "$HETS" 0
check_ge "heterozygous calls on the PAR1 slice (diploid there)" "$PAR_HETS" 1
check "gVCF is readable" vcf_ok "$GVCF"
check "gVCF index exists" nonempty "${GVCF}.tbi"
check_ge "gVCF reference blocks (records with END)" "$(vcf_count -i 'INFO/END>0' "$GVCF")" 10
check "the step log names the PAR BED" has 'par_grch38\.bed' "$(cat "$STEP_LOG")"
check_eq "intermediate and .part files left in ${X}/vcf" \
  "$(find "${GENOME_DIR}/${X}/vcf" -mindepth 1 \( -name 'deepvariant_tmp' -o -name '*.part.*' \) | wc -l | tr -d ' ')" 0
rm -rf "${GENOME_DIR:?}/${X}"

finish
