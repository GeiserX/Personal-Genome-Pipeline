#!/usr/bin/env bash
# A sample with a BAM, a VCF and the gVCF step 03 writes beside it
# (vcf/<sample>.g.vcf.gz): run-all.sh starts from the BAM and the VCF, not
# calling again, and gives the pipeline the gVCF in the gvcf and gvcf_index
# columns, so PharmCAT and PRS read it. Without the gVCF (or without its
# index) the row has no gvcf columns and the note says what is lost. A rerun
# writes the same samplesheet, so -resume finds the tasks again.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
sheet() { cat "${G}/$1/nextflow/samplesheet.csv"; }

seed_sample "$G" s1
B="${G}/s1/aligned/s1_sorted.bam" V="${G}/s1/vcf/s1.vcf.gz" GV="${G}/s1/vcf/s1.g.vcf.gz"
printf 'placeholder\n' > "$GV"
printf 'placeholder\n' > "${GV}.tbi"
want=$(printf 'sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex,gvcf,gvcf_index\ns1,,,%s,%s,%s,%s,male,%s,%s' \
  "$B" "${B}.bai" "$V" "${V}.tbi" "$GV" "${GV}.tbi")
run_expect 0 gvcf "${SCRIPTS}/run-all.sh" s1 male
[ "$(sheet s1)" = "$want" ] || fail "BAM+VCF+gVCF samplesheet: $(sheet s1)"
output_has gvcf 'PharmCAT and PRS read the gVCF beside it'
output_lacks gvcf 'a gVCF beside it is not read'
run_expect 0 rerun "${SCRIPTS}/run-all.sh" s1 male
[ "$(sheet s1)" = "$want" ] || fail "the rerun changed the samplesheet: $(sheet s1)"

# The gVCF without its index is not an input; nor is a missing gVCF.
rm -f "${GV}.tbi"
run_expect 0 notbi "${SCRIPTS}/run-all.sh" s1 male
[ "$(sheet s1)" = "$(printf 'sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex\ns1,,,%s,%s,%s,%s,male' "$B" "${B}.bai" "$V" "${V}.tbi")" ] \
  || fail "a gVCF without its index was used: $(sheet s1)"
output_has notbi 'no gVCF with its index beside it'
echo "A BAM+VCF row carries the gVCF beside the VCF, and only with its index."
