#!/usr/bin/env bash
# Which input run-all.sh gives the pipeline, and that a rerun gives it the same
# one: Nextflow's -resume finds a task again only when its inputs are the same
# files, so the samplesheet of the last run is kept while the files it names
# exist, even after the pipeline published a BAM and a VCF next to the FASTQ.
#   FASTQ only              -> a FASTQ row
#   rerun, BAM and VCF now  -> the same FASTQ row
#   BAM and VCF             -> a BAM+VCF row (not called again)
#   VCF of that row removed -> a BAM row (called again, with a gVCF)
#   nothing                 -> exit 1, nextflow not started
#   VCF only                -> a VCF row; the steps that read a BAM are
#                              skipped ("no BAM"), not recorded ok
#   a file without its index is not an input:
#     VCF without .tbi      -> a BAM row (called again)
#     BAM without .bai      -> the FASTQ row (aligned again)
#     neither index, no FASTQ -> exit 1, nextflow not started
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook   # the reports at the end read the files their containers write
row() { sed -n 2p "${G}/$1/nextflow/samplesheet.csv"; }
calls() { grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG" || true; }

# FASTQ only, then a rerun after the pipeline wrote the BAM and the VCF
mkdir -p "${G}/s1/fastq"
for r in R1 R2; do printf 'placeholder\n' > "${G}/s1/fastq/s1_${r}.fastq.gz"; done
run_expect 0 fastq "${SCRIPTS}/run-all.sh" s1 female
want="s1,${G}/s1/fastq/s1_R1.fastq.gz,${G}/s1/fastq/s1_R2.fastq.gz,,,,,female"
[ "$(row s1)" = "$want" ] || fail "FASTQ row: $(row s1)"
seed_sample "$G" s1
run_expect 0 rerun "${SCRIPTS}/run-all.sh" s1 female
[ "$(row s1)" = "$want" ] || fail "the rerun changed the samplesheet: $(row s1)"
grep -q -- '-resume' <<<"$(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1)" || fail "the rerun has no -resume"
[ "$(calls)" -eq 2 ] || fail "nextflow calls: $(calls), expected 2"

# BAM and VCF; then the VCF removed
seed_sample "$G" s2
B="${G}/s2/aligned/s2_sorted.bam" V="${G}/s2/vcf/s2.vcf.gz"
run_expect 0 bamvcf "${SCRIPTS}/run-all.sh" s2 male
[ "$(row s2)" = "s2,,,${B},${B}.bai,${V},${V}.tbi,male" ] || fail "BAM+VCF row: $(row s2)"
rm -f "$V" "${V}.tbi"
run_expect 0 bamonly "${SCRIPTS}/run-all.sh" s2 male
[ "$(row s2)" = "s2,,,${B},${B}.bai,,,male" ] || fail "BAM row after the VCF was removed: $(row s2)"
output_lacks bamonly 'starting from the existing VCF'

# No input at all
mkdir -p "${G}/s3"
run_expect 1 none "${SCRIPTS}/run-all.sh" s3 male
output_has none 'ERROR: no input for s3'
[ "$(calls)" -eq 4 ] || fail "nextflow started without an input: $(calls) calls"

# VCF only
seed_sample "$G" s4
rm -rf "${G}/s4/aligned"
V="${G}/s4/vcf/s4.vcf.gz"
run_expect 0 vcfonly "${SCRIPTS}/run-all.sh" s4 female
[ "$(row s4)" = "s4,,,,,${V},${V}.tbi,female" ] || fail "VCF row: $(row s4)"
output_has vcfonly '^  07 PharmCAT +runs$'
for s in '04 Manta' '16b mosdepth' '28 MultiQC'; do output_has vcfonly "^  ${s} +skipped +\(no BAM\)$"; done
ST=$(cat "${G}/s4/logs/run_status.tsv")
grep -q $'^step\t07\tok$' <<<"$ST" || fail "step 07 not ok: ${ST}"
for s in 16 16b 04 19 28; do
  grep -q $'^step\t'"${s}"$'\tskipped (no BAM)$' <<<"$ST" || fail "step ${s} is not 'skipped (no BAM)': ${ST}"
done
grep -qF -- " $(printf '%q ' --tools pharmcat,cpic,roh,mito_haplogroup,clinvar)--" <<<"$(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1)" \
  || fail "a BAM step is in --tools: $(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1)"
# The index decides: a BAM or a VCF without its index is not an input
seed_sample "$G" s5
mkdir -p "${G}/s5/fastq"
for r in R1 R2; do printf 'placeholder\n' > "${G}/s5/fastq/s5_${r}.fastq.gz"; done
B="${G}/s5/aligned/s5_sorted.bam" V="${G}/s5/vcf/s5.vcf.gz"
rm -f "${V}.tbi"
run_expect 0 notbi "${SCRIPTS}/run-all.sh" s5 male
[ "$(row s5)" = "s5,,,${B},${B}.bai,,,male" ] || fail "a VCF without its index was used: $(row s5)"
rm -f "${B}.bai"
run_expect 0 nobai "${SCRIPTS}/run-all.sh" s5 male
[ "$(row s5)" = "s5,${G}/s5/fastq/s5_R1.fastq.gz,${G}/s5/fastq/s5_R2.fastq.gz,,,,,male" ] || fail "a BAM without its index was used: $(row s5)"
rm -rf "${G}/s5/fastq"
n=$(calls)
run_expect 1 noindex "${SCRIPTS}/run-all.sh" s5 male
output_has noindex 'ERROR: no input for s5'
[ "$(calls)" -eq "$n" ] || fail "nextflow started on a BAM and a VCF without their indexes"
echo "The samplesheet follows the inputs and stays the same across a rerun."
