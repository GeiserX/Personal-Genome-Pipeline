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
echo "The samplesheet follows the inputs and stays the same across a rerun."
