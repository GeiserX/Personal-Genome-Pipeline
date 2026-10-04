#!/usr/bin/env bash
# The environment switches of run-all.sh map onto nextflow run flags, options
# after the sex reach nextflow unchanged, TOOLS narrows the step list, and the
# script-only steps (GRIDSS=true, EXTRA_CALLERS=...) run after the pipeline,
# their failure making run-all.sh exit 1. Unknown names in TOOLS or
# EXTRA_CALLERS stop it with exit 2 before nextflow.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook   # the reports at the end read the files their containers write
seed_sample "$G" sample1
mkdir -p "${G}/sample1/aligned_bwamem2"
cp "${G}/sample1/aligned/sample1_sorted.bam" "${G}/sample1/aligned/sample1_sorted.bam.bai" "${G}/sample1/aligned_bwamem2/"
rm -f "${G}/sample1/vcf/sample1.vcf.gz" "${G}/sample1/vcf/sample1.vcf.gz.tbi"
last() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1; }
has_args() { grep -qF -- "$(printf '%q ' "$@")" <<<"$(last)" || fail "nextflow lacks '$*': $(last)"; }

run_rc switches env THREADS=4 SKIP_TRIM=true INTERVALS="chr20:1-100 chr22" ALIGN_DIR=aligned_bwamem2 TOOLS=pharmcat,cpic \
  GRIDSS=true EXTRA_CALLERS=gatk MAX_JOBS=3 "${SCRIPTS}/run-all.sh" sample1 male --max_memory 8.GB --sex_check warn
has_args --tools pharmcat,cpic
has_args --intervals "chr20:1-100 chr22"
has_args --max_cpus 4
has_args --skip_trim true
has_args --max_memory 8.GB --sex_check warn
B="${G}/sample1/aligned_bwamem2/sample1_sorted.bam"
[ "$(sed -n 2p "${G}/sample1/nextflow/samplesheet.csv")" = "sample1,,,${B},${B}.bai,,,male" ] || fail "ALIGN_DIR not used: $(cat "${G}/sample1/nextflow/samplesheet.csv")"
output_has switches '^  07 PharmCAT +runs$'
output_has switches '^  11 ROH +skipped +\(not in TOOLS\)$'
output_has switches 'MAX_JOBS is no longer read'
# Steps 04b and 03a ran after the pipeline; they fail on the placeholder data,
# and run-all.sh reports that with exit 1
expect_rc switches 1
for s in 04b-gridss 03a-gatk-haplotypecaller; do [ -f "${G}/sample1/logs/${s}.log" ] || fail "scripts/${s}.sh did not run"; done
grep -qE $'^step\t04b\t(ok|failed)$' "${G}/sample1/logs/run_status.tsv" || fail "step 04b has no status line"
[ -s "${G}/sample1/logs/24-html-report.log" ] || fail "the reports did not run after the script-only steps"

# Without the switches none of those flags is passed
run_expect 0 plain "${SCRIPTS}/run-all.sh" sample1 male
for f in --max_cpus --skip_trim --intervals --max_memory; do
  if grep -qF -- " ${f} " <<<"$(last)"; then fail "${f} passed without its switch: $(last)"; fi
done

n=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 2 badtool env TOOLS=pharmcat,clinvar_screen "${SCRIPTS}/run-all.sh" sample1 male
output_has badtool "unknown step 'clinvar_screen' in TOOLS"
run_expect 2 badcaller env EXTRA_CALLERS=gatk,bogus "${SCRIPTS}/run-all.sh" sample1 male
output_has badcaller "unknown caller 'bogus'"
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$n" ] || fail "nextflow started despite an unknown name"
echo "Switches map onto nextflow flags; script-only steps run after the pipeline."
