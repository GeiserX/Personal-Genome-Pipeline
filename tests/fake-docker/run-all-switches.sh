#!/usr/bin/env bash
# The environment switches of run-all.sh map onto nextflow run flags, options
# after the sex reach nextflow unchanged, TOOLS narrows the step list, and the
# script-only steps (GRIDSS=true, EXTRA_CALLERS=...) run after the pipeline,
# their failure making run-all.sh exit 1. Unknown names in TOOLS or
# EXTRA_CALLERS stop it with exit 2 before nextflow.
#   plain run      the host's CPU count and RAM as --max_cpus and --max_memory;
#                  BENCHMARK=true with one caller VCF is skipped, not failed
#   ALIGN_DIR run  its BAM replaces the plain run's samplesheet row
#   plain again    back to aligned/: a kept row names this call's BAM only
#   publish mode   --publish_dir_mode link when a hard link from the work
#                  directory to GENOME_DIR works; copy (no flag) when ln
#                  fails, when -w names another work directory, or when the
#                  user gave --publish_dir_mode (passed once, as given)
#   -bg            refused with exit 2 before nextflow: run-all.sh would go on
#                  at once and write the step results before the run ended
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
lacks_arg() { if grep -qF -- " $1 " <<<"$(last)"; then fail "$1 passed $2: $(last)"; fi; }
row() { sed -n 2p "${G}/sample1/nextflow/samplesheet.csv"; }
A="${G}/sample1/aligned/sample1_sorted.bam" B="${G}/sample1/aligned_bwamem2/sample1_sorted.bam"

# A plain run: the host's caps, none of the switch flags
run_expect 0 plain env BENCHMARK=true "${SCRIPTS}/run-all.sh" sample1 male
has_args --max_cpus "$(host_cpus)" --max_memory "$(host_mem_gb).GB" --publish_dir_mode link
for f in --skip_trim --intervals; do lacks_arg "$f" "without its switch"; done
[ "$(row)" = "sample1,,,${A},${A}.bai,,,male" ] || fail "plain run row: $(row)"
output_has plain '^  benchmark-variants skipped \(only one caller VCF'
[ ! -e "${G}/sample1/logs/benchmark-variants.log" ] || fail "benchmark-variants.sh ran with one caller VCF"

# Every switch at once, after the plain run
run_rc switches env THREADS=4 SKIP_TRIM=true INTERVALS="chr20:1-100 chr22" ALIGN_DIR=aligned_bwamem2 TOOLS=pharmcat,cpic \
  GRIDSS=true EXTRA_CALLERS=gatk MAX_JOBS=3 "${SCRIPTS}/run-all.sh" sample1 male --max_memory 8.GB --sex_check warn
has_args --tools pharmcat,cpic
has_args --intervals "chr20:1-100 chr22"
has_args --max_cpus 4
has_args --skip_trim true
has_args --max_memory 8.GB --sex_check warn
[ "$(grep -o ' --max_memory ' <<<"$(last)" | wc -l)" -eq 1 ] || fail "--max_memory passed twice: $(last)"
[ "$(row)" = "sample1,,,${B},${B}.bai,,,male" ] || fail "ALIGN_DIR not used, the plain run's row was kept: $(row)"
output_has switches '^  07 PharmCAT +runs$'
output_has switches '^  11 ROH +skipped +\(not in TOOLS\)$'
output_has switches 'MAX_JOBS is no longer read'
# Steps 04b and 03a ran after the pipeline; they fail on the placeholder data,
# and run-all.sh reports that with exit 1
expect_rc switches 1
for s in 04b-gridss 03a-gatk-haplotypecaller; do [ -f "${G}/sample1/logs/${s}.log" ] || fail "scripts/${s}.sh did not run"; done
grep -qE $'^step\t04b\t(ok|failed)$' "${G}/sample1/logs/run_status.tsv" || fail "step 04b has no status line"
[ -s "${G}/sample1/logs/24-html-report.log" ] || fail "the reports did not run after the script-only steps"

# Without ALIGN_DIR again: the row from the ALIGN_DIR run is not reused
run_expect 0 plain2 env TOOLS=pharmcat "${SCRIPTS}/run-all.sh" sample1 male
[ "$(row)" = "sample1,,,${A},${A}.bai,,,male" ] || fail "the ALIGN_DIR run's row was kept without ALIGN_DIR: $(row)"

# The publish mode: copy when no hard link can be made, or when the user chose
mkdir -p "${CASE_WORK}/noln"
printf '#!/bin/sh\necho "ln: failed to create hard link: Invalid cross-device link" >&2\nexit 1\n' > "${CASE_WORK}/noln/ln"
chmod +x "${CASE_WORK}/noln/ln"
run_expect 0 crossdev env TOOLS=pharmcat PATH="${CASE_WORK}/noln:${PATH}" "${SCRIPTS}/run-all.sh" sample1 male
lacks_arg --publish_dir_mode "when ln fails"
run_expect 0 otherwork env TOOLS=pharmcat "${SCRIPTS}/run-all.sh" sample1 male -w "${CASE_WORK}/elsewhere"
lacks_arg --publish_dir_mode "with -w"
run_expect 0 usermode env TOOLS=pharmcat "${SCRIPTS}/run-all.sh" sample1 male --publish_dir_mode copy
has_args --publish_dir_mode copy
[ "$(grep -o ' --publish_dir_mode ' <<<"$(last)" | wc -l)" -eq 1 ] || fail "--publish_dir_mode passed twice: $(last)"
if compgen -G "${G}/.pgp-link-probe*" > /dev/null || compgen -G "${G}/sample1/nextflow/work/.pgp-link-probe*" > /dev/null; then
  fail "the link probe left files behind"
fi

n=$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")
run_expect 2 bg env TOOLS=pharmcat "${SCRIPTS}/run-all.sh" sample1 male -bg
output_has bg "ERROR: run-all.sh does not take -bg"
run_expect 2 badtool env TOOLS=pharmcat,clinvar_screen "${SCRIPTS}/run-all.sh" sample1 male
output_has badtool "unknown step 'clinvar_screen' in TOOLS"
run_expect 2 badcaller env EXTRA_CALLERS=gatk,bogus "${SCRIPTS}/run-all.sh" sample1 male
output_has badcaller "unknown caller 'bogus'"
[ "$(grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG")" -eq "$n" ] || fail "nextflow started despite an unknown name"
echo "Switches map onto nextflow flags; script-only steps run after the pipeline."
