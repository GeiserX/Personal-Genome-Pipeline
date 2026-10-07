#!/usr/bin/env bash
# run-all.sh stops on a failed validation (exit 1, nextflow not started) and
# passes on the pipeline's exit code (the reports do not run, and no step is
# recorded as ok). SKIP_VALIDATION=true skips the validation.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_clinvar "$G"
use_output_hook   # the reports at the end read the files their containers write
seed_sample "$G" sample1
calls() { grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG" || true; }

# No reference: validate-setup.sh fails
run_expect 1 novalidate "${SCRIPTS}/run-all.sh" sample1 male
output_has novalidate '=== Summary ==='
output_has novalidate 'ERROR: setup validation failed'
[ "$(calls)" -eq 0 ] || fail "nextflow started after a failed validation"

# The same, with the validation skipped: nextflow starts
run_expect 0 skipvalidate env SKIP_VALIDATION=true "${SCRIPTS}/run-all.sh" sample1 male
output_lacks skipvalidate '=== Summary ==='
[ "$(calls)" -eq 1 ] || fail "nextflow calls: $(calls), expected 1"

# The pipeline fails with exit 3
rm -rf "${G}/sample1/logs"
run_expect 3 nffail env SKIP_VALIDATION=true FAKE_NEXTFLOW_RC=3 "${SCRIPTS}/run-all.sh" sample1 male
output_has nffail 'ERROR: the pipeline failed \(exit 3\)'
output_has nffail 'nextflow/\.nextflow\.log'
[ ! -e "${G}/sample1/logs/24-html-report.log" ] || fail "the HTML report ran after the pipeline failed"
if grep -q $'\tok$' "${G}/sample1/logs/run_status.tsv"; then fail "a step is recorded ok after the pipeline failed"; fi
echo "A failed validation stops run-all.sh before nextflow; the pipeline's exit code is passed on."
