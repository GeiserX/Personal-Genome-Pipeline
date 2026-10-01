#!/usr/bin/env bash
# run-all.sh with the default settings (SKIP_VALIDATION unset, no optional
# data) on a sample that already has a BAM and a VCF: the pre-flight
# validation passes and the run exits 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

run_expect 0 run-all "${SCRIPTS}/run-all.sh" sample1 male
output_lacks run-all 'Setup validation failed'
output_lacks run-all 'unbound variable'
