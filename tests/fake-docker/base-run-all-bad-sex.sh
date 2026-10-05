#!/usr/bin/env bash
# run-all.sh with a sex that is neither male nor female: it must stop at once
# with exit 2 and a usage line, before validation, any container or nextflow.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

run_rc run-all "${SCRIPTS}/run-all.sh" sample1 x
if awk '/^run |^nextflow / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "run-all.sh started a container or nextflow although the sex 'x' is invalid"
fi
expect_rc run-all 2
output_has run-all '[Uu]sage'
output_has run-all "sex must be 'male' or 'female', got 'x'"
