#!/usr/bin/env bash
# validate-setup.sh with Docker up, every required reference file present and
# a sample with a BAM and a VCF: all critical checks pass, exit 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

run_expect 0 validate "${SCRIPTS}/validate-setup.sh" sample1
output_lacks validate 'unbound variable'
output_has validate 'Docker images are pulled'
