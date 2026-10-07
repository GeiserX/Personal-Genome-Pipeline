#!/usr/bin/env bash
# run-all.sh without Nextflow, without Java or with a Java older than 17: it
# stops with exit 2 and the pinned install line, before validation, any
# container or nextflow.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
use_output_hook   # the reports at the end read the files their containers write
seed_sample "$GENOME_DIR" sample1
PIN=$(sed -n 's/^NEXTFLOW_VERSION="\(.*\)"/\1/p' "${REPO_ROOT}/versions.env")
[ -n "$PIN" ] || fail "no NEXTFLOW_VERSION in versions.env"
started() { awk '/^(run|nextflow|info|version) / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; }

FAKE_JAVA_VERSION=11.0.22 run_expect 2 java11 "${SCRIPTS}/run-all.sh" sample1 male
output_has java11 'needs Java 17 or later \(found: 11\)'
output_has java11 "NXF_VER=${PIN} bash"
FAKE_JAVA_VERSION=1.8.0_402 run_expect 2 java8 "${SCRIPTS}/run-all.sh" sample1 male
output_has java8 '\(found: 8\)'
FAKE_JAVA_VERSION=21.0.2 run_expect 0 java21 env SKIP_VALIDATION=true "${SCRIPTS}/run-all.sh" sample1 male
: > "$FAKE_DOCKER_LOG"

( hide_commands nextflow
  run_expect 2 nonextflow "${SCRIPTS}/run-all.sh" sample1 male )
output_has nonextflow 'Nextflow \(not on PATH\)'
output_has nonextflow "curl -s https://get.nextflow.io \| NXF_VER=${PIN} bash"
( hide_commands java
  run_expect 2 nojava "${SCRIPTS}/run-all.sh" sample1 male )
output_has nojava '\(found: none\)'
if started; then fail "run-all.sh validated or started something without its prerequisites: $(cat "$FAKE_DOCKER_LOG")"; fi
echo "Without Java 17 or Nextflow run-all.sh stops with exit 2 and the pinned install line."
