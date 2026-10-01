#!/usr/bin/env bash
# setup.sh on a machine without wget (stock macOS) but with curl: the
# downloads must go through curl and setup must finish.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

G="${CASE_WORK}/genome"
mkdir -p "$G"
hide_commands wget

run_rc setup "${SCRIPTS}/setup.sh" "$G"
docker_log_has '^curl [^ ]*\.fa' "setup.sh did not download the reference with curl"
docker_log_has '^curl [^ ]*clinvar\.vcf\.gz ' "setup.sh did not download ClinVar with curl"
output_lacks setup 'wget: command not found'
output_lacks setup 'unbound variable'
expect_rc setup 0
