#!/usr/bin/env bash
# setup.sh when only the raw ClinVar file and its index exist (an earlier
# run stopped before the derived files): it must build clinvar_chr.vcf.gz
# and clinvar_pathogenic_chr.vcf.gz instead of reporting ClinVar as done.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/clinvar"
printf 'placeholder\n' > "${GENOME_DIR}/clinvar/clinvar.vcf.gz"
printf 'placeholder\n' > "${GENOME_DIR}/clinvar/clinvar.vcf.gz.tbi"
use_output_hook

run_rc setup "${SCRIPTS}/setup.sh" "$GENOME_DIR"
awk '/^run / && /clinvar_chr\.vcf\.gz/ && !/clinvar_pathogenic/ { found = 1 } END { exit !found }' \
  "$FAKE_DOCKER_LOG" || fail "setup.sh did not build clinvar_chr.vcf.gz from the raw ClinVar file"
docker_log_has '^run .*clinvar_pathogenic_chr\.vcf\.gz' \
  "setup.sh did not build clinvar_pathogenic_chr.vcf.gz"
output_lacks setup 'unbound variable'
expect_rc setup 0
