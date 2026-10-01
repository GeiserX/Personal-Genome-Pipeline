#!/usr/bin/env bash
# setup.sh on an empty GENOME_DIR: the reference and ClinVar downloads go
# through the fake wget, the image checks through the fake docker. It must
# finish with exit 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

G="${CASE_WORK}/genome"
mkdir -p "$G"

run_expect 0 setup "${SCRIPTS}/setup.sh" "$G"
output_lacks setup 'unbound variable'
output_has setup 'Setup complete!'

for f in reference/Homo_sapiens_assembly38.fasta reference/Homo_sapiens_assembly38.fasta.fai \
         clinvar/clinvar.vcf.gz clinvar/clinvar.vcf.gz.tbi; do
  [ -s "${G}/${f}" ] || fail "setup.sh did not download ${f}"
done
grep -q '^image inspect' "$FAKE_DOCKER_LOG" || fail "setup.sh never checked a docker image"
