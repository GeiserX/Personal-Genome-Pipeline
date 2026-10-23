#!/usr/bin/env bash
# validate-setup.sh requires ClinVar only for a run that screens it: with
# every other file present and no clinvar/ folder,
#   TOOLS unset (a default run, step 06 included)   [FAIL], exit 1
#   TOOLS=pharmcat,cpic (run-all.sh's step list)     [WARN], exit 0
#   TOOLS=pharmcat,clinvar                           [FAIL], exit 1
# and with ClinVar installed, TOOLS without clinvar still reports it present.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

run_expect 1 default "${SCRIPTS}/validate-setup.sh" sample1
output_has default '\[FAIL\].*ClinVar VCF not found'
output_has default '\[FAIL\].*ClinVar pathogenic subset not found'

TOOLS=pharmcat,cpic run_expect 0 no-clinvar "${SCRIPTS}/validate-setup.sh" sample1
output_has no-clinvar '\[WARN\].*ClinVar not installed .*TOOLS \(pharmcat,cpic\) has no clinvar'
output_lacks no-clinvar '\[FAIL\].*ClinVar'

TOOLS='pharmcat, clinvar' run_expect 1 with-clinvar "${SCRIPTS}/validate-setup.sh" sample1
output_has with-clinvar '\[FAIL\].*ClinVar VCF not found'

seed_clinvar "$GENOME_DIR"
TOOLS=pharmcat run_expect 0 installed "${SCRIPTS}/validate-setup.sh" sample1
output_has installed 'ClinVar VCF \(chr-prefixed\): present'
output_has installed 'ClinVar pathogenic subset: present'
