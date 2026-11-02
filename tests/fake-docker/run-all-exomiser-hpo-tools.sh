#!/usr/bin/env bash
# run-all.sh with HPO set and a TOOLS list that leaves exomiser out: a set
# HPO turns step 48 (Exomiser) on, as --hpo does in main.nf, so the plan
# lists it as running and nextflow gets --hpo, the Exomiser data and
# exomiser in --tools. Without the Exomiser data it is skipped as not
# installed, never as "not in TOOLS".
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
ED="${G}/reference/exomiser-${EXOMISER_DATA_VERSION}"
NFL() { awk '/^nextflow :: / {l = $0} END {print l}' "$FAKE_DOCKER_LOG"; }

# 1. no Exomiser data: skipped for the data, not for TOOLS
run_expect 0 hpo-no-data env TOOLS=hla_typing HPO=HP:0001250 "${SCRIPTS}/run-all.sh" sample1 male
output_lacks hpo-no-data 'unbound variable'
output_has hpo-no-data '^  48 Exomiser \(HPO terms\) +skipped +\(data not installed: '
output_lacks hpo-no-data '48 Exomiser.*not in TOOLS'

# 2. with the data: step 48 runs
mkdir -p "${ED}/${EXOMISER_DATA_VERSION}_hg38" "${ED}/${EXOMISER_DATA_VERSION}_phenotype"
: > "$FAKE_DOCKER_LOG"
run_expect 0 hpo-tools env TOOLS=hla_typing HPO=HP:0001250 "${SCRIPTS}/run-all.sh" sample1 male
output_has hpo-tools '^  48 Exomiser \(HPO terms\) +runs$'
L=$(NFL)
[[ "$L" == *" --hpo HP:0001250 --exomiser_data ${ED} "* ]] || fail "nextflow was not given --hpo and the Exomiser data: ${L}"
[[ "$L" =~ --tools\ ([^ ]*) ]] && [[ ",${BASH_REMATCH[1]}," == *",exomiser,"* ]] \
  || fail "exomiser is not in the --tools nextflow got: ${L}"
grep -q $'^step\t48\tskipped' "${G}/sample1/logs/run_status.tsv" && fail "run_status.tsv records step 48 as skipped"
echo "HPO with a TOOLS list that leaves exomiser out still plans and passes step 48."
