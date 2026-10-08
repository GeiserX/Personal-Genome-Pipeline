#!/usr/bin/env bash
# run-all.sh passes the ancestry panel and the pgsc_calc checkout setup.sh
# installs: --ancestry_ref (step 26 inside the PRS run) and --pgsc_calc.
# ANCESTRY_PANEL=none and a missing panel skip step 26 with the reason;
# PGSC_CALC_DIR points at another checkout. ANCESTRY=true is no longer read
# and runs no script. Step 37 (y_haplogroup) is opt-in: listed as skipped by
# default, skipped without Yleaf's marker tables (setup.sh --yleaf-data), and
# passed in --tools with --yleaf_data when TOOLS names it and they are there.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
use_output_hook
seed_sample "$G" sample1
. "${REPO_ROOT}/versions.env"
PANEL="${G}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst"
PC="${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}"
mkdir -p "${G}/prs_scores" "$(dirname "$PANEL")" "$PC" "${CASE_WORK}/other-pgsc"
printf 'placeholder\n' | gzip -c > "${G}/prs_scores/PGS000018.txt.gz"
printf 'placeholder\n' > "$PANEL"
printf 'chr1\t1\tA\tG\n' > "${PANEL%.tar.zst}_GRCh38_sites.tsv"
: > "${PC}/main.nf"
: > "${CASE_WORK}/other-pgsc/main.nf"
last() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1; }
has_args() { grep -qF -- "$(printf '%q ' "$@")" <<<"$(last)" || fail "nextflow lacks '$*': $(last)"; }
lacks_arg() { if grep -qF -- " $1 " <<<"$(last)"; then fail "$1 passed $2: $(last)"; fi; }

run_expect 0 panel env ANCESTRY=true "${SCRIPTS}/run-all.sh" sample1 male
has_args --ancestry_ref "$PANEL"
has_args --pgsc_calc "$PC"
output_has panel '^  26 Ancestry \(pgsc_calc\) +runs$'
output_has panel 'NOTE: ANCESTRY is no longer read'
[ ! -e "${G}/sample1/logs/26-ancestry.log" ] || fail "ANCESTRY=true still ran scripts/26-ancestry.sh after the pipeline"
output_has panel '^  37 Y haplogroup \(Yleaf\) +skipped +\(opt-in: add y_haplogroup to TOOLS\)$'
grep -qE -- '--tools [^ ]*ancestry' <<<"$(last)" || fail "ancestry is not in --tools: $(last)"
if grep -qE -- '--tools [^ ]*y_haplogroup' <<<"$(last)"; then fail "y_haplogroup ran without being named"; fi

run_expect 0 panel-none env ANCESTRY_PANEL=none PGSC_CALC_DIR="${CASE_WORK}/other-pgsc" "${SCRIPTS}/run-all.sh" sample1 male
lacks_arg --ancestry_ref "with ANCESTRY_PANEL=none"
has_args --pgsc_calc "${CASE_WORK}/other-pgsc"
output_has panel-none '^  26 Ancestry \(pgsc_calc\) +skipped +\(ANCESTRY_PANEL=none\)$'

rm -f "${PANEL%.tar.zst}_GRCh38_sites.tsv"
run_expect 0 no-sites "${SCRIPTS}/run-all.sh" sample1 male
lacks_arg --ancestry_ref "without the panel's site list"
output_has no-sites '^  26 Ancestry \(pgsc_calc\) +skipped +\(data not installed: reference/pgsc_calc/'

run_expect 0 y-no-data env TOOLS=pharmcat,y_haplogroup "${SCRIPTS}/run-all.sh" sample1 male
output_has y-no-data '^  37 Y haplogroup \(Yleaf\) +skipped +\(data not installed: reference/yleaf-'
lacks_arg --yleaf_data "without Yleaf's marker tables"
YD="${G}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
mkdir -p "${YD}/hg38"
printf 'chry\tM1\tA-M1\t2781480\tC->T\tC\tT\n' > "${YD}/hg38/new_positions.txt"
run_expect 0 y env TOOLS=pharmcat,y_haplogroup "${SCRIPTS}/run-all.sh" sample1 male
has_args --tools pharmcat,y_haplogroup
has_args --yleaf_data "$YD"
output_has y '^  37 Y haplogroup \(Yleaf\) +runs$'
echo "The panel, pgsc_calc and step 37 reach nextflow as asked."
