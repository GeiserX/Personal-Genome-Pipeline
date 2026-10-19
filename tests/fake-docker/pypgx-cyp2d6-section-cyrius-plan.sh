#!/usr/bin/env bash
# run-all.sh plans step 21 (Cyrius) in a default run when setup.sh --cyrius
# installed this version, and says why; a stamp of another image or lock file
# leaves it out with the reason, as does no install at all (that one is in
# base-run-all-default); TOOLS without cyrius leaves it out even when it is
# installed.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
case "${*:2}" in
  "cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json") echo '[{"LocusId": "HTT"}]' ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"
CY="${G}/tools/cyrius-${CYRIUS_VERSION}"
nf_args() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -n 1; }

run_expect 0 install "${SCRIPTS}/setup.sh" --cyrius "$G"
[ -s "${CY}/INSTALLED" ] || fail "setup.sh --cyrius wrote no INSTALLED stamp"

# --- installed: a default run plans it -----------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 planned "${SCRIPTS}/run-all.sh" sample1 male
output_has planned '^  21 Cyrius CYP2D6 +runs$'
output_has planned 'planned because setup\.sh --cyrius installed it'
nf_args | grep -qE -- "--tools [^ ]*,cyrius( |,)" || fail "--tools lacks cyrius: $(nf_args)"
nf_args | grep -qF -- "--cyrius_install ${CY} " || fail "nextflow did not get --cyrius_install ${CY}: $(nf_args)"
grep -q $'^step\t21\t' "${G}/sample1/logs/run_status.tsv" || fail "run_status.tsv has no step 21 row"
if grep -q $'^step\t21\tskipped' "${G}/sample1/logs/run_status.tsv"; then fail "step 21 recorded as skipped while installed"; fi

# --- TOOLS without cyrius: left out, even installed -----------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 tools env TOOLS=pharmcat,cpic "${SCRIPTS}/run-all.sh" sample1 male
output_has tools '^  21 Cyrius CYP2D6 +skipped +\(not in TOOLS\)$'
output_lacks tools 'planned because'
if nf_args | grep -qE -- '--tools [^ ]*cyrius|--cyrius_install'; then fail "TOOLS without cyrius still passed it: $(nf_args)"; fi

# --- a stamp of another lock file: left out, with the reason --------------------------
echo "python=${PYTHON_IMAGE} lock=0000" > "${CY}/INSTALLED"
: > "$FAKE_DOCKER_LOG"
run_expect 0 stale "${SCRIPTS}/run-all.sh" sample1 male
output_has stale "^  21 Cyrius CYP2D6 +skipped +\\(data not installed: tools/cyrius-${CYRIUS_VERSION//./\\.} for this version \\(setup\\.sh --cyrius\\)\\)$"
if nf_args | grep -qE -- '--tools [^ ]*cyrius|--cyrius_install'; then fail "a stale install still passed cyrius: $(nf_args)"; fi
echo "run-all.sh plans Cyrius when setup.sh --cyrius installed this version, and leaves it out otherwise."
