#!/usr/bin/env bash
# The PGx steps of package pgx-outside-calls-and-paralogs, against the fake docker:
#   setup.sh --cyrius   installs Cyrius from scripts/cyrius-constraints.txt
#                       with hashes (--require-hashes --no-deps
#                       --only-binary), with network, into tools/cyrius-<v>
#   step 21             refuses to run before that install, and after it runs
#                       every container with --network none (Cyrius from the
#                       install, not from PyPI)
#   step 36             writes the outside calls with bin/pgx_outside_calls.py
#   step 07             gives PharmCAT -po with a non-empty outside-call file,
#                       and no -po with an empty one
#   step 35             refuses to run before its data is installed, and says
#                       how to install it
#   steps 21 and 32     mark a CYP2D6 call when the depth check fails, and a
#                       caller that stops after writing its result leaves no
#                       call where step 36 reads it
# Files in CASE_WORK steer the hook: depth-status (the depth check's verdict,
# ok by default) and interrupt (the caller writes its result, then fails).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_sample "$G" sample1
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args="${*:2}"
opt() {
  local -a w
  read -r -d '' -a w <<<"$args" || true
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]//[\'\"]/}"; return 0; fi
  done
  return 1
}
put() {
  local h
  h=$(host_path "$1")
  mkdir -p "$(dirname "$h")"
  printf '%b' "$2" > "$h"
}
case "$args" in
  *"cyp2d6_depth_check.py check"*)
    put "$(opt --out)" "metric\\tvalue\\nstatus\\t$(cat "${CASE_WORK}/depth-status" 2>/dev/null || echo ok)\\nmessage\\tfake\\n" ;;
  *"-m cyrius"*)
    # bash -c SCRIPT _ INSTALL BAM PREFIX OUTDIR
    if [ -e "${CASE_WORK}/interrupt" ]; then
      called='*1/*2'
      put "${@: -1}${@: -2:1}.tsv" "Sample\\tGenotype\\tFilter\\nsample1\\t${called}\\tPASS\\n"
      exit 1
    fi
    put "${@: -1}${@: -2:1}.tsv" 'Sample\tGenotype\tFilter\nsample1\tNone\tNot_assigned_to_haplotypes\n' ;;
  *"summary_path = "*)
    # Step 32's summary, at the path its Python opens.
    suffix=""
    case "$args" in *"summary_path + '.partial'"*) suffix=.partial ;; esac
    put "/genome/sample1/pypgx/sample1_pypgx_summary.tsv${suffix}" \
      'Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\n'
    [ ! -e "${CASE_WORK}/interrupt" ] || exit 1 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"
run_lines() { grep "^run image=[^ ]*$1" "$FAKE_DOCKER_LOG" || true; }

# --- Cyrius: nothing before the install, then the hash-locked install --------------
run_expect 1 noinstall "${SCRIPTS}/21-cyrius.sh" sample1
output_has noinstall "Cyrius ${CYRIUS_VERSION} is not installed"
output_has noinstall 'scripts/setup.sh --cyrius '
[ -z "$(run_lines python)" ] || fail "step 21 started a container before Cyrius was installed"

: > "$FAKE_DOCKER_LOG"
run_expect 0 install "${SCRIPTS}/setup.sh" --cyrius "$G"
output_has install 'PolyForm Strict'
INSTALL=$(run_lines python)
[ "$(grep -c . <<<"$INSTALL")" -eq 1 ] || fail "setup.sh --cyrius ran $(grep -c . <<<"$INSTALL") python containers, expected 1"
grep -qE -- '--require-hashes --no-deps --only-binary :all: --target /genome/tools/cyrius-[0-9.]+\.part -r /lock\.txt' <<<"$INSTALL" \
  || fail "the install is not hash-locked from the lock file: ${INSTALL}"
grep -qE -- "-v ${REPO_ROOT}/scripts/cyrius-constraints\.txt:/lock\.txt:ro" <<<"$INSTALL" \
  || fail "the install does not read scripts/cyrius-constraints.txt: ${INSTALL}"
! grep -q -- '--network none' <<<"$INSTALL" || fail "the install has no network to download with"
[ -s "${G}/tools/cyrius-${CYRIUS_VERSION}/INSTALLED" ] || fail "no INSTALLED stamp"
grep -c -- '--hash=sha256:' "${REPO_ROOT}/scripts/cyrius-constraints.txt" >/dev/null || fail "the lock file has no hashes"
run_expect 0 again "${SCRIPTS}/setup.sh" --cyrius "$G"
output_has again 'already installed'

: > "$FAKE_DOCKER_LOG"
run_expect 0 cyrius "${SCRIPTS}/21-cyrius.sh" sample1
RUNS=$(grep '^run image=' "$FAKE_DOCKER_LOG")
[ "$(grep -c . <<<"$RUNS")" -ge 4 ] || fail "step 21 ran fewer containers than the depth check and Cyrius need: ${RUNS}"
if grep -v -- '--network none' <<<"$RUNS"; then fail "a container of step 21 had network"; fi
grep -q 'mosdepth' <<<"$RUNS" || fail "step 21 ran no depth check before Cyrius"
grep -q -- 'python3 -m cyrius' <<<"$RUNS" || fail "step 21 did not run Cyrius from the install"
if grep -q 'pip install' <<<"$RUNS"; then fail "step 21 still installs with pip"; fi

# --- step 36 and step 07 ----------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 consensus "${SCRIPTS}/36-pgx-consensus.sh" sample1
grep -q '^run image=[^ ]*python[^ ]* :: .*pgx_outside_calls\.py --calls /genome/sample1/pgx_consensus/sample1_outside_calls\.tsv' \
  "$FAKE_DOCKER_LOG" || fail "step 36 did not run bin/pgx_outside_calls.py"
grep -q -- '--cyrius /genome/sample1/cyrius/sample1_cyp2d6.tsv' "$FAKE_DOCKER_LOG" || fail "step 36 did not read Cyrius's TSV"

OUT="${G}/sample1/pgx_consensus/sample1_outside_calls.tsv"
: > "$OUT"
: > "$FAKE_DOCKER_LOG"
run_expect 0 pharmcat-empty "${SCRIPTS}/07-pharmacogenomics.sh" sample1
if grep 'pharmcat.jar' "$FAKE_DOCKER_LOG" | grep -q -- ' -po '; then fail "an empty outside-call file reached PharmCAT"; fi
printf 'HLA-A\t*02:01/*24:02\n' > "$OUT"
: > "$FAKE_DOCKER_LOG"
run_expect 0 pharmcat-po "${SCRIPTS}/07-pharmacogenomics.sh" sample1
grep 'pharmcat.jar' "$FAKE_DOCKER_LOG" | grep -q -- "-v ${OUT}:/outside_calls.tsv:ro .*-po /outside_calls.tsv" \
  || fail "PharmCAT did not get the outside calls: $(grep pharmcat.jar "$FAKE_DOCKER_LOG")"

# --- opt-in data not installed ---------------------------------------------------
run_expect 1 paralogs "${SCRIPTS}/35-paralogs.sh" sample1
output_has paralogs 'scripts/setup.sh --parascopy-data '
# --- a failed depth check marks the call; a caller that stops publishes nothing -------
CY_TSV="${G}/sample1/cyrius/sample1_cyp2d6.tsv"
echo unreliable > "${CASE_WORK}/depth-status"
run_expect 0 cyrius-unreliable "${SCRIPTS}/21-cyrius.sh" sample1
[ "$(awk -F'\t' 'NR == 2 {print $3}' "$CY_TSV" 2>/dev/null)" = CYP2D6_depth_unreliable ] \
  || fail "step 21 did not mark the Cyrius call after a failed depth check: $(cat "$CY_TSV" 2>/dev/null)"
[ ! -e "${G}/sample1/cyrius/partial" ] || fail "step 21 left its partial directory behind"
touch "${CASE_WORK}/interrupt"
run_rc cyrius-interrupted "${SCRIPTS}/21-cyrius.sh" sample1
[ "$RC" -ne 0 ] || fail "step 21 exited 0 when Cyrius failed"
[ ! -e "$CY_TSV" ] || fail "Cyrius stopped after writing its TSV and step 21 left that unchecked call for step 36: $(cat "$CY_TSV")"
rm -f "${CASE_WORK}/interrupt"

BUNDLE="${G}/reference/pypgx-bundle"
mkdir -p "$BUNDLE"
git -C "$BUNDLE" -c user.name=case -c user.email=case@example.invalid -c init.defaultBranch=main init -q
git -C "$BUNDLE" -c user.name=case -c user.email=case@example.invalid commit -q --allow-empty -m bundle
git -C "$BUNDLE" tag "$PYPGX_BUNDLE_VERSION"
SUMMARY="${G}/sample1/pypgx/sample1_pypgx_summary.tsv"
run_expect 0 pypgx-unreliable "${SCRIPTS}/32-pypgx.sh" sample1
[ "$(awk -F'\t' '$1 == "CYP2D6" {print $2}' "$SUMMARY" 2>/dev/null)" = Indeterminate ] \
  || fail "step 32 did not mark CYP2D6 Indeterminate after a failed depth check: $(cat "$SUMMARY" 2>/dev/null)"
[ ! -e "${SUMMARY}.partial" ] || fail "step 32 left its partial summary behind"
touch "${CASE_WORK}/interrupt"
run_rc pypgx-interrupted "${SCRIPTS}/32-pypgx.sh" sample1
[ "$RC" -ne 0 ] || fail "step 32 exited 0 when its summary failed"
[ ! -e "$SUMMARY" ] || fail "step 32 stopped after writing its summary and left that unchecked CYP2D6 call for step 36: $(cat "$SUMMARY")"
rm -f "${CASE_WORK}/interrupt" "${CASE_WORK}/depth-status"

echo "Cyrius installs hash-locked with network and runs without; step 07 passes only non-empty outside calls."
