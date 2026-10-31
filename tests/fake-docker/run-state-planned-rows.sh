#!/usr/bin/env bash
# run-all.sh writes the run's plan before the pipeline starts, so
# bin/run_state.py can show the run as it goes:
#   1. while nextflow runs, logs/run_status.tsv has a 'planned' row for every
#      step that runs (a GRIDSS=true script step too), next to the skipped
#      rows, and no 'meta finished_utc'; bin/run_state.py reads them as
#      waiting steps and a running run;
#   2. after a clean run each step has one row, its result (the planned row
#      gave way), and 'meta finished_utc' closes the file; run_state.py then
#      says the run finished;
#   3. after a failed pipeline (exit 3) the script step run-all.sh no longer
#      starts is 'not finished', not left 'planned', and the end time is there.
# A fake nextflow stands in; a wrapper in front of it copies the status file
# and runs bin/run_state.py with the host's python3 while "the pipeline runs".
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
S="${G}/sample1"
ST="${S}/logs/run_status.tsv"
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
seed_sample "$G" sample1

mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
# The launch directory is <sample>/nextflow: the status file is ../logs/run_status.tsv.
cp ../logs/run_status.tsv "${CASE_WORK}/during.tsv"
python3 "${REPO_ROOT}/bin/run_state.py" --sample sample1 --sample-dir .. --out "${CASE_WORK}/during.json" > /dev/null
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"
DURING="${CASE_WORK}/during.tsv"
rows() { grep -c $'^step\t'"$1"$'\t' "$2" || true; }
# jq-free reads of run_state.json
rs() { python3 -c 'import json, sys; d = json.load(open(sys.argv[1])); s = {x["step"]: x for x in d["steps"]}
print(d["run"]["state"] if sys.argv[2] == "run" else s[sys.argv[2]]["state"])' "$@"; }

# --- 1 and 2: a clean run with a script step ----------------------------------------------
# GRIDSS fails on the fake docker, so run-all.sh exits 1; the case reads the rows, not the exit code.
run_rc clean env PATH="$NFPATH" GRIDSS=true "${SCRIPTS}/run-all.sh" sample1 male
[ "$RC" -le 1 ] || fail "clean exited ${RC}, expected 0 or 1"
[ -s "$DURING" ] || fail "the wrapper did not see the status file"
for s in 06 07 16 04 04b; do
  grep -q $'^step\t'"${s}"$'\tplanned$' "$DURING" || fail "no 'planned' row for step ${s} while the pipeline ran: $(cat "$DURING")"
done
grep -q $'^step\t13\tskipped (data not installed: ' "$DURING" || fail "the skipped row of step 13 is missing: $(cat "$DURING")"
if grep -q $'^step\t13\tplanned$' "$DURING"; then fail "a skipped step is planned: $(cat "$DURING")"; fi
if grep -q $'^meta\tfinished_utc\t' "$DURING"; then fail "finished_utc is written before the run ended"; fi
# One planned row per --tools name, plus step 16 (indexcov, no --tools name) and 04b
tools=$(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1 | sed -n 's/.* --tools \([^ ]*\) .*/\1/p' | tr -d '\\')
want=$(( $(tr ',' '\n' <<<"$tools" | grep -c .) + 2 ))
n=$(grep -c $'\tplanned$' "$DURING"); [ "$want" -gt 10 ] && [ "$n" -eq "$want" ] || fail "${n} planned rows while the pipeline ran, want ${want} (--tools ${tools}, 16, 04b): $(cat "$DURING")"
[ "$(rs "${CASE_WORK}/during.json" run)" = running ] || fail "run_state.py: the run is not running while nextflow runs"
[ "$(rs "${CASE_WORK}/during.json" 07)" = waiting ] || fail "run_state.py: step 07 is not waiting while nextflow runs"
[ "$(rs "${CASE_WORK}/during.json" 13)" = skipped ] || fail "run_state.py: step 13 is not skipped"

if grep -q $'\tplanned$' "$ST"; then fail "a planned row is left after the run: $(cat "$ST")"; fi
for s in 06 07 16 04 04b; do [ "$(rows "$s" "$ST")" -eq 1 ] || fail "step ${s} has $(rows "$s" "$ST") rows after the run, want 1: $(cat "$ST")"; done
grep -q $'^step\t07\tok$' "$ST" || fail "step 07 is not ok: $(cat "$ST")"
grep -qE $'^step\t04b\t(ok|failed)$' "$ST" || fail "step 04b has no result: $(cat "$ST")"
grep -qE $'^meta\tfinished_utc\t[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]{8} UTC$' "$ST" || fail "no finished_utc row: $(cat "$ST")"
[ "$(tail -n 1 "$ST" | cut -f2)" = finished_utc ] || fail "finished_utc is not the last row: $(cat "$ST")"
python3 "${REPO_ROOT}/bin/run_state.py" --sample sample1 --sample-dir "$S" --out "${CASE_WORK}/after.json" > /dev/null
[ "$(rs "${CASE_WORK}/after.json" run)" = finished ] || fail "run_state.py: the run is not finished after run-all.sh ended"
[ "$(rs "${CASE_WORK}/after.json" 07)" = done ] || fail "run_state.py: step 07 is not done after the run"

# --- 3. a failed pipeline: the script step is not finished -----------------------------------
rm -f "$DURING"
run_expect 3 failed env PATH="$NFPATH" GRIDSS=true FAKE_NEXTFLOW_RC=3 "${SCRIPTS}/run-all.sh" sample1 male
grep -q $'^step\t04b\tplanned$' "$DURING" || fail "no 'planned' row for step 04b in the failed run: $(cat "$DURING")"
grep -q $'^step\t04b\tnot finished$' "$ST" || fail "step 04b is not 'not finished' after the failed pipeline: $(cat "$ST")"
[ "$(rows 04b "$ST")" -eq 1 ] || fail "step 04b has $(rows 04b "$ST") rows, want 1"
grep -q $'^meta\tfinished_utc\t' "$ST" || fail "no finished_utc after the failed run"
echo "run-all.sh wrote the plan before the pipeline, one result per step after it, and the end time."
