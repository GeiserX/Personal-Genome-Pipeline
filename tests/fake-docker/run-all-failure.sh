#!/usr/bin/env bash
# run-all.sh stops on a failed validation (exit 1, nextflow not started). After
# the pipeline it reads <sample>/nextflow/failed_tasks.tsv, which main.nf's
# completion handler writes when a task failed:
#   - a report-only tool that failed and was ignored (the pipeline went on and
#     exited 0): its step is 'failed' in run_status.tsv, the others 'ok', both
#     reports are written, the console names the tool and the TOOLS list
#     without it, and run-all.sh exits 1;
#   - a failed pipeline (exit 3): the exit code is passed on, the reports are
#     still written, no step is recorded ok, the failed tool's step is
#     'failed', and the error names the process and the TOOLS entry to drop;
#   - a later clean run does not read the earlier run's file.
# A fake nextflow (logs its arguments) stands in; a wrapper in front of it
# writes the file the completion handler would.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_clinvar "$G"
use_output_hook   # the reports at the end read the files their containers write
seed_sample "$G" sample1
calls() { grep -c '^nextflow :: ' "$FAKE_DOCKER_LOG" || true; }
ST="${G}/sample1/logs/run_status.tsv"

# No reference: validate-setup.sh fails
run_expect 1 novalidate "${SCRIPTS}/run-all.sh" sample1 male
output_has novalidate '=== Summary ==='
output_has novalidate 'ERROR: setup validation failed'
[ "$(calls)" -eq 0 ] || fail "nextflow started after a failed validation"

# The same, with the validation skipped: nextflow starts
run_expect 0 skipvalidate env SKIP_VALIDATION=true "${SCRIPTS}/run-all.sh" sample1 male
output_lacks skipvalidate '=== Summary ==='
[ "$(calls)" -eq 1 ] || fail "nextflow calls: $(calls), expected 1"

# A nextflow that leaves FAILED_TASKS (process<TAB>status lines) where the
# completion handler writes them, in its launch directory
mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
[ -z "${FAILED_TASKS:-}" ] || printf 'process\tstatus\n%b' "$FAILED_TASKS" > failed_tasks.tsv
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"
TOOLS_ALL='pharmcat,cpic,roh,mito_haplogroup,mosdepth,telomere_hunter,mito_variants,manta,delly,duphold,survivor_merge,multiqc,clinvar'

# A report-only tool failed and was ignored; the pipeline exited 0
rm -rf "${G}/sample1/logs"
run_expect 1 ignored env SKIP_VALIDATION=true PATH="$NFPATH" FAILED_TASKS='BAM_ANALYSIS:TELOMERE_HUNTER\tignored\n' \
  "${SCRIPTS}/run-all.sh" sample1 male
grep -q $'^step\t10\tfailed$' "$ST" || fail "step 10 (TelomereHunter) is not 'failed': $(cat "$ST")"
for s in 07 16 16b 19 28; do grep -q $'^step\t'"${s}"$'\tok$' "$ST" || fail "step ${s} is not ok: $(cat "$ST")"; done
[ "$(grep -c $'\tfailed$' "$ST")" -eq 1 ] || fail "more than one step is failed: $(cat "$ST")"
for r in sample1_report.html sample1_report.txt; do [ -f "${G}/sample1/${r}" ] || fail "${r} was not written"; done
output_has ignored 'telomere_hunter \(TELOMERE_HUNTER\) failed'
output_has ignored "TOOLS=${TOOLS_ALL/telomere_hunter,/}"

# The pipeline failed in a tool's task (exit 3)
rm -rf "${G}/sample1/logs" "${G}/sample1/sample1_report.html" "${G}/sample1/sample1_report.txt"
run_expect 3 nffail env SKIP_VALIDATION=true PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 FAILED_TASKS='SV:DELLY\tfailed\n' \
  "${SCRIPTS}/run-all.sh" sample1 male
output_has nffail 'ERROR: the pipeline failed \(exit 3\)'
output_has nffail 'nextflow/\.nextflow\.log'
output_has nffail 'delly \(DELLY\) failed'
output_has nffail "TOOLS=${TOOLS_ALL/delly,/}"
for r in sample1_report.html sample1_report.txt; do [ -f "${G}/sample1/${r}" ] || fail "${r} was not written after the pipeline failed"; done
if grep -q $'\tok$' "$ST"; then fail "a step is recorded ok after the pipeline failed: $(cat "$ST")"; fi
grep -q $'^step\t19\tfailed$' "$ST" || fail "step 19 (Delly) is not 'failed': $(cat "$ST")"
grep -q $'^step\t07\tnot finished$' "$ST" || fail "step 07 is not 'not finished': $(cat "$ST")"

# The pipeline failed before any task (no file): the exit code is passed on
rm -rf "${G}/sample1/logs"
run_expect 3 nffail-early env SKIP_VALIDATION=true FAKE_NEXTFLOW_RC=3 "${SCRIPTS}/run-all.sh" sample1 male
output_has nffail-early 'ERROR: the pipeline failed \(exit 3\)'
output_lacks nffail-early 'DELLY'
if grep -qE $'\t(ok|failed)$' "$ST"; then fail "a step is recorded ok or failed: $(cat "$ST")"; fi

# A clean run after them: nothing is failed
run_expect 0 clean env SKIP_VALIDATION=true "${SCRIPTS}/run-all.sh" sample1 male
if grep -q $'\tfailed$' "$ST"; then fail "the clean run read an earlier run's failed tasks: $(cat "$ST")"; fi
echo "A failed validation stops run-all.sh before nextflow; a failed or ignored task is named, its step is failed, and the reports are written."
