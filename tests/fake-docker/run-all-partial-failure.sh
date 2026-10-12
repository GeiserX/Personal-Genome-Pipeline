#!/usr/bin/env bash
# run-all.sh after a pipeline that one tool stopped (exit 3). main.nf's
# completion handler leaves <sample>/nextflow/process_status.tsv: one row per
# process with its tasks succeeded, cached, failed, ignored and unfinished.
# run-all.sh must record each step from those rows, not from the exit code:
#   1. a -resume run where PRS fails, some tools' tasks completed and others
#      came from the cache: completed and cached steps are 'ok', PRS and the
#      ancestry step inside it 'failed', a report-only tool whose task was
#      ignored 'failed', and only the steps the stop cut short 'not finished'
#      (a process that got no task, or one left running);
#   2. the reports of that run: a cached step's file from an earlier run is
#      shown as current in both, and a step that did not finish keeps the STALE
#      mark on its old file;
#   3. the same failure on a first run (nothing cached): the same statuses;
#   4. no process_status.tsv (the pipeline died before its handler ran): no
#      step is 'ok';
#   5. every tool run-all.sh can select has a process named after it in
#      modules/local, which is how a row finds its step.
# A fake nextflow stands in; a wrapper in front of it writes the two files the
# completion handler would. The reports are the real ones: the hook runs
# bin/render_report.py with the host's python3.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
S="${G}/sample1"
ST="${S}/logs/run_status.tsv"
seed_clinvar "$G"
seed_sample "$G" sample1
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
PANEL="${G}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst"
mkdir -p "${G}/prs_scores" "$(dirname "$PANEL")"
printf 'placeholder\n' | gzip -c > "${G}/prs_scores/PGS000018.txt.gz"
printf 'placeholder\n' > "$PANEL"
printf 'chr1\t1\tA\tG\n' > "${PANEL%.tar.zst}_GRCh38_sites.tsv"

# Results of an earlier run, 30 days old: ROH and Manta (cached in this run)
# and Delly (cut short in this run).
mkdir -p "${S}/roh" "${S}/manta/results/variants" "${S}/delly"
printf 'RG\tsample1\tchr1\t1000000\t3000000\t2000000\t500\t40.0\n' > "${S}/roh/sample1_roh.txt"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > "${S}/manta/results/variants/diploidSV.vcf.gz"
cp "${S}/manta/results/variants/diploidSV.vcf.gz" "${S}/delly/sample1_sv.vcf.gz"
touch -d '30 days ago' "${S}/roh/sample1_roh.txt" "${S}/manta/results/variants/diploidSV.vcf.gz" "${S}/delly/sample1_sv.vcf.gz"

use_output_hook
cat > "${CASE_WORK}/render-hook" <<'HOOK'
#!/usr/bin/env bash
# The real renderer, on the host paths of every /... argument.
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
case "${*:2}" in
  "python3 /pgp-bin/render_report.py"*)
    shift
    mapped=()
    for a in "$@"; do
      if [[ "$a" == /* ]] && h=$(host_path "$a"); then mapped+=("$h"); else mapped+=("$a"); fi
    done
    exec "${mapped[@]}" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/render-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/render-hook"

# A nextflow that leaves PROCESS_STATUS and FAILED_TASKS (files of rows) where
# the completion handler writes them, in its launch directory.
mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
[ -z "${PROCESS_STATUS:-}" ] || { printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\n'; cat "$PROCESS_STATUS"; } > process_status.tsv
[ -z "${FAILED_TASKS:-}" ] || { printf 'process\tstatus\n'; cat "$FAILED_TASKS"; } > failed_tasks.tsv
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"

# process, succeeded, cached, failed (not retried, not ignored), ignored, unfinished
row() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$@"; }
{
  row UPSTREAM:INDEXCOV 0 1 0 0 0
  row VCF_PRECHECK 0 1 0 0 0
  row PGX:CLINVAR_SCREEN 0 1 0 0 0
  row PGX:PHARMCAT_PREPROCESS 0 1 0 0 0
  row PGX:PHARMCAT 1 0 0 0 0
  row PGX:CPIC_LOOKUP 1 0 0 0 0
  row CLINICAL:ROH 0 1 0 0 0
  row CLINICAL:MITO_EXTRACT_CHRM 1 0 0 0 0
  row CLINICAL:MITO_HAPLOGROUP 1 0 0 0 0
  row CLINICAL:HAPLOCHECK 0 0 0 0 0
  row BAM_ANALYSIS:MOSDEPTH 1 0 0 0 0
  row BAM_ANALYSIS:TELOMERE_HUNTER 0 0 0 1 0
  row BAM_ANALYSIS:MITO_VARIANTS 0 1 0 0 0
  row SV:MANTA 0 1 0 0 0
  row SV:DELLY 1 0 0 0 0
  row SV:DELLY_BCF2VCF 0 0 0 0 0
  row SV:DUPHOLD 0 0 0 0 1
  row SV:DUPHOLD_FILTER 0 0 0 0 0
  row SV:SURVIVOR_PREP 1 0 0 0 0
  row SV:SURVIVOR_MERGE 0 0 0 0 0
  row SV:SURVIVOR_SORT 0 0 0 0 0
  row CLINICAL:PRS_PREPARE 0 1 0 0 0
  row CLINICAL:PRS_SCORE_SITES 0 1 0 0 0
  row CLINICAL:PRS 0 0 1 0 0
  row CLINICAL:PRS_SUMMARY 0 0 0 0 0
  row REPORTING:MULTIQC 0 0 0 0 0
} > "${CASE_WORK}/resume.tsv"
printf 'BAM_ANALYSIS:TELOMERE_HUNTER\tignored\nCLINICAL:PRS\tfailed\n' > "${CASE_WORK}/failed.tsv"

# A wrong result is counted, not fatal, so one run lists every step and report line that is wrong.
MISSES=0
miss() { echo "ASSERT FAIL: $*" >&2; MISSES=$((MISSES + 1)); }
# want NAME STEP RESULT: the step's line in run_status.tsv
want() { grep -q $'^step\t'"$2"$'\t'"$3"'$' "$ST" || miss "$1: step $2 is not '$3' but '$(awk -F'\t' -v s="$2" '$2 == s {print $3}' "$ST")'"; }
statuses() {
  local s
  for s in 16 06 07 27 11 12 16b 20 04; do want "$1" "$s" ok; done   # completed or cached
  for s in 25 26 10; do want "$1" "$s" failed; done                   # PRS, the ancestry step inside it, the ignored tool
  for s in 19 15 22 28; do want "$1" "$s" 'not finished'; done        # no task, or one left running
  [ "$(grep -c $'\tok$' "$ST")" -eq 9 ] || miss "$1: $(grep -c $'\tok$' "$ST") steps are ok, want 9"
  [ "$(grep -c $'\tnot finished$' "$ST")" -eq 4 ] || miss "$1: $(grep -c $'\tnot finished$' "$ST") steps are not finished, want 4"
}

# --- 1. a -resume run that PRS stops ------------------------------------------------------
run_expect 3 resume env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 PROCESS_STATUS="${CASE_WORK}/resume.tsv" \
  FAILED_TASKS="${CASE_WORK}/failed.tsv" "${SCRIPTS}/run-all.sh" sample1 male
output_has resume '^  26 Ancestry \(pgsc_calc\) +runs$'
output_has resume 'prs \(PRS\) failed: '
output_has resume 'telomere_hunter \(TELOMERE_HUNTER\) failed, skipped'
statuses resume

# --- 2. the reports of that run -----------------------------------------------------------
TXT="${S}/sample1_report.txt"
HTML="${S}/sample1_report.html"
[ -s "$TXT" ] || fail "no text report: $(cat "${S}/logs/generate-report.log" 2>/dev/null)"
[ -s "$HTML" ] || fail "no HTML report: $(cat "${S}/logs/24-html-report.log" 2>/dev/null)"
grep -q '^## Runs of Homozygosity' "$TXT" || fail "the text report has no ROH section, so its stale check would prove nothing"
grep -q '^## Structural Variants' "$TXT" || fail "the text report has no SV section, so its stale check would prove nothing"
for r in "$TXT" "$HTML"; do
  for s in 11 04; do
    if grep -q "STALE: from an earlier run.*step ${s} was" "$r"; then
      miss "${r##*/} marks step ${s} stale, but its task came from the cache in this run: $(grep -o -m1 "STALE: from an earlier run[^<]*step ${s} was[^<]*" "$r")"
    fi
  done
  grep -q 'STALE: from an earlier run.*step 19 was not finished' "$r" \
    || miss "${r##*/} does not mark step 19's old result stale, but Delly did not finish in this run"
done

# --- 3. the same failure on a first run: nothing cached ------------------------------------
awk -F'\t' -v OFS='\t' '{$2 += $3; $3 = 0; print}' "${CASE_WORK}/resume.tsv" > "${CASE_WORK}/first.tsv"
if awk -F'\t' '$3 != 0 {f = 1} END {exit !f}' "${CASE_WORK}/first.tsv"; then fail "first.tsv still has a cached task"; fi
rm -rf "${S}/logs"
run_expect 3 first env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 PROCESS_STATUS="${CASE_WORK}/first.tsv" \
  FAILED_TASKS="${CASE_WORK}/failed.tsv" "${SCRIPTS}/run-all.sh" sample1 male
statuses first

# --- 4. no process table: nothing is ok ---------------------------------------------------
rm -rf "${S}/logs"
run_expect 3 notable env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 FAILED_TASKS="${CASE_WORK}/failed.tsv" "${SCRIPTS}/run-all.sh" sample1 male
if grep -q $'\tok$' "$ST"; then miss "notable: a step is ok with no process_status.tsv: $(grep $'\tok$' "$ST" | tr '\n' '|')"; fi
want notable 25 failed
want notable 11 'not finished'
[ ! -e "${S}/nextflow/process_status.tsv" ] || miss "notable: run-all.sh kept an earlier run's process_status.tsv"

# --- 5. every selectable tool has a process named after it --------------------------------
run_expect 2 known env TOOLS=no_such_tool "${SCRIPTS}/run-all.sh" sample1 male
KNOWN=$(sed -n "s/^ERROR: unknown step 'no_such_tool' in TOOLS\. Known: //p" "${CASE_WORK}/known.out")
PROCS=$(sed -n 's/^process \([A-Z0-9_]*\).*/\1/p' "${REPO_ROOT}"/modules/local/*/main.nf)
[ "$(wc -w <<<"$KNOWN")" -ge 25 ] || fail "run-all.sh names too few tools to check: '${KNOWN}'"
[ "$(wc -l <<<"$PROCS")" -ge 40 ] || fail "too few processes found in modules/local to check"
for t in $KNOWN; do
  case "$t" in
    ancestry) continue ;;   # step 26 is computed inside the PRS run
    y_haplogroup) re='^Y_' ;;
    survivor_merge) re='^SURVIVOR_' ;;
    *) re="^${t^^}(_|\$)" ;;
  esac
  grep -qE "$re" <<<"$PROCS" || miss "tool '${t}' has no process named ${t^^} or ${t^^}_<part> in modules/local: run-all.sh could never record its step ok after a failed run"
done
if grep -qE '^NO_SUCH_TOOL(_|$)' <<<"$PROCS"; then fail "control: a made-up tool name found a process"; fi

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "After a failed pipeline each step is recorded from its own processes: completed and cached steps are ok and current in the reports, and only a step that did not finish is stale."
