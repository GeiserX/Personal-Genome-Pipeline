#!/usr/bin/env bash
# run-all.sh with HLA_QC=true passes --hla_qc true to the pipeline when
# hla_typing is selected, and records step 08 (HLA typing) from OptiType's two
# processes as well: HLA_QC_READS and HLA_QC are not named after the step's
# tool (hla_typing), so each must be counted under it, or an OptiType task that
# failed leaves step 08 'ok'. Runs with TOOLS=hla_typing and the pipeline exit 3:
#   1. HLA_QC=true, every task completed: --hla_qc true reaches nextflow and
#      step 08 is 'ok' (the control);
#   2. HLA_QC=true, HLA_QC failed: step 08 is 'failed';
#   3. HLA_QC=true, HLA_QC_READS failed: step 08 is 'failed';
#   4. HLA_QC unset, then HLA_QC=false: no --hla_qc on the command line.
# A fake nextflow stands in; a wrapper in front of it writes the two files
# main.nf's completion handler would (run-all-kir-step.sh is the KIR twin).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
# Step 08's data where run-all.sh looks for it
for name in hla_dat gencode_genes; do
  # shellcheck source=../../scripts/lib/common.sh
  path=$(. "${REPO_ROOT}/scripts/lib/common.sh"; data_file "$name" || true)
  mkdir -p "$(dirname "$path")"; printf 'placeholder\n' > "$path"
done

mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
[ -z "${PROCESS_STATUS:-}" ] || printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\n%b' "$PROCESS_STATUS" > process_status.tsv
[ -z "${FAILED_TASKS:-}" ] || printf 'process\tstatus\n%b' "$FAILED_TASKS" > failed_tasks.tsv
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"

# row PROCESS SUCCEEDED CACHED FAILED IGNORED UNFINISHED, as one escaped line for printf %b
row() { printf 'BAM_ANALYSIS:%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$@"; }
MISSES=0
miss() { echo "ASSERT FAIL: $*" >&2; MISSES=$((MISSES + 1)); }
# nf_line NAME: the nextflow command line of the run NAME (the last one logged)
nf_line() { awk '/^nextflow :: / {l = $0} END {print l}' "$FAKE_DOCKER_LOG"; }
# qc NAME HLA_QC WANT FAILED_TASKS ROWS: a run that stops with ROWS as its process table; step 08 must be WANT
qc() {
  local name=$1 hq=$2 want=$3 failed=$4 rows=$5 got
  rm -rf "${G}/sample1/logs"
  run_expect 3 "$name" env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 TOOLS=hla_typing ${hq:+HLA_QC=$hq} \
    PROCESS_STATUS="$rows" FAILED_TASKS="$failed" "${SCRIPTS}/run-all.sh" sample1 male
  output_has "$name" '^  08 HLA typing \(T1K\) +runs$'
  got=$(awk -F'\t' '$2 == "08" {print $3}' "$ST")
  [ "$got" = "$want" ] || miss "${name}: step 08 is '${got}', want '${want}'"
}
T1K="$(row T1K_BUILD 1 0 0 0 0)$(row HLA_TYPING 1 0 0 0 0)"

qc qc-done true ok '' "${T1K}$(row HLA_QC_READS 1 0 0 0 0)$(row HLA_QC 1 0 0 0 0)"
nf_line | grep -q -- ' --hla_qc true' || miss "qc-done: nextflow was not given --hla_qc true"
output_has qc-done 'HLA_QC=true: OptiType checks'
qc qc-failed true failed 'BAM_ANALYSIS:HLA_QC\tfailed\n' "${T1K}$(row HLA_QC_READS 1 0 0 0 0)$(row HLA_QC 0 0 1 0 0)"
grep -q 'hla_typing (HLA_QC) failed' "${CASE_WORK}/qc-failed.out" || miss "qc-failed: the console does not name hla_typing for HLA_QC"
qc qc-reads-failed true failed 'BAM_ANALYSIS:HLA_QC_READS\tfailed\n' "${T1K}$(row HLA_QC_READS 0 0 1 0 0)$(row HLA_QC 0 0 0 0 0)"
for hq in '' false; do
  qc "qc-off${hq:+-$hq}" "$hq" ok '' "$T1K"
  if nf_line | grep -q -- '--hla_qc'; then miss "qc-off${hq:+-$hq}: nextflow was given --hla_qc without HLA_QC=true"; fi
  output_lacks "qc-off${hq:+-$hq}" 'OptiType checks'
done

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "HLA_QC=true reaches the pipeline as --hla_qc true, and step 08 counts OptiType's two processes."
