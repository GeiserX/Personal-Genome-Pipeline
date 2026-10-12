#!/usr/bin/env bash
# run-all.sh after a failed pipeline records step 08 (HLA typing) from every
# process of that step: T1K_BUILD and HLA_TYPING, and with KIR=true KIR_BUILD
# and KIR_TYPING, the second T1K pass over the KIR genes. None of the other
# three is named after the step's tool (hla_typing), so each must be counted
# under it, or a KIR or index task that failed or did not finish leaves
# step 08 'ok'. Runs with TOOLS=hla_typing KIR=true and the pipeline exit 3:
#   1. every task of the four completed: step 08 is 'ok' (the control);
#   2. KIR_TYPING left running: 'not finished';
#   3. KIR_TYPING failed: 'failed', and the console names hla_typing;
#   4. T1K_BUILD failed, so HLA_TYPING got no task: 'failed'.
# A fake nextflow stands in; a wrapper in front of it writes the two files
# main.nf's completion handler would (run-all-partial-failure.sh has the rest).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
# Step 08's data where run-all.sh looks for it, and the KIR data KIR=true needs
for name in hla_dat gencode_genes; do
  # shellcheck source=../../scripts/lib/common.sh
  path=$(. "${REPO_ROOT}/scripts/lib/common.sh"; data_file "$name" || true)
  mkdir -p "$(dirname "$path")"; printf 'placeholder\n' > "$path"
done
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
mkdir -p "${G}/kir/IPD-KIR_${KIR_DB_RELEASE}"
printf 'placeholder\n' > "${G}/kir/IPD-KIR_${KIR_DB_RELEASE}/kir.dat"

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
# kir NAME WANT FAILED_TASKS ROWS: a run that stops with ROWS as its process table; step 08 must be WANT
kir() {
  local name=$1 want=$2 failed=$3 rows=$4 got
  rm -rf "${G}/sample1/logs"
  run_expect 3 "$name" env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 TOOLS=hla_typing KIR=true \
    PROCESS_STATUS="$rows" FAILED_TASKS="$failed" "${SCRIPTS}/run-all.sh" sample1 male
  output_has "$name" '^  08 HLA typing \(T1K\) +runs$'
  output_lacks "$name" 'KIR typing skipped'
  got=$(awk -F'\t' '$2 == "08" {print $3}' "$ST")
  [ "$got" = "$want" ] || miss "${name}: step 08 is '${got}', want '${want}'"
}

kir kir-done ok '' "$(row T1K_BUILD 1 0 0 0 0)$(row HLA_TYPING 1 0 0 0 0)$(row KIR_BUILD 1 0 0 0 0)$(row KIR_TYPING 1 0 0 0 0)"
docker_log_has '^nextflow :: .* --kir true --kir_dat ' "nextflow was not given --kir true and the KIR data"
kir kir-running 'not finished' '' "$(row T1K_BUILD 1 0 0 0 0)$(row HLA_TYPING 1 0 0 0 0)$(row KIR_BUILD 1 0 0 0 0)$(row KIR_TYPING 0 0 0 0 1)"
kir kir-failed failed 'BAM_ANALYSIS:KIR_TYPING\tfailed\n' \
  "$(row T1K_BUILD 1 0 0 0 0)$(row HLA_TYPING 1 0 0 0 0)$(row KIR_BUILD 1 0 0 0 0)$(row KIR_TYPING 0 0 1 0 0)"
grep -q 'hla_typing (KIR_TYPING) failed' "${CASE_WORK}/kir-failed.out" || miss "kir-failed: the console does not name hla_typing for KIR_TYPING"
kir t1k-failed failed 'BAM_ANALYSIS:T1K_BUILD\tfailed\n' \
  "$(row T1K_BUILD 0 0 1 0 0)$(row HLA_TYPING 0 0 0 0 0)$(row KIR_BUILD 0 0 0 0 0)$(row KIR_TYPING 0 0 0 0 0)"

[ "$MISSES" -eq 0 ] || fail "${MISSES} result(s) above are wrong"
echo "Step 08 is recorded from all four of its processes: a KIR or T1K index task that failed or did not finish makes it failed or not finished."
