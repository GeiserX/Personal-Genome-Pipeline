#!/usr/bin/env bash
# Stub runs of main.nf (profile test_all, no containers) on the stub VCF+BAM row
# under short or bad sample ids, with a failure planted through beforeScript:
#   1. id S1 (CPSR 2.3.2 takes 3 to 40 characters), TELOMERE_HUNTER failing:
#      the run ends with exit 0, CPSR's files carry S1's name (its task runs
#      bin/cpsr_sample_id), TELOMERE_HUNTER is ignored and named at the end
#      and in failed_tasks.tsv, and the report is written;
#   2. the same with HLA_TYPING failing instead: HLA feeds PharmCAT, so its
#      failure still stops the run, and failed_tasks.tsv names it as failed;
#   3. id '-x': the samplesheet check stops the run before any task, saying
#      why (tools read '-x' as an option).
# The stubs only touch files; nothing here needs a database or an image.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

ST="${REPO}/assets/stub"
# stub_run NAME ID [ARGS...]: one stub run in its own launch directory; sets RC, LOG and DIR
stub_run() {
  local name=$1 id=$2
  shift 2
  DIR="${CASE_TMP}/${name}"
  rm -rf "$DIR"
  mkdir -p "$DIR"
  printf 'sample,bam,bam_index,vcf,vcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$id" \
    "${ST}/sample.bam" "${ST}/sample.bam.bai" "${ST}/sample.vcf.gz" "${ST}/sample.vcf.gz.tbi" > "${DIR}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile test_all -stub (${name}, sample '${id}') $*"
  ( cd "$DIR" && nextflow run "${REPO}/main.nf" -profile test_all -stub -ansi-log false \
      -work-dir "${DIR}/work" --input "${DIR}/samplesheet.csv" --outdir "${DIR}/out" \
      --tools cpsr,telomere_hunter,hla_typing,pharmcat,cpic,html_report "$@" ) > "${DIR}/run.log" 2>&1
  RC=$?
  LOG=$(cat "${DIR}/run.log")
  echo "$LOG"
  cp "${DIR}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
}
# plant PROCESS: a config whose beforeScript makes every task of PROCESS exit 1
plant() {
  printf "process { withName: '%s' { beforeScript = 'echo planted failure >&2; exit 1' } }\n" "$1" > "${CASE_TMP}/plant-$1.config"
  echo "${CASE_TMP}/plant-$1.config"
}
# tasks_file_has PROCESS STATUS: a row of the launch directory's failed_tasks.tsv
tasks_file_has() { awk -F'\t' -v p="$1" -v s="$2" '$1 ~ ("(^|:)" p "$") && $2 == s {f = 1} END {exit !f}' "${DIR}/failed_tasks.tsv" 2>/dev/null; }

# --- 1. S1, a report-only tool failing -------------------------------------------------
stub_run leaf S1 -c "$(plant TELOMERE_HUNTER)"
check_eq "S1 with TELOMERE_HUNTER failing: the run exits 0" "$RC" 0
O="${DIR}/out/S1"
check "CPSR's HTML is named after the sample" test -e "${O}/cpsr/S1.cpsr.grch38.html"
check "CPSR's classification is named after the sample" test -e "${O}/cpsr/S1.cpsr.grch38.classification.tsv.gz"
check "no CPSR file is named after the padded id" test -z "$(find "${O}/cpsr" -name 'S1_cpsr*' 2>/dev/null)"
CPSR_SH=$(find "${DIR}/work" -name .command.sh -exec grep -l 'cpsr_sample_id' {} + 2>/dev/null | head -n 1)
check "the CPSR task ran bin/cpsr_sample_id" test -n "$CPSR_SH"
check "the planted failure ran (TelomereHunter's task failed)" has 'TELOMERE_HUNTER.*(Error is ignored|terminated with an error)' "$LOG"
check "the end of the run names the ignored task" has 'failed and (was|were) skipped: .*TELOMERE_HUNTER' "$LOG"
check "it does not claim a clean success" lacks 'Pipeline completed successfully' "$LOG"
check "failed_tasks.tsv lists TELOMERE_HUNTER as ignored" tasks_file_has TELOMERE_HUNTER ignored
check "the report was written" test -e "${O}/S1_report.html"

# --- 2. the same with a tool that feeds another ------------------------------------------
stub_run nonleaf S1 -c "$(plant HLA_TYPING)"
check "S1 with HLA_TYPING failing: the run stops (exit ${RC})" test "$RC" -ne 0
check "failed_tasks.tsv lists HLA_TYPING as failed" tasks_file_has HLA_TYPING failed
check "the end of the run names the failed process" has 'Pipeline failed.*HLA_TYPING' "$LOG"

# --- 3. an id that starts with '-' ------------------------------------------------------
stub_run dash -x
check "sample '-x': the run stops (exit ${RC})" test "$RC" -ne 0
check "the message says why" has "Sample name '-x' starts with '-'" "$LOG"
check_eq "no task ran" "$(find "${DIR}/work" -name .command.run 2>/dev/null | wc -l | tr -d ' ')" 0

finish
