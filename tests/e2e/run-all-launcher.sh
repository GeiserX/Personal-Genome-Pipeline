#!/usr/bin/env bash
# scripts/run-all.sh, the launcher, from the fixture's FASTQ pair as sample
# HG002P, in a GENOME_DIR of its own that holds only what case
# nextflow-from-fastq-2-nextflow gives the pipeline directly: the reference
# without its minimap2 index, ClinVar, the HLA data and the PRS score files.
# TOOLS is that case's list without vcfanno (run-all.sh runs it only with VEP)
# and html_report (run-all.sh renders the full report with step 24 instead).
#   1. It publishes the same files as the direct `nextflow run`, and the seven
#      items of scripts/ci/parity-diff.sh agree with the bash steps' outputs.
#   2. A second run finds every task in Nextflow's cache (-resume): no task
#      runs again, and it ends in minutes.
#   3. The reports and logs/run_status.tsv are written.
#   4. It runs as the guide documents it, without THREADS or --max_memory:
#      the launcher passes the runner's CPU count and RAM, and Nextflow starts
#      the 8-CPU and 32 GB tasks capped to them instead of refusing them.
#   5. Its work directory and GENOME_DIR share a filesystem, so it publishes by
#      hard link: the BAM is one file with two names (work/ and aligned/), and
#      no published file is a symbolic link into work/.
# Both comparisons get a negative control in the same run: a planted missing
# file must be reported, and the first run's trace (tasks COMPLETED) must fail
# the all-cached check.
# SKIP_VALIDATION=true: the fixture's reference is a slice that validate-setup.sh
# rejects (case 10); the fake-docker cases check that run-all.sh calls it.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
P=HG002P
G2="${E2E_WORK}/genome-runall"
DIRECT="${GENOME_DIR}/nf-fastq/${P}"
check "case nextflow-from-fastq-2-nextflow published its outputs" test -d "$DIRECT"
rm -rf "$G2"
mkdir -p "${G2}/reference" "${G2}/clinvar" "${G2}/prs_scores" "${G2}/${P}/fastq"
for e in fasta fasta.fai dict; do ln -f "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.${e}" "${G2}/reference/"; done
ln -f "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" "${G2}/clinvar/"
for r in R1 R2; do ln -f "${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz" "${G2}/${P}/fastq/${P}_${r}.fastq.gz"; done
ln -f "${E2E_WORK}/parity-pgs/"*.txt.gz "${G2}/prs_scores/"
for d in hla_dat gencode_genes; do
  from=$(GENOME_DIR="$GENOME_DIR" bash -c '. "$1/scripts/lib/common.sh" && data_file "$2"' _ "$REPO" "$d")
  to=$(GENOME_DIR="$G2" bash -c '. "$1/scripts/lib/common.sh" && data_file "$2"; true' _ "$REPO" "$d")
  mkdir -p "$(dirname "$to")" && ln -f "$from" "$to"
done
INTERVALS=$(awk '{printf "%s%s:%d-%d", (NR > 1 ? " " : ""), $1, $2 + 1, $3}' "${FIXTURE_DIR}/regions.bed")
TOOLS=pharmcat,cpic,roh,prs,mito_haplogroup,telomere_hunter,mosdepth,mito_variants,multiqc,clinvar,hla_typing

launch() {  # launch NAME: run-all.sh, output to the case log and logs/<case>.<NAME>.log; exit code in RC, seconds in SECS
  local t0=$SECONDS
  echo "+ scripts/run-all.sh ${P} male --sex_check warn, THREADS unset (${1})"
  env -u THREADS GENOME_DIR="$G2" SKIP_VALIDATION=true TOOLS="$TOOLS" INTERVALS="$INTERVALS" \
    "${REPO}/scripts/run-all.sh" "$P" male --sex_check warn 2>&1 | tee "${E2E_WORK}/logs/${CASE_NAME}.${1}.log"
  RC=${PIPESTATUS[0]} SECS=$((SECONDS - t0))
  echo "+ exit ${RC} after ${SECS}s"
}
latest_trace() { find "${G2}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}'; }
# status_count TRACE STATUS: tasks with that status
status_count() {
  awk -F'\t' -v s="$2" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next} $c["status"] == s {k++} END {print k + 0}' "$1" 2>/dev/null || echo 0
}
all_cached() { [ "$(status_count "$1" CACHED)" -ge 10 ] && [ "$(status_count "$1" COMPLETED)" -eq 0 ] && [ "$(status_count "$1" FAILED)" -eq 0 ]; }
not_all_cached() { ! all_cached "$1"; }
# published DIR: the files under DIR, relative, without the launcher's own (inputs, launch directory,
# logs, manifest, bash reports) and the two outputs the direct case makes and run-all.sh does not.
published() {
  (cd "$1" && find . -type f | sed 's|^\./||' | LC_ALL=C sort) \
    | grep -vE "^(fastq|nextflow|logs|vep)/|^run_manifest\.tsv$|^summary\.json$|^${P}_report\.(html|txt)$"
}

# --- first run ------------------------------------------------------------------------
launch first
check_eq "run-all.sh exits 0" "$RC" 0
FIRST_TRACE=$(latest_trace)
check "the first run wrote a trace" test -s "${FIRST_TRACE:-/dev/null}"
check_ge "tasks completed in the first run" "$(status_count "${FIRST_TRACE:-/dev/null}" COMPLETED)" 10
check "the samplesheet is a FASTQ row" has "^${P},${G2}/${P}/fastq/${P}_R1.fastq.gz,${G2}/${P}/fastq/${P}_R2.fastq.gz,,,,,male\$" \
  "$(sed -n 2p "${G2}/${P}/nextflow/samplesheet.csv")"

# --- 4. the caps are the runner's ---------------------------------------------------------
NFCMD=$(grep -m1 -F '$> nextflow run' "${G2}/${P}/nextflow/.nextflow.log" 2>/dev/null || true)
echo "+ ${NFCMD#*\$> }"
check "nextflow got --max_cpus with the runner's CPU count" has " --max_cpus $(getconf _NPROCESSORS_ONLN) " "${NFCMD} "
check "nextflow got --max_memory with the runner's RAM" has " --max_memory $(awk '/^MemTotal:/ {print int($2 / 1048576)}' /proc/meminfo)\.GB " "${NFCMD} "

# --- 5. published by hard link ------------------------------------------------------------
check "nextflow got --publish_dir_mode link" has " --publish_dir_mode link " "${NFCMD} "
PBAM="${G2}/${P}/aligned/${P}_sorted.bam"
check_eq "the published BAM has two links (work/ and aligned/)" "$(stat -c %h "$PBAM" 2>/dev/null)" 2
cp "$PBAM" "${CASE_TMP}/copied.bam" 2>/dev/null
check_eq "negative control: a copy of it has one link" "$(stat -c %h "${CASE_TMP}/copied.bam" 2>/dev/null)" 1
rm -f "${CASE_TMP}/copied.bam"
# symlinks DIR: the symbolic links under DIR (not its nextflow/ launch directory), relative
symlinks() { (cd "$1" 2>/dev/null && find . -path ./nextflow -prune -o -type l -print | LC_ALL=C sort); }
SYMLINKS=$(LC_ALL=C comm -13 <(symlinks "$DIRECT") <(symlinks "${G2}/${P}"))
[ -z "$SYMLINKS" ] || printf 'symbolic links the direct (copy) run does not have:\n%s\n' "$SYMLINKS"
check_eq "published symbolic links the direct run does not have" "$(grep -c . <<<"$SYMLINKS" || true)" 0

# --- 1. the same published files as the direct run, and parity with the bash steps -------
published "$DIRECT" > "${CASE_TMP}/direct.txt"
published "${G2}/${P}" > "${CASE_TMP}/launcher.txt"
check_ge "files the direct run published" "$(wc -l < "${CASE_TMP}/direct.txt" | tr -d ' ')" 20
DIFF=$(diff "${CASE_TMP}/direct.txt" "${CASE_TMP}/launcher.txt")
[ -z "$DIFF" ] || printf 'published files that differ (< direct, > launcher):\n%s\n' "$DIFF"
check_eq "published files that differ from the direct run" "$(grep -c '^[<>]' <<<"$DIFF" || true)" 0
sed 1d "${CASE_TMP}/launcher.txt" > "${CASE_TMP}/planted.txt"
check "negative control: a planted missing file is reported" \
  test "$(diff "${CASE_TMP}/direct.txt" "${CASE_TMP}/planted.txt" | grep -c '^<' || true)" -eq 1
echo "+ scripts/ci/parity-diff.sh ${GENOME_DIR}/${P} ${G2}/${P} ${P}"
GITHUB_STEP_SUMMARY=/dev/null "${REPO}/scripts/ci/parity-diff.sh" "${GENOME_DIR}/${P}" "${G2}/${P}" "$P"
check_eq "parity of the launcher's outputs with the bash steps (exit)" "$?" 0

# --- 3. reports and step status -----------------------------------------------------------
check "the HTML report (step 24) is written" test -s "${G2}/${P}/${P}_report.html"
check "the text report is written" test -s "${G2}/${P}/${P}_report.txt"
check "summary.json validates" python3 "${REPO}/tests/schema/validate.py" "${REPO}/tests/schema/summary.schema.json" "${G2}/${P}/summary.json"
STATUS=$(cat "${G2}/${P}/logs/run_status.tsv" 2>/dev/null)
check "run_status.tsv records step 07 ok" has $'^step\t07\tok$' "$STATUS"
check "run_status.tsv records step 13 skipped" has $'^step\t13\tskipped ' "$STATUS"

# --- 2. the second run finds every task cached ---------------------------------------------
launch second
check_eq "the second run exits 0" "$RC" 0
SECOND_TRACE=$(latest_trace)
check "the second run wrote its own trace" test "${SECOND_TRACE:-}" != "${FIRST_TRACE:-}"
echo "second run: $(status_count "${SECOND_TRACE:-/dev/null}" CACHED) cached, $(status_count "${SECOND_TRACE:-/dev/null}" COMPLETED) completed, ${SECS}s"
check "every task of the second run came from the cache" all_cached "${SECOND_TRACE:-/dev/null}"
check "negative control: the first run's trace fails the all-cached check" not_all_cached "${FIRST_TRACE:-/dev/null}"
check "the second run takes under 10 minutes (${SECS}s)" test "$SECS" -lt 600
check_eq "the samplesheet is the same in the second run" "$(sed -n 2p "${G2}/${P}/nextflow/samplesheet.csv")" \
  "${P},${G2}/${P}/fastq/${P}_R1.fastq.gz,${G2}/${P}/fastq/${P}_R2.fastq.gz,,,,,male"

finish
