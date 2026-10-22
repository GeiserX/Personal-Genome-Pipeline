#!/usr/bin/env bash
# pgsc_calc gets the outer run's docker flags, and its own trace is kept.
#   - Nextflow, the HG002P row of case nextflow-from-fastq-1-bash (its VCF and
#     gVCF, male) with --tools prs, under an extra config whose
#     docker.runOptions adds --cpu-shares 300 and a label and whose
#     process.containerOptions adds another label. The images.config PRS
#     writes for pgsc_calc must carry both after -u and --network none (and
#     the sex mount), without the outer --user root:root, and Docker's create
#     events must show pgsc_calc's own containers (its pgscatalog images)
#     with both labels.
#   - Nextflow and step 25: pgsc_calc's trace, report and timeline are
#     published as prs/pgsc_calc/pipeline_info/execution_{trace.txt,
#     report.html,timeline.html}, and nothing else is there. The trace has
#     Nextflow's default columns only, a FORMAT_SCOREFILES row with its
#     peak_rss, and no sample label (nothing matching HG002).
# Each check also runs on a broken copy of its input and must fail there.
# Removes what it made.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

P=HG002P
G="$GENOME_DIR"
D="${G}/${P}"
for f in "vcf/${P}.vcf.gz" "vcf/${P}.vcf.gz.tbi" "vcf/${P}.g.vcf.gz" "vcf/${P}.g.vcf.gz.tbi"; do
  check "case nextflow-from-fastq-1-bash wrote ${f}" nonempty "${P}/${f}"
done
PGS="${E2E_WORK}/parity-pgs"
check "case nextflow-from-fastq-1-bash wrote the score files" test -s "${PGS}/PGS000018.txt.gz"

refused() { ! "$@"; }
COLS=$'task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar'
# trace_ok FILE: Nextflow's default trace columns and nothing else, a
# FORMAT_SCOREFILES row with a peak_rss, and no sample label.
trace_ok() {
  [ -s "$1" ] && [ "$(head -n 1 "$1")" = "$COLS" ] &&
    awk -F'\t' 'NR > 1 && $4 ~ /FORMAT_SCOREFILES/ && $11 != "" && $11 != "-" {f = 1} END {exit !f}' "$1" &&
    ! grep -q 'HG002' "$1"
}
# info_ok DIR: the three files under fixed names, nothing else, and the
# report and the timeline name FORMAT_SCOREFILES.
info_ok() {
  [ "$(cd "$1" 2>/dev/null && find . -type f | LC_ALL=C sort | tr '\n' ' ')" = \
    "./execution_report.html ./execution_timeline.html ./execution_trace.txt " ] &&
    trace_ok "$1/execution_trace.txt" &&
    grep -q FORMAT_SCOREFILES "$1/execution_report.html" &&
    grep -q FORMAT_SCOREFILES "$1/execution_timeline.html"
}
# check_info LABEL DIR: info_ok on DIR, and its controls on broken copies.
check_info() {
  local c="${CASE_TMP}/broken"
  echo "${1}: $(cd "$2" 2>/dev/null && find . -type f | LC_ALL=C sort | tr '\n' ' ')"
  head -n 3 "$2/execution_trace.txt" 2>/dev/null
  check "${1}: pgsc_calc's trace, report and timeline are published in ${2#"${G}/"}" info_ok "$2"
  rm -rf "$c" && cp -r "$2" "$c" 2>/dev/null
  sed -i '2 s/\t/ HG002\t/4' "${c}/execution_trace.txt" 2>/dev/null
  check "${1}: control: a trace naming the sample is refused" refused trace_ok "${c}/execution_trace.txt"
  rm -rf "$c" && cp -r "$2" "$c" 2>/dev/null
  sed -i '1 s/$/\tworkdir/' "${c}/execution_trace.txt" 2>/dev/null
  check "${1}: control: a trace with another column is refused" refused trace_ok "${c}/execution_trace.txt"
  rm -rf "$c" && cp -r "$2" "$c" 2>/dev/null
  sed -i '/FORMAT_SCOREFILES/d' "${c}/execution_trace.txt" 2>/dev/null
  check "${1}: control: a trace without FORMAT_SCOREFILES is refused" refused trace_ok "${c}/execution_trace.txt"
  rm -rf "$c" && cp -r "$2" "$c" 2>/dev/null
  touch "${c}/pipeline_dag_2026-01-01_00-00-00.html"
  check "${1}: control: another file in pipeline_info is refused" refused info_ok "$c"
  rm -rf "$c"
}

# --- Nextflow -------------------------------------------------------------------------
NF_OUT="${G}/nf-prs-trace"
WORK="${E2E_WORK}/nf-work-prs-trace"
SHEET="${CASE_TMP}/samplesheet.csv"
rm -rf "$NF_OUT" "$WORK"
printf 'sample,vcf,vcf_index,gvcf,gvcf_index,sex\n%s,%s,%s,%s,%s,male\n' "$P" \
  "${D}/vcf/${P}.vcf.gz" "${D}/vcf/${P}.vcf.gz.tbi" "${D}/vcf/${P}.g.vcf.gz" "${D}/vcf/${P}.g.vcf.gz.tbi" > "$SHEET"
# The outer flags: a low-priority share and a label each way. The docker
# profile's --user root:root and --network none are repeated, as a user's
# config would, and must not reach pgsc_calc.
cat > "${CASE_TMP}/outer.config" <<'NFCONF'
docker.runOptions = '--user root:root --cpu-shares 300 --label pgp.runopts=e2e'
process.containerOptions = '--network none --label pgp.containeropts=e2e'
NFCONF

# Every container created during the run: image, then the two labels.
EVENTS="${CASE_TMP}/events.tsv"
docker events --filter type=container --filter event=create \
  --format $'{{.Actor.Attributes.image}}\t{{index .Actor.Attributes "pgp.runopts"}}\t{{index .Actor.Attributes "pgp.containeropts"}}' \
  > "$EVENTS" 2> "${CASE_TMP}/events.err" &
EV_PID=$!
trap 'kill "$EV_PID" 2>/dev/null' EXIT

echo "+ nextflow run main.nf -profile docker -c outer.config --tools prs (VCF + gVCF row)"
(
  cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -c "${CASE_TMP}/outer.config" \
    -work-dir "$WORK" \
    --input "$SHEET" \
    --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --sex_check warn \
    --tools prs \
    --pgs_scoring "$PGS" \
    --outdir "$NF_OUT" \
    --max_cpus 4 --max_memory 14.GB
) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
sleep 2
kill "$EV_PID" 2>/dev/null
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
TRACE=$(find "${NF_OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
check "the trace lists the tasks" test -s "${TRACE:-/dev/null}"

# images.config of the PRS task, found by its hash in the trace.
HASH=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
  { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == "PRS") {print $c["hash"]; exit} }' "${TRACE:-/dev/null}" 2>/dev/null)
CONFIG=$(find "$WORK/${HASH%%/*}" -maxdepth 2 -path "${WORK}/${HASH:-none}*/images.config" 2>/dev/null | head -n 1)
echo "PRS task ${HASH:-none}: $(grep '^docker' "${CONFIG:-/dev/null}" 2>/dev/null)"
# runopts_ok FILE: -u and --network none, the sex mount, then the outer flags
# with their user and network flags dropped.
runopts_ok() {
  grep -Eq "^docker\.runOptions = '-u [0-9]+:[0-9]+ --network none -v \"[^\"]+/sex:[^\"]+/sex:ro\" --cpu-shares 300 --label pgp\.runopts=e2e --label pgp\.containeropts=e2e'\$" "$1"
}
check "pgsc_calc's docker.runOptions carries the outer flags after -u and --network none" runopts_ok "${CONFIG:-/dev/null}"
sed "s/^docker\.runOptions = .*/docker.runOptions = '-u 1001:118 --network none -v \"\/w\/sex:\/w\/sex:ro\"'/" "${CONFIG:-/dev/null}" > "${CASE_TMP}/main.config"
check "control: the line before this package (no outer flags) is refused" refused runopts_ok "${CASE_TMP}/main.config"
sed "s/'\$/ --user root:root'/" "${CONFIG:-/dev/null}" > "${CASE_TMP}/root.config"
check "control: a line that keeps the outer --user root:root is refused" refused runopts_ok "${CASE_TMP}/root.config"

# Docker's view: pgsc_calc's containers got both labels.
echo "containers created: $(wc -l < "$EVENTS" | tr -d ' '); events errors: $(head -c 300 "${CASE_TMP}/events.err")"
awk -F'\t' '{print $1 "\t" $2 "\t" $3}' "$EVENTS" | LC_ALL=C sort | uniq -c
inner()    { awk -F'\t' '$1 ~ /pgscatalog/ {n++} END {print n + 0}' "$1"; }
labelled() { awk -F'\t' '$1 ~ /pgscatalog/ && $2 == "e2e" && $3 == "e2e" {n++} END {print n + 0}' "$1"; }
INNER=$(inner "$EVENTS")
check_ge "pgsc_calc's own containers seen (pgscatalog images)" "$INNER" 3
check_eq "pgsc_calc's containers with both outer labels" "$(labelled "$EVENTS")" "$INNER"
check_ge "positive control: the outer tasks' containers show the labels" \
  "$(awk -F'\t' '$1 !~ /pgscatalog/ && $2 == "e2e" && $3 == "e2e" {n++} END {print n + 0}' "$EVENTS")" 1
awk -F'\t' -v OFS='\t' '{$2 = ""; $3 = ""; print}' "$EVENTS" > "${CASE_TMP}/events_unlabelled.tsv"
check "control: containers without the labels are counted as unlabelled" \
  test "$(labelled "${CASE_TMP}/events_unlabelled.tsv")" -eq 0

check_info "nextflow" "${NF_OUT}/${P}/prs/pgsc_calc/pipeline_info"
check_eq "PRS scores in the summary (one per id of assets/pgs_scores.tsv)" "$(awk 'NR > 1' "${NF_OUT}/${P}/prs/${P}_prs_summary.tsv" 2>/dev/null | wc -l | tr -d ' ')" \
  "$(awk -F'\t' '$1 ~ /^PGS[0-9]+$/' "${REPO}/assets/pgs_scores.tsv" | wc -l | tr -d ' ')"

# --- step 25 ---------------------------------------------------------------------------
# A copy of HG002P's VCF and gVCF under another label, so HG002P's own
# outputs, which later cases read, stay as they are.
T=HG002T
rm -rf "${G:?}/${T}" "${G}/prs_scores"
mkdir -p "${G}/${T}/vcf" "${G}/prs_scores"
for s in vcf.gz vcf.gz.tbi g.vcf.gz g.vcf.gz.tbi; do cp "${D}/vcf/${P}.${s}" "${G}/${T}/vcf/${T}.${s}"; done
cp "${PGS}"/*.txt.gz "${G}/prs_scores/"
ANCESTRY_PANEL=none run_step 25-prs.sh "$T"
check_step_exit 25-prs.sh
check_info "step 25" "${G}/${T}/prs/pgsc_calc/pipeline_info"
check "step 25: the run reports are no longer under results/" test ! -e "${G}/${T}/prs/pgsc_calc/results/pipeline_info"

printf '#### pgsc_calc: outer docker flags and its own trace\n\n%s of %s pgsc_calc containers had both outer labels; trace published by both entry points\n\n' \
  "$(labelled "$EVENTS")" "$INNER" >> "$E2E_NOTES"
rm -rf "${G:?}/${T}" "${G}/prs_scores" "$NF_OUT" "$WORK"
finish
