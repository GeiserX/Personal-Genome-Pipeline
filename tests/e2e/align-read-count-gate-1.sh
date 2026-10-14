#!/usr/bin/env bash
# The alignment read-count gate on real tools. minimap2 writes every read it
# is given, mapped or not, once, so the BAM's primary records (samtools view
# -c -F 0x900) must equal the FASTQ records that went in, R1 plus R2, which is
# also fastp's summary.after_filtering.total_reads.
#   1. Case 20's BAM (the raw pair): its table holds two equal counts.
#   2. The fixture through fastp and step 02: fastp's count is R1 plus R2
#      records (reads, not pairs), and the BAM holds exactly that many.
#   3. validate-setup.sh passes that BAM against fastp's report beside it,
#      and fails the same BAM subsampled with samtools view -s 0.6.
#   4. R2 cut to 12,000 of R1's 20,000 records: minimap2 stops at the end of
#      R2 with only a warning, so the BAM misses 8,000 reads. Step 02 stops
#      and keeps no BAM; so does the ALIGN_MINIMAP2 task, and with
#      ext.read_count_check = 'warn' the task keeps the BAM, the table says
#      mismatch and the report shows the missing reads.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
FQ="${G}/${SAMPLE}/fastq"
records() { gzip -cdf "$@" | awk 'END { printf "%.0f\n", NR / 4 }'; }
val() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "${G}/$2" 2>/dev/null; }
primary() { sam view -c -F 0x900 "$1" 2>/dev/null; }

# --- 1. case 20 ------------------------------------------------------------------------
T20="${SAMPLE}/qc/${SAMPLE}_read_count.tsv"
RAW=$(records "${FQ}/${SAMPLE}_R1.fastq.gz" "${FQ}/${SAMPLE}_R2.fastq.gz")
check "case 20's step 02 wrote ${T20}" nonempty "$T20"
check_eq "case 20: input_reads is R1 plus R2 records of the raw pair" "$(val input_reads "$T20")" "$RAW"
check_eq "case 20: primary_records is samtools view -c -F 0x900 of the BAM" \
  "$(val primary_records "$T20")" "$(primary "${SAMPLE}/aligned/${SAMPLE}_sorted.bam")"
check_eq "case 20: read_count_check" "$(val read_count_check "$T20")" ok

# --- 2. fastp, then step 02 on the trimmed reads ------------------------------------------
T=${SAMPLE}T
rm -rf "${G:?}/${T}"
mkdir -p "${G}/${T}/fastq"
for r in R1 R2; do ln -f "${FQ}/${SAMPLE}_${r}.fastq.gz" "${G}/${T}/fastq/${T}_${r}.fastq.gz"; done
run_step 01b-fastp-qc.sh "$T"
check_step_exit 01b-fastp-qc.sh
JSON="${G}/${T}/fastq_trimmed/${T}_fastp.json"
FASTP=$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["summary"]["after_filtering"]["total_reads"])' "$JSON" 2>/dev/null)
TRIM1=$(records "${G}/${T}/fastq_trimmed/${T}_R1.fastq.gz")
TRIM2=$(records "${G}/${T}/fastq_trimmed/${T}_R2.fastq.gz")
check_ge "fastp kept reads" "${FASTP:-0}" 100000
check_eq "fastp's total_reads is R1 plus R2 records of its output (reads, not pairs)" "$FASTP" "$((TRIM1 + TRIM2))"
run_step 02-alignment.sh "$T"
check_step_exit 02-alignment.sh
check "step 02 aligned the trimmed reads" has 'Using trimmed FASTQs from fastp' "$(cat "$STEP_LOG")"
TT="${T}/qc/${T}_read_count.tsv"
check_eq "input_reads equals fastp's total_reads" "$(val input_reads "$TT")" "$FASTP"
check_eq "the BAM's primary records equal fastp's total_reads" "$(primary "${T}/aligned/${T}_sorted.bam")" "$FASTP"
check_eq "read_count_check" "$(val read_count_check "$TT")" ok
python3 "${REPO}/bin/collect_summary.py" --sample "$T" --sample-dir "${G}/${T}" --out "${CASE_TMP}/T.json" >/dev/null
python3 "${REPO}/bin/render_report.py" --summary "${CASE_TMP}/T.json" -o "${CASE_TMP}/T.txt" >/dev/null
check "the text report shows both counts" \
  has "^  Alignment: ${FASTP} reads in, ${FASTP} primary records in the BAM\$" "$(cat "${CASE_TMP}/T.txt" 2>/dev/null)"

# --- 3. validate-setup.sh: the same BAM, then the BAM subsampled ------------------------
run_step validate-setup.sh "$T"
check "validate-setup passes the BAM against fastp's report" \
  has "\[OK\].*BAM holds the ${FASTP} reads fastp wrote as primary records" "$(cat "$STEP_LOG")"
S=${SAMPLE}S
rm -rf "${G:?}/${S}"
mkdir -p "${G}/${S}/aligned" "${G}/${S}/fastq_trimmed"
sam view -b -s 0.6 -o "${S}/aligned/${S}_sorted.bam" "${T}/aligned/${T}_sorted.bam"
sam index "${S}/aligned/${S}_sorted.bam"
cp "$JSON" "${G}/${S}/fastq_trimmed/${S}_fastp.json"
NS=$(primary "${S}/aligned/${S}_sorted.bam")
check_ge "the subsampled BAM still has reads" "${NS:-0}" 1
check "the subsampled BAM holds fewer primary records than fastp's count" test "${NS:-0}" -lt "$FASTP"
run_step validate-setup.sh "$S"
OUT=$(cat "$STEP_LOG")
check "validate-setup fails the subsampled BAM" \
  has "\[FAIL\].*BAM read count: ${NS} primary records, but fastp wrote ${FASTP} reads" "$OUT"
check "validate-setup says how many reads are missing" has "$((FASTP - NS)) reads \([0-9.]+%\) are missing" "$OUT"

# --- 4. R2 cut short: step 02 ---------------------------------------------------------------
X=${SAMPLE}X
rm -rf "${G:?}/${X}"
mkdir -p "${G}/${X}/fastq"
{ gzip -cd "${FQ}/${SAMPLE}_R1.fastq.gz" || true; } | head -n 80000 | gzip -c > "${G}/${X}/fastq/${X}_R1.fastq.gz"
{ gzip -cd "${FQ}/${SAMPLE}_R2.fastq.gz" || true; } | head -n 48000 | gzip -c > "${G}/${X}/fastq/${X}_R2.fastq.gz"
check_eq "R1 records of the cut pair" "$(records "${G}/${X}/fastq/${X}_R1.fastq.gz")" 20000
check_eq "R2 records of the cut pair" "$(records "${G}/${X}/fastq/${X}_R2.fastq.gz")" 12000
run_step 02-alignment.sh "$X"
OUT=$(cat "$STEP_LOG")
check "step 02 exits non-zero on the cut pair (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "minimap2 only warned about the cut pair" has 'query files have different number of records' "$OUT"
check "step 02 names the 8000 missing reads" \
  has 'ERROR: the BAM holds 24000 primary records, but 32000 reads went in \(.*\): 8000 reads \(25\.0%\) are missing\.' "$OUT"
check "step 02 kept no BAM" test ! -e "${G}/${X}/aligned/${X}_sorted.bam"
check "step 02 kept no read-count table" test ! -e "${G}/${X}/qc/${X}_read_count.tsv"

# --- 4b. The same pair through the ALIGN_MINIMAP2 task --------------------------------------
# Only the alignment matters here, so INDEXCOV and DEEPVARIANT are switched
# off (ext.when); the reference keeps the .sr.mmi case 20 built beside it.
nf_align() {
  local name=$1 conf=$2
  printf 'sample,fastq_1,fastq_2,sex\n%s,%s,%s,male\n' "$X" \
    "${G}/${X}/fastq/${X}_R1.fastq.gz" "${G}/${X}/fastq/${X}_R2.fastq.gz" > "${CASE_TMP}/${name}.csv"
  mkdir -p "${CASE_TMP}/${name}"
  (
    cd "${CASE_TMP}/${name}" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -c "$conf" \
      -work-dir "${E2E_WORK}/nf-work-${name}" \
      --input "${CASE_TMP}/${name}.csv" \
      --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
      --skip_trim \
      --tools mosdepth \
      --outdir "${G}/nf-${name}" \
      --max_cpus 4 --max_memory 14.GB
  ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
  cp "${CASE_TMP}/${name}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}-${name}.nextflow.log" 2>/dev/null || true
}
command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
cat > "${CASE_TMP}/off.config" <<'NFCONF'
process {
    withName: 'INDEXCOV'    { ext.when = false }
    withName: 'DEEPVARIANT' { ext.when = false }
}
NFCONF
nf_align readcount-fail "${CASE_TMP}/off.config"
OUT=$(cat "$STEP_LOG" "${CASE_TMP}/readcount-fail/.nextflow.log" 2>/dev/null)
check "the pipeline fails on the cut pair (exit ${NF_RC})" test "$NF_RC" -ne 0
check "ALIGN_MINIMAP2 is the failed task" \
  has '^ALIGN_MINIMAP2|:ALIGN_MINIMAP2' "$(awk -F'\t' '$2 == "failed" {print $1}' "${CASE_TMP}/readcount-fail/failed_tasks.tsv" 2>/dev/null)"
check "ALIGN_MINIMAP2 names the 8000 missing reads" \
  has 'ERROR: the BAM holds 24000 primary records, but 32000 reads went in .*: 8000 reads \(25\.0%\) are missing\.' "$OUT"
check "no BAM was published" test ! -e "${G}/nf-readcount-fail/${X}/aligned/${X}_sorted.bam"

{ cat "${CASE_TMP}/off.config"
  printf "process {\n    withName: 'ALIGN_MINIMAP2' { ext.read_count_check = 'warn' }\n}\n"; } > "${CASE_TMP}/warn.config"
nf_align readcount-warn "${CASE_TMP}/warn.config"
check_eq "with ext.read_count_check = 'warn' the pipeline exits 0" "$NF_RC" 0
R="nf-readcount-warn/${X}"
check "the BAM was published" nonempty "${R}/aligned/${X}_sorted.bam"
check_eq "the published table says mismatch" "$(val read_count_check "${R}/qc/${X}_read_count.tsv")" mismatch
check_eq "the published table: input_reads" "$(val input_reads "${R}/qc/${X}_read_count.tsv")" 32000
check_eq "the published table: primary_records" "$(val primary_records "${R}/qc/${X}_read_count.tsv")" 24000
check_eq "the published table: missing_reads" "$(val missing_reads "${R}/qc/${X}_read_count.tsv")" 8000
check_eq "the published BAM's primary records" "$(primary "${R}/aligned/${X}_sorted.bam")" 24000
check "the table is not published beside the BAM" test ! -e "${G}/${R}/aligned/${X}_read_count.tsv"
python3 "${REPO}/bin/collect_summary.py" --sample "$X" --sample-dir "${G}/${R}" --out "${CASE_TMP}/X.json" >/dev/null
python3 "${REPO}/bin/render_report.py" --summary "${CASE_TMP}/X.json" -o "${CASE_TMP}/X.txt" >/dev/null
python3 "${REPO}/bin/render_report.py" --summary "${CASE_TMP}/X.json" -o "${CASE_TMP}/X.html" >/dev/null
check "the text report flags the missing reads" \
  has '^  Alignment: 32000 reads in, 24000 primary records in the BAM: 8000 reads missing .*READS LOST IN ALIGNMENT' \
  "$(cat "${CASE_TMP}/X.txt" 2>/dev/null)"
check "the HTML report shows the missing reads in red" \
  has 'Reads kept by alignment</span><span class="value"><span class="badge badge-red">32000 reads in, 24000 primary records in the BAM: 8000 reads missing' \
  "$(cat "${CASE_TMP}/X.html" 2>/dev/null)"

rm -rf "${G:?}/${T}" "${G:?}/${S}" "${G:?}/${X}" "${G:?}/nf-readcount-fail" "${G:?}/nf-readcount-warn"
finish
