#!/usr/bin/env bash
# A VCF-first Nextflow run (no BAM) with the tools a vendor-VCF user picks:
#   - the HTML report shows ROH, the mito haplogroup and CPIC, with the same
#     ROH numbers and haplogroup bin/collect_summary.py (the bash report) reads
#     from the same files, and it starts after the steps it shows;
#   - multiqc is selected but has nothing to read: the log says so in one line;
#   - the completion message prints, with no onComplete handler error;
#   - CLINVAR_SCREEN writes the hits TSV, lists only carried alleles, and adds
#     no command line of its own (with ClinVar's absolute path) to the hits VCF.
# Its PharmCAT, ROH and haplogroup outputs are the reference the later
# vendor-vcf-intake cases compare against.
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

sheet "${CASE_TMP}/chr.csv" "$CHR_VCF"
nf_run chr "${CASE_TMP}/chr.csv" clinvar,pharmcat,cpic,roh,mito_haplogroup,html_report,multiqc \
  --clinvar "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
  --clinvar_index "${G}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi"
check_eq "nextflow run exits 0" "$NF_RC" 0
LOG=$(cat "$NF_LOG")
R="${NF_OUT}/${SAMPLE}"

# --- Completion handler and MultiQC --------------------------------------------
check "the completion message prints" has 'Pipeline completed successfully!' "$LOG"
check "no onComplete handler error" lacks 'Failed to invoke .workflow.onComplete. event handler' "$LOG"
check_eq "one log line says MultiQC was skipped" "$(grep -c 'multiqc skipped: ' <<< "$LOG" || true)" 1
check "the skip line names the reason (mosdepth)" has 'multiqc skipped: .*mosdepth' "$LOG"
check "no multiqc folder" test ! -e "${NF_OUT}/multiqc"

# --- HTML report: ROH, mito haplogroup and CPIC ----------------------------------
HTML="${R}/${SAMPLE}_report.html"
check "the report exists" test -s "$HTML"
for h in 'Runs of Homozygosity' 'Mitochondrial Haplogroup' 'CPIC Drug Recommendations'; do
  check_eq "report card '${h}' is filled" "$(html_stat "$HTML" "$h" Status)" Complete
done
python3 "${REPO}/bin/collect_summary.py" --sample "$SAMPLE" --sample-dir "$R" --out "${CASE_TMP}/summary.json" \
  > "${CASE_TMP}/summary.log" 2>&1 || cat "${CASE_TMP}/summary.log"
want() { python3 -c "import json, sys; d = json.load(open('${CASE_TMP}/summary.json'))['sections']$1; print($2)" 2>/dev/null; }
ROH_TOTAL=$(want "['roh']['data']['total_mb']" "'%.1f MB' % d")
ROH_LARGEST=$(want "['roh']['data']['largest_mb']" "'%.1f MB' % d")
HG=$(want "['haplogroup']['data']['haplogroup']" "d")
echo "collect_summary.py: ROH total ${ROH_TOTAL:-?}, largest ${ROH_LARGEST:-?}, haplogroup ${HG:-?}"
check "collect_summary.py read the ROH file" test -n "$ROH_TOTAL"
check_eq "ROH total matches the bash report's" "$(html_stat "$HTML" 'Runs of Homozygosity' 'ROH total')" "$ROH_TOTAL"
check_eq "ROH largest segment matches the bash report's" \
  "$(html_stat "$HTML" 'Runs of Homozygosity' 'ROH largest segment')" "$ROH_LARGEST"
check "collect_summary.py read the haplogroup" test -n "$HG"
check_eq "haplogroup matches the bash report's" "$(html_stat "$HTML" 'Mitochondrial Haplogroup' Haplogroup)" "$HG"

# The report waits for the steps it shows.
REPORT_SUBMIT=$(trace_col HTML_REPORT submit)
for p in ROH MITO_HAPLOGROUP CPIC_LOOKUP CLINVAR_SCREEN; do
  done_at=$(trace_col ":${p} " complete)
  if [ -n "$REPORT_SUBMIT" ] && [ -n "$done_at" ] && [[ ! "$REPORT_SUBMIT" < "$done_at" ]]; then
    pass "HTML_REPORT submitted (${REPORT_SUBMIT}) after ${p} completed (${done_at})"
  else
    fail "HTML_REPORT submitted at '${REPORT_SUBMIT}', before ${p} completed at '${done_at}'"
  fi
done

# --- ClinVar screen ----------------------------------------------------------------
HITS="${R}/clinvar/${SAMPLE}_clinvar_hits.vcf"
TSV="${R}/clinvar/${SAMPLE}_clinvar_hits.tsv"
check_ge "CLINVAR_SCREEN reports the planted hit" "$(sample_side_hits "${R}/clinvar" "")" 1
check_eq "hits TSV has one row per hit in the VCF" \
  "$(awk 'NR > 1' "$TSV" 2>/dev/null | grep -c . || true)" "$(grep -vc '^#' "$HITS" 2>/dev/null || true)"
check_ge "the planted hit is a TSV row with its gene" \
  "$(awk -F'\t' -v p="$(planted pos)" -v g="$(planted gene)" '$2 == p && index($7, g)' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 1
check_eq "every hit's genotype carries an ALT allele" \
  "$(grep -v '^#' "$HITS" 2>/dev/null | awk -F'\t' '{split($10, f, ":"); if (f[1] !~ /[1-9]/) n++} END {print n + 0}')" 0
check_eq "no isec/annotate command line in the hits VCF" \
  "$(grep -cE '^##bcftools_(isec|annotate)Command' "$HITS" || true)" 0
# Any bcftools command line with an absolute path came with the input VCF.
gzip -dc "$CHR_VCF" | grep '^##' > "${CASE_TMP}/input_header.txt"
NEW_PATHS=$(grep -E '^##bcftools_.*Command=.* /' "$HITS" | grep -vxF -f "${CASE_TMP}/input_header.txt" | grep -c . || true)
check_eq "command lines with an absolute path the screen added" "$NEW_PATHS" 0

# --- Kept for the later cases -----------------------------------------------------------
mkdir -p "${INTAKE}/chr-run"
cp "${R}/roh/${SAMPLE}_roh.txt" "${R}/mito/${SAMPLE}_haplogroup.txt" "${INTAKE}/chr-run/" 2>/dev/null
cp "${R}/pharmcat/${SAMPLE}.match.json" "${R}/pharmcat/${SAMPLE}.phenotype.json" "${INTAKE}/chr-run/" 2>/dev/null
check "ROH, haplogroup and PharmCAT outputs kept for the later cases" \
  test -s "${INTAKE}/chr-run/${SAMPLE}_roh.txt" -a -s "${INTAKE}/chr-run/${SAMPLE}_haplogroup.txt" \
    -a -s "${INTAKE}/chr-run/${SAMPLE}.match.json"

finish
