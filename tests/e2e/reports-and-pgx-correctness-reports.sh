#!/usr/bin/env bash
# Both reports from one summary (bin/collect_summary.py, bin/render_report.py)
# and the run manifest (bin/write_manifest.sh), on the outputs of the cases
# before this one:
#   1. HG002: summary.json validates against tests/schema/summary.schema.json;
#      the text and HTML reports show the same ClinVar count, PRS rows and mean
#      depth as the summary; the HTML ClinVar table has a Stars column;
#   2. run_manifest.tsv holds a digest for every image on this machine that
#      versions.env names, and the ClinVar file date;
#   3. a report over the Nextflow output folder of case 60 shows ROH (the old
#      reports looked for it under vcf/ only);
#   4. after a run where step 31 was skipped (logs/run_status.tsv), the old
#      slivar result is marked stale in both reports, and step 06's is not;
#   5. scripts/run-all.sh with the real Nextflow and stub tasks (-stub: each
#      task only touches its outputs, inside its image), in a GENOME_DIR of its
#      own. A clean run of ROH and the ClinVar screen, then the same command
#      with Manta, duphold and AnnotSV added and duphold made to fail. The
#      engine's process_status.tsv says ROH, ClinVar and indexcov came from the
#      cache, Manta completed, duphold failed and AnnotSV got no task, and
#      run_status.tsv records them ok, ok, ok, ok, failed and not finished. An
#      old ROH file of the cached step is then shown as current in both
#      reports. Two negative controls in the same run: with step 11 rewritten
#      to 'not finished' the report marks that file stale, and a run whose
#      process_status.tsv is removed records no cached step ok.
# The PRS summary is planted (step 25 needs downloads the fixture lacks); it
# is removed at the end.
. "$(dirname "$0")/lib.sh"

D="${GENOME_DIR}/${SAMPLE}"
PRS="${D}/prs/${SAMPLE}_prs_summary.tsv"
STATUS="${D}/logs/run_status.tsv"
cleanup() { rm -f "$PRS" "$STATUS"; rmdir "${D}/prs" 2>/dev/null || true; }
trap cleanup EXIT
mkdir -p "${D}/prs"
printf 'Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\nSynthetic trait A\tPGS999901\t0.123\t10\t20\nSynthetic trait B\tPGS999902\t-0.456\t5\t9\n' > "$PRS"
rm -f "${D}/run_manifest.tsv"

run_step 24-html-report.sh "$SAMPLE"
check_step_exit 24-html-report.sh
run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh
HTML="${D}/${SAMPLE}_report.html"
TXT="${D}/${SAMPLE}_report.txt"
JSON="${D}/summary.json"
cp "$JSON" "${E2E_WORK}/logs/${SAMPLE}.summary.json" 2>/dev/null || true
cp "$HTML" "${E2E_WORK}/logs/${SAMPLE}_report.html" 2>/dev/null || true
cp "$TXT" "${E2E_WORK}/logs/${SAMPLE}_report.txt" 2>/dev/null || true

# 1. one summary, two reports
check "summary.json validates against the schema" \
  python3 "${REPO}/tests/schema/validate.py" "${REPO}/tests/schema/summary.schema.json" "$JSON"
jq_py() { python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print($1)" "$JSON" 2>/dev/null; }
N=$(jq_py "d['sections']['clinvar']['data'].get('count', 'none')")
echo "ClinVar hits in summary.json: ${N}"
check_ge "summary.json counts the planted ClinVar hit" "${N:-0}" 1
check "the text report shows the same ClinVar count" grep -q "^  Pathogenic/Likely Pathogenic hits: ${N}\$" "$TXT"
check "the HTML report shows the same ClinVar count" grep -qE "class=\"badge badge-(green|yellow)\">${N}</span>" "$HTML"
check "the HTML ClinVar table has a Stars column" grep -q '<th>Review status</th><th>Stars</th>' "$HTML"
for id in PGS999901 PGS999902; do
  check "PRS row ${id} in the text report" grep -q "${id}" "$TXT"
  check "PRS row ${id} in the HTML report" grep -q "<td>${id}</td>" "$HTML"
done
MOS=$(awk -F'\t' '$1 == "total" {printf "%.1f", $4}' "${D}/mosdepth/${SAMPLE}.mosdepth.summary.txt" 2>/dev/null)
echo "mosdepth mean depth: ${MOS:-none}"
check "the text report shows mosdepth's mean depth (${MOS:-none}x)" grep -q "Mean depth: ${MOS:-none}x" "$TXT"
check "the HTML report shows the same mean depth" grep -q ">${MOS:-none}x<" "$HTML"
check "the text report has the 'not assessed' list" grep -q '^## Not Assessed by This Pipeline' "$TXT"

# 2. the manifest
MAN="${D}/run_manifest.tsv"
check "24-html-report.sh wrote run_manifest.tsv" grep -q $'^run\twritten_by\t24-html-report.sh' "$MAN"
missing=0 present=0
while IFS= read -r line; do
  [[ "$line" =~ ^[A-Z0-9_]+_IMAGE= ]] || continue
  ref=$(cut -d= -f2- <<< "$line" | sed 's/[[:space:]]*#.*//; s/^"//; s/"$//')
  docker image inspect "$ref" >/dev/null 2>&1 || continue
  present=$((present + 1))
  if ! awk -F'\t' -v r="$ref" '$1 == "image" && $2 == r && $3 ~ /sha256:[0-9a-f]{64}/ {f = 1} END {exit !f}' "$MAN"; then
    echo "  no digest for ${ref}"; missing=$((missing + 1))
  fi
done < "${REPO}/versions.env"
check_ge "images of versions.env on this machine" "$present" 2
check_eq "of those, images without a digest in the manifest" "$missing" 0
DATE=$(awk -F'\t' '$1 == "data" && $2 == "clinvar_file_date" {print $3}' "$MAN")
echo "ClinVar file date: ${DATE:-none}"
check "the manifest has the ClinVar file date" has '^[0-9]{8}|^[0-9]{4}-[0-9]{2}-[0-9]{2}' "${DATE:-none}"
check "the HTML footer prints the manifest" grep -q 'Run manifest (run_manifest.tsv)' "$HTML"
check "the text report shows the ClinVar file date" grep -qF "ClinVar file date: ${DATE:-none}" "$TXT"

# 3. a Nextflow output folder
NF="${GENOME_DIR}/nf-results"
GENOME_DIR="$NF" run_step generate-report.sh "$SAMPLE"
check_step_exit "generate-report.sh (GENOME_DIR=nf-results)"
NFTXT="${NF}/${SAMPLE}/${SAMPLE}_report.txt"
grep -A4 '^## Runs of Homozygosity' "$NFTXT" 2>/dev/null
check "the report over the Nextflow folder shows ROH" grep -q '^## Runs of Homozygosity' "$NFTXT"
check "with its segments" grep -qE '^  Segments: [0-9]+ ' "$NFTXT"
check "and reads the CPIC table from cpic/" grep -q '^## CPIC Drug-Gene Recommendations' "$NFTXT"
check "the Nextflow summary validates" \
  python3 "${REPO}/tests/schema/validate.py" "${REPO}/tests/schema/summary.schema.json" "${NF}/${SAMPLE}/summary.json"

# 4. a step skipped in the latest run
SLIVAR="${D}/slivar/${SAMPLE}_slivar_summary.tsv"
check "case 42 left a slivar summary" test -s "$SLIVAR"
touch -d '2 days ago' "$SLIVAR"
mkdir -p "${D}/logs"
{
  printf 'meta\tstarted_epoch\t%s\nmeta\tdeclared_sex\tmale\n' "$(date -d '1 day ago' +%s)"
  for s in 06 07 11 12 16 16b 20 21 27 32; do printf 'step\t%s\tok\n' "$s"; done
  printf 'step\t31\tskipped (needs VEP, step 13 failed)\n'
} > "$STATUS"
check "the run status file is in place" grep -q $'^step\t31\tskipped' "$STATUS"
touch -d '3 days ago' "${D}/clinvar/${SAMPLE}_clinvar_hits.vcf"
run_step generate-report.sh "$SAMPLE"
check_step_exit generate-report.sh
run_step 24-html-report.sh "$SAMPLE"
check_step_exit 24-html-report.sh
grep -A3 '^## Variant Prioritization (slivar)' "$TXT"
check "text: step 31's old result is marked stale" \
  grep -q 'STALE: from an earlier run.*step 31 was skipped (needs VEP, step 13 failed)' "$TXT"
check "HTML: step 31's old result is marked stale" grep -q '<div class="stale">STALE: from an earlier run.*step 31 was skipped' "$HTML"
check "text: step 06 (ok in that run) is not stale" lacks 'STALE' "$(grep -A3 '^## ClinVar Pathogenic Screen' "$TXT")"

# 5. a pipeline one tool stops: each step is recorded from its own processes
G3="${E2E_WORK}/genome-status"
S3="${G3}/${SAMPLE}"
ST3="${S3}/logs/run_status.tsv"
PT3="${S3}/nextflow/process_status.tsv"
TXT3="${S3}/${SAMPLE}_report.txt"
HTML3="${S3}/${SAMPLE}_report.html"
ROH3="${S3}/roh/${SAMPLE}_roh.txt"
PLANT="${CASE_TMP}/plant.config"
check "nextflow is on PATH" command -v nextflow
rm -rf "$G3"
mkdir -p "${G3}/reference" "${G3}/clinvar" "${S3}/aligned" "${S3}/vcf" "${G3}/annotsv_annotations/Annotations_Human/Genes/GRCh38"
for e in fasta fasta.fai dict; do ln -f "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.${e}" "${G3}/reference/"; done
ln -f "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" "${G3}/clinvar/"
ln -f "${D}/aligned/${SAMPLE}_sorted.bam" "${D}/aligned/${SAMPLE}_sorted.bam.bai" "${S3}/aligned/"
ln -f "${D}/vcf/${SAMPLE}.vcf.gz" "${D}/vcf/${SAMPLE}.vcf.gz.tbi" "${S3}/vcf/"
# status_run NAME TOOLS: run-all.sh on G3 with stub tasks and the planted config; exit code in RC
status_run() {
  echo "+ scripts/run-all.sh ${SAMPLE} male -stub -ansi-log false -c ${PLANT} (${1}, TOOLS=${2})"
  env GENOME_DIR="$G3" SKIP_VALIDATION=true TOOLS="$2" "${REPO}/scripts/run-all.sh" "$SAMPLE" male -stub -ansi-log false -c "$PLANT" 2>&1 \
    | tee "${E2E_WORK}/logs/${CASE_NAME}.${1}.log"
  RC=${PIPESTATUS[0]}
  echo "+ exit ${RC}"
  cp "$ST3" "${E2E_WORK}/logs/${CASE_NAME}.${1}.run_status.tsv" 2>/dev/null
  cp "$PT3" "${E2E_WORK}/logs/${CASE_NAME}.${1}.process_status.tsv" 2>/dev/null
  echo "process_status.tsv:"; cat "$PT3" 2>/dev/null
  echo "run_status.tsv:"; grep '^step' "$ST3" 2>/dev/null | grep -v $'\tskipped'
}
# pt PROCESS COLUMN: a count of the process's row in process_status.tsv
#   (2 succeeded, 3 cached, 4 failed, 5 ignored, 6 unfinished; 0 for all five summed)
pt() {
  awk -F'\t' -v p="$1" -v c="$2" '$1 ~ ("(^|:)" p "$") {print (c == 0 ? $2 + $3 + $4 + $5 + $6 : $c); f = 1} END {if (!f) print "no row"}' "$PT3" 2>/dev/null
}
st() { awk -F'\t' -v s="$1" '$1 == "step" && $2 == s {print $3}' "$ST3" 2>/dev/null; }

printf '// nothing planted\n' > "$PLANT"
status_run clean roh,clinvar
check_eq "clean run: run-all.sh exits 0" "$RC" 0
check_eq "clean run: ROH ran one task" "$(pt ROH 2)" 1
check_eq "clean run: ROH has no cached task" "$(pt ROH 3)" 0
check_eq "clean run: step 11 is ok" "$(st 11)" ok
echo "ROH file after the clean run: $(stat -c '%y' "$ROH3" 2>/dev/null)"

printf "process { withName: 'DUPHOLD' { beforeScript = 'sleep 10; echo planted failure >&2; exit 1' } }\n" > "$PLANT"
status_run failed roh,clinvar,manta,duphold,annotsv
check "failed run: run-all.sh passes the pipeline's failure on (exit ${RC})" test "$RC" -ne 0
echo "ROH file after the failed run: $(stat -c '%y' "$ROH3" 2>/dev/null)"
check_eq "failed run: ROH came from the cache" "$(pt ROH 3)" 1
check_eq "failed run: ROH ran no task again" "$(pt ROH 2)" 0
check_eq "failed run: the ClinVar screen came from the cache" "$(pt CLINVAR_SCREEN 3)" 1
check_eq "failed run: indexcov came from the cache" "$(pt INDEXCOV 3)" 1
check_eq "failed run: Manta completed" "$(pt MANTA 2)" 1
check_eq "failed run: duphold failed" "$(pt DUPHOLD 4)" 1
check_eq "failed run: AnnotSV has a row and no task" "$(pt ANNOTSV 0)" 0
check "failed run: failed_tasks.tsv names DUPHOLD" grep -qE $'(^|:)DUPHOLD\tfailed$' "${S3}/nextflow/failed_tasks.tsv"
for s in 16 11 06 04; do check_eq "failed run: step ${s} (completed or cached) is ok" "$(st "$s")" ok; done
check_eq "failed run: step 15 (duphold) is failed" "$(st 15)" failed
check_eq "failed run: step 05 (AnnotSV, never started) is not finished" "$(st 05)" "not finished"

# The cached step's file, made older than the run: current in both reports
touch -d '2 days ago' "$ROH3"
GENOME_DIR="$G3" run_step generate-report.sh "$SAMPLE"
check_step_exit "generate-report.sh (GENOME_DIR=genome-status)"
GENOME_DIR="$G3" run_step 24-html-report.sh "$SAMPLE"
check_step_exit "24-html-report.sh (GENOME_DIR=genome-status)"
check "text: the ROH section is there" grep -q '^## Runs of Homozygosity' "$TXT3"
check "text: ROH (cached in the failed run) is not stale" lacks 'step 11 was' "$(cat "$TXT3" 2>/dev/null)"
check "HTML: ROH (cached in the failed run) is not stale" lacks 'step 11 was' "$(cat "$HTML3" 2>/dev/null)"
cp "$ST3" "${CASE_TMP}/run_status.kept"
sed -i $'s/^step\t11\tok$/step\t11\tnot finished/' "$ST3"
GENOME_DIR="$G3" run_step generate-report.sh "$SAMPLE"
check "negative control: with step 11 'not finished' the same file is marked stale" \
  grep -q 'STALE: from an earlier run.*step 11 was not finished' "$TXT3"
cp "${CASE_TMP}/run_status.kept" "$ST3"

# Negative control: the same failed run without the engine's process table
mkdir -p "${CASE_TMP}/nfbin"
# shellcheck disable=SC2016  # the dollar signs go into the wrapper as written
printf '#!/usr/bin/env bash\n# the real nextflow, then its process table removed\n"%s" "$@"\nrc=$?\nrm -f process_status.tsv\nexit "$rc"\n' \
  "$(command -v nextflow)" > "${CASE_TMP}/nfbin/nextflow"
chmod +x "${CASE_TMP}/nfbin/nextflow"
PATH_KEPT=$PATH
PATH="${CASE_TMP}/nfbin:${PATH}"
status_run notable roh,clinvar,manta,duphold,annotsv
PATH=$PATH_KEPT
check "negative control, no process table: the run fails as before (exit ${RC})" test "$RC" -ne 0
check "negative control, no process table: the file is gone" test ! -e "$PT3"
check_eq "negative control, no process table: step 11 is not finished" "$(st 11)" "not finished"
check_eq "negative control, no process table: step 15 is still failed" "$(st 15)" failed

finish
