#!/usr/bin/env bash
# The Nextflow HTML report and three tables the bash report already reads:
# step 06's cDNA read check (<id>_clinvar_cdna_check.tsv), step 32's CYP2D6
# depth check (<id>_cyp2d6_depth_check.tsv) and step 31's phased
# compound-het pairs (<id>_compound_het_pairs.tsv).
#   1. HTML_REPORT, run for real in its image through a two-line pipeline, on
#      tables earlier cases wrote: HG002 (case 38: pypgx summary and depth
#      check), HG002cdna (case clinvar-clip-check-1: hits with one flagged)
#      and CHPH (case comphet-phasing-1-reads: prioritized VCF and the pairs).
#      It must place each table where bin/collect_summary.py reads it (a name
#      it has no place for stops the task), and the reports must carry the
#      depth-check line, the cDNA column with its flagged count, and the
#      pairs by phase. Negative control: the same samples without the three
#      tables; the same checks must find none of it.
#   2. main.nf passes the cDNA table on to HTML_REPORT: a stub run of every
#      tool (-profile test_all, no containers) stages <id>_clinvar_cdna_check.tsv
#      in each HTML_REPORT task that got ClinVar hits. Negative control: a copy
#      of the pipeline with that join removed stages none.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

C="${SAMPLE}cdna"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
PGX_SUM="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv"
PGX_DC="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv"
CV_HITS="${GENOME_DIR}/${C}/clinvar/${C}_clinvar_hits.vcf"
CV_CDNA="${GENOME_DIR}/${C}/clinvar/${C}_clinvar_cdna_check.tsv"
SL_VCF="${GENOME_DIR}/CHPH/slivar/CHPH_prioritized.vcf.gz"
SL_PAIRS="${GENOME_DIR}/CHPH/slivar/CHPH_compound_het_pairs.tsv"
for f in "$VCF" "$PGX_SUM" "$PGX_DC" "$CV_HITS" "$CV_CDNA" "$SL_VCF" "$SL_PAIRS"; do
  check "an earlier case left ${f#"$GENOME_DIR"/}" nonempty "${f#"$GENOME_DIR"/}"
done

# --- 1. HTML_REPORT on the three tables, and without them ----------------------------
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# report_run NAME WITH: HTML_REPORT on the three samples, with the tables (WITH=1) or without them
report_run() {
  local name=$1 with=$2 out="${CASE_TMP}/$1"
  rm -rf "$out"
  local dc="" cdna="" pairs=""
  [ "$with" = 1 ] && dc=", file('${PGX_DC}')" cdna=", file('${CV_CDNA}')" pairs=", file('${SL_PAIRS}')"
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    HTML_REPORT(channel.of(
        [[id: '${SAMPLE}', sex: 'male'], file('${VCF}'), [file('${PGX_SUM}')${dc}]],
        [[id: '${C}', sex: 'male'], file('${VCF}'), [file('${CV_HITS}')${cdna}]],
        [[id: 'CHPH', sex: 'male'], file('${VCF}'), [file('${SL_VCF}')${pairs}]]))
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${out}/results'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  echo "+ nextflow run (HTML_REPORT, ${name})"
  (
    cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
      -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work-${name}"
  ) > "${CASE_TMP}/${name}.log" 2>&1
  local rc=$?
  check_eq "HTML_REPORT (${name}) exits 0" "$rc" 0
  [ "$rc" = 0 ] || { grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/${name}.log" | head -40; }
}
# html_has OUT SAMPLE PATTERN: the sample's HTML report holds PATTERN (grep -E)
html_count() { grep -Ec "$3" "${CASE_TMP}/$1/results/$2/$2_report.html" 2>/dev/null || true; }

report_run with 1
for s in "$SAMPLE" "$C" CHPH; do
  check "HTML_REPORT wrote ${s}'s report" test -s "${CASE_TMP}/with/results/${s}/${s}_report.html"
done
check_ge "the CYP2D6 depth-check line from step 32's table" "$(html_count with "$SAMPLE" 'CYP2D6 depth check: (ok|unreliable);')" 1
check_ge "the cDNA read check column of the ClinVar hits" "$(html_count with "$C" '<th>cDNA read check</th>')" 1
check_ge "the flagged count of the cDNA read check, 1 for the planted hit" \
  "$(html_count with "$C" 'Flagged by the cDNA read check</span><span class="value"><span class="badge badge-red">1<')" 1
check_ge "the compound-het pairs by phase from step 31's table" "$(html_count with CHPH 'Pairs in trans / cis / unknown')" 1
check "the summary of CHPH counts the pairs table's rows" python3 - "${CASE_TMP}/with/results/CHPH/summary.json" "$SL_PAIRS" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))["sections"]["slivar"]["data"]
rows = sum(1 for _ in open(sys.argv[2])) - 1
print(d.get("compound_het_pairs"), rows)
sys.exit(d.get("compound_het_pairs") != rows)
PY

report_run without 0
check_eq "control: without the depth check there is no depth-check line" \
  "$(html_count without "$SAMPLE" 'CYP2D6 depth check: (ok|unreliable);')" 0
check_eq "control: without the cDNA table there is no cDNA column" "$(html_count without "$C" '<th>cDNA read check</th>')" 0
check_eq "control: without the pairs table there are no pairs by phase" \
  "$(html_count without CHPH 'Pairs in trans / cis / unknown')" 0

# --- 2. main.nf stages the cDNA table in HTML_REPORT (stub run) -----------------------
# stub_staged DIR: in the stub run under DIR, "<HTML_REPORT tasks with ClinVar hits> <of them with the cDNA table>"
stub_staged() {
  local hits=0 cdna=0 d id
  while IFS= read -r d; do
    d=$(dirname "$d")
    id=$(sed -n "s/^### name: [\"'].*HTML_REPORT (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
    [ -e "${d}/input/${id}_clinvar_hits.vcf" ] || continue
    hits=$((hits + 1))
    [ ! -e "${d}/input/${id}_clinvar_cdna_check.tsv" ] || cdna=$((cdna + 1))
  done < <(grep -l "^### name: [\"'].*HTML_REPORT (" -r "$1/work" --include=.command.run 2>/dev/null)
  echo "$hits $cdna"
}
# stub_run NAME: the test_all stub run in a copy of the pipeline under CASE_TMP/NAME
stub_run() {
  local d="${CASE_TMP}/$1"
  echo "+ nextflow run main.nf -profile test_all -stub ($1)"
  (cd "$d" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${d}/out" -work-dir "${d}/work") \
    > "${d}.log" 2>&1
  check_eq "stub run ($1) exits 0" "$?" 0
}
for n in stub stub-nojoin; do
  rm -rf "${CASE_TMP:?}/${n}"
  mkdir -p "${CASE_TMP}/${n}"
  cp -R "${REPO}/main.nf" "${REPO}/nextflow.config" "${REPO}/nextflow_schema.json" "${REPO}/versions.env" \
    "${REPO}/conf" "${REPO}/modules" "${REPO}/workflows" "${REPO}/bin" "${REPO}/assets" "${REPO}/scripts" "${CASE_TMP}/${n}/"
  mkdir -p "${CASE_TMP}/${n}/tests/fixtures"  # test_all's consumer_sites
  cp -R "${REPO}/tests/fixtures/consumer_export" "${CASE_TMP}/${n}/tests/fixtures/"
done
JOIN='PGX.out.clinvar_cdna.map'
check_eq "main.nf joins the cDNA table into the report inputs once" "$(grep -c "$JOIN" "${CASE_TMP}/stub/main.nf")" 1
sed -i "/${JOIN}/d" "${CASE_TMP}/stub-nojoin/main.nf"
check_eq "control copy: the join is gone" "$(grep -c "$JOIN" "${CASE_TMP}/stub-nojoin/main.nf")" 0
stub_run stub
read -r HITS CDNA <<< "$(stub_staged "${CASE_TMP}/stub")"
echo "stub run: ${HITS} HTML_REPORT tasks with ClinVar hits, ${CDNA} with the cDNA table"
check_ge "HTML_REPORT tasks with ClinVar hits in the stub run" "$HITS" 1
check_eq "each of them gets the cDNA table" "$CDNA" "$HITS"
stub_run stub-nojoin
read -r HITS0 CDNA0 <<< "$(stub_staged "${CASE_TMP}/stub-nojoin")"
echo "control stub run: ${HITS0} HTML_REPORT tasks with ClinVar hits, ${CDNA0} with the cDNA table"
check_ge "control: HTML_REPORT tasks with ClinVar hits" "$HITS0" 1
check_eq "control: without the join no task gets the cDNA table" "$CDNA0" 0

rm -rf "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-nojoin" "${CASE_TMP:?}/work-with" "${CASE_TMP:?}/work-without"
finish
