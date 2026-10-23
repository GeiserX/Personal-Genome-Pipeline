#!/usr/bin/env bash
# The Nextflow HTML report and step 35's HBA1/HBA2 table
# (<id>_hba_copy_number.tsv, bin/hba_copy_number.py), which the bash report
# already reads beside the SMN table.
#   1. HTML_REPORT, run for real in its image through a two-line pipeline, on
#      HG002's VCF with an SMN table and an HBA table that
#      bin/hba_copy_number.py writes from made-up depths (four genes). It
#      must place the HBA table in paralogs/ (a name it has no place for
#      stops the task), the report must carry the HBA line with the count,
#      and the not-assessed list must leave HBA1/HBA2 out. Negative control:
#      the same sample without the HBA table has no HBA line and lists
#      HBA1/HBA2 as not assessed.
#   2. main.nf passes the table on: a stub run of every tool (-profile
#      test_all, no containers) stages <id>_hba_copy_number.tsv in each
#      HTML_REPORT task that got the SMN table. Negative control: a copy of
#      the pipeline with that join removed stages none.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "an earlier case left ${VCF#"$GENOME_DIR"/}" nonempty "${VCF#"$GENOME_DIR"/}"

# --- 1. HTML_REPORT with the HBA table, and without it -------------------------------
T="${CASE_TMP}/tables"
rm -rf "$T"
mkdir -p "$T"
printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n' \
  > "${T}/${SAMPLE}_smn_copy_number.tsv"
# Made-up depths: 30x everywhere, a quarter of the gene depth at MAPQ 0
python3 "${REPO}/bin/hba_copy_number.py" bed \
  | awk -F'\t' -v OFS='\t' '{print $0, ($3 - $2) * 30}' > "${T}/bedcov.tsv"
python3 "${REPO}/bin/hba_copy_number.py" bed \
  | awk -F'\t' -v OFS='\t' '{print $0, ($3 - $2) * ($4 ~ /^HBA/ ? 22.5 : 30)}' > "${T}/bedcov_q1.tsv"
python3 "${REPO}/bin/hba_copy_number.py" check --bedcov "${T}/bedcov.tsv" --bedcov-mapq1 "${T}/bedcov_q1.tsv" \
  --out "${T}/${SAMPLE}_hba_copy_number.tsv"
cat "${T}/${SAMPLE}_hba_copy_number.tsv"
check_eq "the made-up HBA table counts four" \
  "$(awk -F'\t' '$1 == "alpha_globin_genes" {print $2}' "${T}/${SAMPLE}_hba_copy_number.tsv")" 4

NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# report_run NAME WITH: HTML_REPORT on HG002 with the SMN table, and the HBA table when WITH=1
report_run() {
  local name=$1 with=$2 out="${CASE_TMP}/$1" hba=""
  rm -rf "$out"
  [ "$with" = 1 ] && hba=", file('${T}/${SAMPLE}_hba_copy_number.tsv')"
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    HTML_REPORT(channel.of(
        [[id: '${SAMPLE}', sex: 'male'], file('${VCF}'), [file('${T}/${SAMPLE}_smn_copy_number.tsv')${hba}]]))
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
LINE='HBA1 and HBA2 together (alpha-globin genes, from depth)'
html_count() { grep -Fc "$2" "${CASE_TMP}/$1/results/${SAMPLE}/${SAMPLE}_report.html" 2>/dev/null || true; }
# hba_not_assessed NAME: 1 when summary.json lists HBA1/HBA2 as not assessed, else 0
hba_not_assessed() {
  python3 - "${CASE_TMP}/$1/results/${SAMPLE}/summary.json" <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except (OSError, ValueError):
    print("no summary"); sys.exit(0)
print(int(any(x.startswith("HBA1/HBA2 copy number") for x in d.get("not_assessed", []))))
PY
}

report_run with 1
check "HTML_REPORT wrote the report" test -s "${CASE_TMP}/with/results/${SAMPLE}/${SAMPLE}_report.html"
check_ge "the HBA line with the count of four" "$(html_count with "${LINE}</span><span class=\"value\">4 (four alpha-globin genes")" 1
check_eq "the not-assessed list leaves HBA1/HBA2 out" "$(hba_not_assessed with)" 0

report_run without 0
check "control: HTML_REPORT wrote the report without the HBA table" \
  test -s "${CASE_TMP}/without/results/${SAMPLE}/${SAMPLE}_report.html"
check_eq "control: without the HBA table there is no HBA line" "$(html_count without "$LINE")" 0
check_eq "control: without the HBA table HBA1/HBA2 are not assessed" "$(hba_not_assessed without)" 1

# --- 2. main.nf stages the HBA table in HTML_REPORT (stub run) ------------------------
# stub_staged DIR: in the stub run under DIR, "<HTML_REPORT tasks with the SMN table> <of them with the HBA table>"
stub_staged() {
  local smn=0 hba=0 d id
  while IFS= read -r d; do
    d=$(dirname "$d")
    id=$(sed -n "s/^### name: [\"'].*HTML_REPORT (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
    [ -e "${d}/input/${id}_smn_copy_number.tsv" ] || continue
    smn=$((smn + 1))
    [ ! -e "${d}/input/${id}_hba_copy_number.tsv" ] || hba=$((hba + 1))
  done < <(grep -l "^### name: [\"'].*HTML_REPORT (" -r "$1/work" --include=.command.run 2>/dev/null)
  echo "$smn $hba"
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
JOIN='BAM_ANALYSIS.out.hba_copy_number.map'
check_eq "main.nf joins the HBA table into the report inputs once" "$(grep -c "$JOIN" "${CASE_TMP}/stub/main.nf")" 1
sed -i "/${JOIN}/d" "${CASE_TMP}/stub-nojoin/main.nf"
check_eq "control copy: the join is gone" "$(grep -c "$JOIN" "${CASE_TMP}/stub-nojoin/main.nf")" 0
stub_run stub
read -r SMN HBA <<< "$(stub_staged "${CASE_TMP}/stub")"
echo "stub run: ${SMN} HTML_REPORT tasks with the SMN table, ${HBA} with the HBA table"
check_ge "HTML_REPORT tasks with the SMN table in the stub run" "$SMN" 1
check_eq "each of them gets the HBA table" "$HBA" "$SMN"
stub_run stub-nojoin
read -r SMN0 HBA0 <<< "$(stub_staged "${CASE_TMP}/stub-nojoin")"
echo "control stub run: ${SMN0} HTML_REPORT tasks with the SMN table, ${HBA0} with the HBA table"
check_ge "control: HTML_REPORT tasks with the SMN table" "$SMN0" 1
check_eq "control: without the join no task gets the HBA table" "$HBA0" 0

rm -rf "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-nojoin" "${CASE_TMP:?}/work-with" "${CASE_TMP:?}/work-without"
finish
