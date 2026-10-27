#!/usr/bin/env bash
# The ACMG SF rules (assets/acmg_sf.tsv) and the carrier list
# (assets/carrier_genes.tsv) reach every report, bash and Nextflow. Without
# the assets the reports still render and only say "not assessed" or "not
# made", so a missing mount or staging passes every other case.
#   1. Bash: the text report and summary.json that generate-report.sh left
#      (case 51, after step 24) and the HTML of step 24 (case 50) have the
#      SF line with the 84 genes and the carrier list, and no "not assessed"
#      or "not made" line. Control: the same render without the assets
#      mount (what generate-report.sh did before it mounted assets/) fails
#      the same checks.
#   2. Nextflow: HTML_REPORT on the sample's ClinVar hits and clinical VCF,
#      with the two assets staged as main.nf stages them, applies the rules,
#      and its carrier list equals the bash one. Control: the same task
#      without the assets says "unavailable". And main.nf itself, run with
#      ClinVar and the report, made the carrier list, so main.nf stages them.
. "$(dirname "$0")/lib.sh"

SD="${GENOME_DIR}/${SAMPLE}"
TXT="${SD}/${SAMPLE}_report.txt"
HTML="${SD}/${SAMPLE}_report.html"
SUM="${SD}/summary.json"
for f in "$TXT" "$HTML" "$SUM" "${SD}/clinvar/${SAMPLE}_clinvar_hits.vcf" "${SD}/clinical/${SAMPLE}_clinical.vcf.gz"; do
  check "an earlier case left ${f#"$GENOME_DIR"/}" test -s "$f"
done

# report_ok FILE: the SF line with its 84 genes, the carrier list, and no
# line saying either could not be made.
report_ok() {
  local r
  r=$(cat "$1" 2>/dev/null)
  has 'Secondary-findings genes \(ACMG SF v3\.3, 84 genes\)' "$r" &&
    has 'Carrier list \([0-9]+ genes\)' "$r" &&
    lacks 'not assessed, acmg_sf\.tsv' "$r" && lacks 'Carrier list: not made' "$r"
}
# summary_ok FILE: the clinical section's SF tier has the 84 rules, and the
# carrier section is ok with its gene list; neither is "unavailable".
summary_ok() {
  python3 - "$1" <<'PY'
import json, sys
try:
    s = json.load(open(sys.argv[1]))["sections"]
    sf = s["clinical"]["data"]["acmg_sf"]
    c = s["carrier"]
except (OSError, ValueError, KeyError, TypeError) as e:
    sys.exit(f"  unreadable: {e!r}")
print(f"  acmg_sf: genes_on_list {sf.get('genes_on_list')}, unavailable {sf.get('unavailable')!r}")
print(f"  carrier: state {c['state']}, genes_on_list {c['data'].get('genes_on_list')}, "
      f"unavailable {c['data'].get('unavailable')!r}")
ok = (sf.get("genes_on_list") == 84 and "unavailable" not in sf and c["state"] == "ok"
      and (c["data"].get("genes_on_list") or 0) > 0 and "unavailable" not in c["data"])
sys.exit(0 if ok else 1)
PY
}

# --- 1. bash ----------------------------------------------------------------------
check "generate-report.sh's text report applies the SF rules and has the carrier list" report_ok "$TXT"
check "generate-report.sh's summary.json has the 84 SF rules and the carrier list" summary_ok "$SUM"
H=$(cat "$HTML" 2>/dev/null)
check "step 24's HTML has the carrier list card" has 'Carrier List \(this sample\)' "$H"
check "step 24's HTML has no 'not assessed' or 'not made' line" \
  lacks 'not assessed, acmg_sf\.tsv|Carrier list: not made' "$H"

CTL="${CASE_TMP}/control"
rm -rf "$CTL"
mkdir -p "$CTL"
docker run --rm --network none --user "$(id -u):$(id -g)" \
  -v "${REPO}/bin:/pgp-bin:ro" -v "${GENOME_DIR}:/genome:ro" -v "${CTL}:/out" "$PYTHON_IMAGE" \
  python3 /pgp-bin/render_report.py --sample "$SAMPLE" --sample-dir "/genome/${SAMPLE}" \
    --json /out/summary.json -o /out/report.txt > "${CTL}/render.log" 2>&1
check_eq "control: the render without assets/ exits 0" "$?" 0
check "control: its text report is written" test -s "${CTL}/report.txt"
if [ -s "${CTL}/report.txt" ] && ! report_ok "${CTL}/report.txt"; then CT=red; else CT=green; fi
check_eq "control: the text check fails without assets/" "$CT" red
if [ -s "${CTL}/summary.json" ] && ! summary_ok "${CTL}/summary.json" > /dev/null; then CS=red; else CS=green; fi
check_eq "control: the summary check fails without assets/" "$CS" red

# --- 2. Nextflow --------------------------------------------------------------------
if ! command -v nextflow >/dev/null; then fail "nextflow is not on PATH"; finish; fi
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# nf_report NAME ASSETS: HTML_REPORT on the sample's ClinVar hits and clinical
# VCF, plus the two assets when ASSETS is 1; output under CASE_TMP/NAME
nf_report() {
  local name=$1 out="${CASE_TMP}/$1" extra="" rc
  rm -rf "$out"
  [ "$2" = 1 ] && extra=", file('${REPO}/assets/acmg_sf.tsv'), file('${REPO}/assets/carrier_genes.tsv')"
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
    HTML_REPORT(channel.of([[id: '${SAMPLE}', sex: 'male'], file('${SD}/vcf/${SAMPLE}.vcf.gz'), [
        file('${SD}/clinvar/${SAMPLE}_clinvar_hits.vcf'), file('${SD}/clinical/${SAMPLE}_clinical.vcf.gz')${extra}]]))
}
NF
  {
    echo "params.outdir = '${out}/results'"
    echo "params.publish_dir_mode = 'copy'"
    echo "docker.enabled = true"
    echo "docker.runOptions = '--user $(id -u):$(id -g) --network none'"
  } > "${NP}/nextflow.config"
  echo "+ nextflow run (${name})"
  (
    cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
      -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work-${name}"
  ) > "${CASE_TMP}/${name}.log" 2>&1
  rc=$?
  check_eq "nextflow run (${name}) exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/${name}.log" | head -40
}

nf_report report 1
NF_SUM="${CASE_TMP}/report/results/${SAMPLE}/summary.json"
check "HTML_REPORT with the assets applies the 84 SF rules and makes the carrier list" summary_ok "$NF_SUM"
check "HTML_REPORT's carrier list equals the bash one" python3 - "$SUM" "$NF_SUM" <<'PY'
import json, sys
try:
    b, n = (json.load(open(p))["sections"]["carrier"]["data"] for p in sys.argv[1:3])
except (OSError, ValueError, KeyError) as e:
    sys.exit(f"  unreadable: {e!r}")
print(f"  bash {len(b.get('carriers', []))} carrier(s), nextflow {len(n.get('carriers', []))}")
sys.exit(0 if b == n else 1)
PY
NH=$(cat "${CASE_TMP}/report/results/${SAMPLE}/${SAMPLE}_report.html" 2>/dev/null)
check "the Nextflow HTML has the carrier list card" has 'Carrier List \(this sample\)' "$NH"

# main.nf itself on a VCF-only sheet with ClinVar and the report: main.nf
# stages the assets for HTML_REPORT, so its carrier list is made.
MN="${CASE_TMP}/mainnf"
rm -rf "$MN"
mkdir -p "$MN/run"
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" "${SD}/vcf/${SAMPLE}.vcf.gz" "${SD}/vcf/${SAMPLE}.vcf.gz.tbi" \
  > "${MN}/sheet.csv"
echo "+ nextflow run main.nf --tools clinvar,html_report"
(
  cd "${MN}/run" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${CASE_TMP}/work-mainnf" \
    --input "${MN}/sheet.csv" --reference "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" \
    --tools clinvar,html_report \
    --clinvar "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
    --clinvar_index "${GENOME_DIR}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
    --outdir "${MN}/results" --max_cpus 4 --max_memory 14.GB
) > "${MN}/nextflow.log" 2>&1
MRC=$?
check_eq "nextflow run main.nf (clinvar,html_report) exits 0" "$MRC" 0
[ "$MRC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${MN}/nextflow.log" | head -40
# No clinical VCF in this run, so no clinical section: the carrier list
# shows the staging (main.nf passes both assets in one list).
check "main.nf's report made the carrier list (main.nf stages the assets)" \
  python3 - "${MN}/results/${SAMPLE}/summary.json" <<'PY2'
import json, sys
try:
    c = json.load(open(sys.argv[1]))["sections"]["carrier"]
except (OSError, ValueError, KeyError) as e:
    sys.exit(f"  unreadable: {e!r}")
print(f"  carrier: state {c['state']}, genes_on_list {c['data'].get('genes_on_list')}, "
      f"unavailable {c['data'].get('unavailable')!r}")
sys.exit(0 if c["state"] == "ok" and (c["data"].get("genes_on_list") or 0) > 0
         and "unavailable" not in c["data"] else 1)
PY2

nf_report no-assets 0
NA_SUM="${CASE_TMP}/no-assets/results/${SAMPLE}/summary.json"
check "control: HTML_REPORT without the assets still writes the summary" test -s "$NA_SUM"
if [ -s "$NA_SUM" ] && ! summary_ok "$NA_SUM" > /dev/null; then CN=red; else CN=green; fi
check_eq "control: the summary check fails without the staged assets" "$CN" red

finish
