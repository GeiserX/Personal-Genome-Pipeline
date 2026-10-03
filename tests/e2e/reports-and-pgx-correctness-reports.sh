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
#      slivar result is marked stale in both reports, and step 06's is not.
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

finish
