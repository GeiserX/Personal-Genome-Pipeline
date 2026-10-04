#!/usr/bin/env bash
# sample_qc in the Nextflow pipeline (SOMALIER, SOMALIER_RELATE, VERIFYBAMID2,
# SAMPLE_QC) on the VCF+BAM row of cases 20 and 21, with the sites case qc-1
# wrote (somalier's plus the slice's own chrX calls) and the panel setup.sh
# installed there:
#   - declared female: INDEXCOV reads the slices as female (case 34) and lets
#     the row through, but somalier finds HG002 male from the reads and stops
#     the run, naming both sexes; the report is never written;
#   - declared male with --sex_check warn (for INDEXCOV's reading): the run
#     ends, and the report, rendered by bin/render_report.py, shows somalier's
#     sex, FREEMIX and the ROH and haplogroup cards, with the numbers of the
#     summary it was rendered from.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

G="$GENOME_DIR"
XS="${G}/reference/somalier/sites_slice_chrX.vcf"
PANEL_DIR="${G}/reference/verifybamid2"
check "case qc-1 wrote the sites with the slice's chrX calls" test -s "$XS"

# nf_run NAME SEX [ARGS...]: one run on a VCF+BAM row declaring SEX. Sets RC,
# LOG, OUT and TRACE.
nf_run() {
  local name=$1 sex=$2 dir="${CASE_TMP}/$1"
  shift 2
  rm -rf "$dir"
  mkdir -p "$dir"
  printf 'sample,vcf,vcf_index,bam,bam_index,sex\n%s,%s,%s,%s,%s,%s\n' "$SAMPLE" \
    "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" \
    "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${G}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" \
    "$sex" > "${dir}/samplesheet.csv"
  OUT="${dir}/out"
  echo "+ nextflow run main.nf (${name}, declared ${sex}) $*"
  ( cd "$dir" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${dir}/work" \
      --input "${dir}/samplesheet.csv" \
      --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" \
      --tools sample_qc,roh,mito_haplogroup,mosdepth,html_report \
      --somalier_sites "$XS" --verifybamid2_panel "$PANEL_DIR" \
      --outdir "$OUT" --max_cpus 4 --max_memory 14.GB "$@" ) > "${dir}/run.log" 2>&1
  RC=$?
  LOG="${dir}/run.log"
  grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "$LOG"
  cp "${dir}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.${name}.nextflow.log" 2>/dev/null
  TRACE=$(find "${OUT}/pipeline_info" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
}
# ran NAME: tasks of process NAME in the trace, in any state
ran() {
  awk -F'\t' -v p="$1" 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    { n = $c["name"]; sub(/ \(.*/, "", n); sub(/.*:/, "", n); if (n == p) k++ }
    END {print k + 0}' "${TRACE:-/dev/null}" 2>/dev/null
}

nf_run female female
check "declared female: the run stops (exit ${RC})" test "$RC" -ne 0
check "the message names the sample and both sexes" \
  has "Sample '${SAMPLE}': the samplesheet says sex female, but somalier infers male from the reads" "$(cat "$LOG")"
check "INDEXCOV let the row through first" has "Sample '${SAMPLE}': indexcov infers female" "$(cat "$LOG")"
check "the message says how to go on" has 'rerun with --sex_check warn' "$(cat "$LOG")"
check_eq "declared female: SAMPLE_QC tasks" "$(ran SAMPLE_QC)" 1
check_eq "declared female: no HTML_REPORT task" "$(ran HTML_REPORT)" 0

nf_run male male --sex_check warn
check_eq "declared male, --sex_check warn: the run exits 0" "$RC" 0
for p in SOMALIER SOMALIER_RELATE VERIFYBAMID2 SAMPLE_QC HTML_REPORT; do
  check_eq "${p} tasks" "$(ran "$p")" 1
done
R="${OUT}/${SAMPLE}"
T="${R}/qc/${SAMPLE}_sample_qc.tsv"
cat "$T" 2>/dev/null
check_eq "SAMPLE_QC: somalier infers male" "$(awk -F'\t' '$1 == "inferred_sex" {print $2}' "$T" 2>/dev/null)" male
check_eq "SAMPLE_QC: the sex check passes" "$(awk -F'\t' '$1 == "sex_check" {print $2}' "$T" 2>/dev/null)" ok
check "VERIFYBAMID2 published its selfSM" test -s "${R}/qc/verifybamid2/${SAMPLE}.selfSM"
check "SOMALIER_RELATE published its tables" test -s "${OUT}/somalier/somalier.samples.tsv"

HTML="${R}/${SAMPLE}_report.html"
JSON="${R}/${SAMPLE}_summary.json"
check "the report exists" test -s "$HTML"
check "summary.json validates against the schema" \
  python3 "${REPO}/tests/schema/validate.py" "${REPO}/tests/schema/summary.schema.json" "$JSON"
FREEMIX=$(python3 -c "import json; print(json.load(open('${JSON}'))['sections']['sample_qc']['data']['freemix'])" 2>/dev/null)
check "the summary has FREEMIX (${FREEMIX:-none})" test -n "$FREEMIX"
grep -o '<div class="stat"><span class="label">\(Sex from the reads (somalier)\|Contamination (FREEMIX)\)</span>.*</div>' "$HTML" 2>/dev/null
check "the QC card shows somalier's sex" \
  has 'Sex from the reads \(somalier\)</span><span class="value"><span class="badge badge-green">male<' "$(cat "$HTML" 2>/dev/null)"
check "the QC card shows FREEMIX" \
  has "Contamination \\(FREEMIX\\)</span><span class=\"value\"><span class=\"badge badge-green\">${FREEMIX} \\(warning above 0.03" "$(cat "$HTML" 2>/dev/null)"
check "the QC card shows mosdepth's depth" has 'Mean depth</span><span class="value">[0-9.]+x<' "$(cat "$HTML" 2>/dev/null)"
check "the ROH card is filled" has '<h2>Runs of Homozygosity</h2>' "$(cat "$HTML" 2>/dev/null)"
check "the haplogroup card is filled" has '<h2>Mitochondrial Haplogroup</h2>' "$(cat "$HTML" 2>/dev/null)"
check "the Variant Calling card counts the VCF" has 'Total variants</span><span class="value">[1-9][0-9]*<' "$(cat "$HTML" 2>/dev/null)"

finish
