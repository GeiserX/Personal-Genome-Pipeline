#!/usr/bin/env bash
# Step 16b's gene callability, on both paths, with the real mosdepth image.
#   1. MAPQ 0: a copy of the fixture's chr20 slice with every read of
#      chr20:10.39-10.49 Mb set to MAPQ 0, and a two-gene table: a gene there
#      is not read, a gene over SNAP25 (unique sequence, 30x) is read. Both
#      by scripts/16b-mosdepth.sh (CALLABILITY_GENES) and by the MOSDEPTH
#      module, run alone through a two-line pipeline. Control: on the
#      unchanged BAM the same gene is read, so the MAPQ 0 is what drops it.
#   2. the reports: the base cases' step 16b and step 24 (bash) and the
#      Nextflow run's MOSDEPTH table, read by HTML_REPORT run alone, show the
#      section, and the genes it lists as not read with no owner are as many
#      as a query of assets/callability_genes.tsv and the per-gene table
#      gives (awk, not the report's code).
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "an earlier case left the fixture BAM" test -s "$BAM"

# --- 1. a gene of MAPQ 0 reads, and a unique gene ------------------------------------
CAT_DIR="${CASE_TMP}/catalog"
mkdir -p "$CAT_DIR"
CAT="${CAT_DIR}/callability_genes.tsv"
printf 'gene\thgnc_id\tchrom\tstart\tend\treasons\towner\n' > "$CAT"
printf 'TEST_UNIQUE\tHGNC:0\tchr20\t10220000\t10300000\tsegdup\t.\n' >> "$CAT"
printf 'TEST_MAPQ0\tHGNC:1\tchr20\t10400000\t10480000\tsegdup\t.\n' >> "$CAT"
CAL="${SAMPLE}-cal" RAW="${SAMPLE}-raw"
for s in "$CAL" "$RAW"; do rm -rf "${GENOME_DIR:?}/${s}"; mkdir -p "${GENOME_DIR}/${s}/aligned"; done
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$SAMTOOLS_IMAGE" sh -c "samtools view -h '${SAMPLE}/aligned/${SAMPLE}_sorted.bam' chr20 \
  | awk -F'\t' -v OFS='\t' '!/^@/ && \$3 == \"chr20\" && \$4 >= 10390000 && \$4 <= 10490000 {\$5 = 0} 1' \
  | samtools view -b -o '${CAL}/aligned/${CAL}_sorted.bam' - && samtools index '${CAL}/aligned/${CAL}_sorted.bam' \
  && samtools view -b -o '${RAW}/aligned/${RAW}_sorted.bam' '${SAMPLE}/aligned/${SAMPLE}_sorted.bam' chr20 \
  && samtools index '${RAW}/aligned/${RAW}_sorted.bam'"
check_ge "the MAPQ 0 copy has reads of MAPQ 0 in the window" \
  "$(sam view -c -q 0 "${CAL}/aligned/${CAL}_sorted.bam" chr20:10400000-10480000 2>/dev/null \
     | awk -v q="$(sam view -c -q 1 "${CAL}/aligned/${CAL}_sorted.bam" chr20:10400000-10480000 2>/dev/null)" '{print $1 - q}')" 1000

# reading SAMPLE_DIR_TABLE: "gene status" lines of a per-gene table read with the test catalog
reading() {
  python3 - "$CAT" "$1" "${REPO}/bin" <<'PY'
import csv, sys
sys.path.insert(0, sys.argv[3])
import collect_summary
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
d = collect_summary.callability(rows, sys.argv[2])
bad = {g["gene"]: g for g in d["not_read"]}
for r in rows:
    g = bad.get(r["gene"])
    print(r["gene"], "not_read" if g else "read", g["share"] if g else "")
PY
}
for s in "$CAL" "$RAW"; do
  CALLABILITY_GENES="$CAT" run_step 16b-mosdepth.sh "$s"
  check_step_exit "16b-mosdepth.sh (${s})"
done
T_CAL="${GENOME_DIR}/${CAL}/mosdepth/${CAL}_callability.bed.gz"
R_CAL=$(reading "$T_CAL" 2>&1); echo "bash, MAPQ 0 copy: ${R_CAL}"
R_RAW=$(reading "${GENOME_DIR}/${RAW}/mosdepth/${RAW}_callability.bed.gz" 2>&1); echo "bash, unchanged: ${R_RAW}"
check "bash: the gene of MAPQ 0 reads is not read" has '^TEST_MAPQ0 not_read' "$R_CAL"
check "bash: the unique gene is read" has '^TEST_UNIQUE read' "$R_CAL"
check "bash control: on the unchanged BAM the same gene is read" has '^TEST_MAPQ0 read' "$R_RAW"

NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# nf_run NAME WORKFLOW_BODY: a two-line pipeline with the real containers
nf_run() {
  local name=$1
  cat > "${NP}/main.nf" <<NF
include { MOSDEPTH } from '${REPO}/modules/local/mosdepth/main'
include { HTML_REPORT } from '${REPO}/modules/local/html_report/main'

workflow {
$2
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${CASE_TMP}/${name}'
params.publish_dir_mode = 'copy'
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  rm -rf "${CASE_TMP:?}/${name}"
  echo "+ nextflow run (${name})"
  (cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
     -work-dir "${CASE_TMP}/work-${name}") > "${CASE_TMP}/${name}.log" 2>&1
  local rc=$?
  check_eq "nextflow run (${name}) exits 0" "$rc" 0
  [ "$rc" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/${name}.log" | head -40
}
CB="${GENOME_DIR}/${CAL}/aligned/${CAL}_sorted.bam"
nf_run mosdepth "    MOSDEPTH(channel.of([[id: '${CAL}'], file('${CB}'), file('${CB}.bai')]), file('${CAT}'))"
N_CAL="${CASE_TMP}/mosdepth/${CAL}/coverage/${CAL}_callability.bed.gz"
R_NF=$(reading "$N_CAL" 2>&1); echo "Nextflow, MAPQ 0 copy: ${R_NF}"
check "Nextflow: the gene of MAPQ 0 reads is not read" has '^TEST_MAPQ0 not_read' "$R_NF"
check "Nextflow: the unique gene is read" has '^TEST_UNIQUE read' "$R_NF"
check_eq "MOSDEPTH's table is step 16b's" "$(gzip -dc "$N_CAL" 2>/dev/null | md5sum)" "$(gzip -dc "$T_CAL" 2>/dev/null | md5sum)"

# --- 2. the reports, with the real catalog -------------------------------------------
REAL="${REPO}/assets/callability_genes.tsv"
# expected TABLE: the genes with no owner that are ALT-only or under 90% at 10x, by awk
expected() {
  gzip -dc "$1" | awk -F'\t' 'NR == FNR {if (!/^#/) share[$4] = ($3 > $2 ? $5 / ($3 - $2) : 0); next}
    !/^#/ && $1 != "gene" && $7 == "." && ($4 == "." || share[$1] < 0.9) {n++} END {print n + 0}' - "$REAL"
}
# listed SUMMARY_JSON: the section's state and its count of genes not read with no owner
listed() {
  python3 -c 'import json, sys; s = json.load(open(sys.argv[1]))["sections"]["callability"]; print(s["state"], len(s["data"].get("not_read", [])))' "$1" 2>&1
}
B_T="${GENOME_DIR}/${SAMPLE}/mosdepth/${SAMPLE}_callability.bed.gz"
check "bash: step 16b wrote the per-gene table" test -s "$B_T"
check_eq "bash: its rows are the catalog genes with a position" "$(gzip -dc "$B_T" 2>/dev/null | grep -vc '^#')" \
  "$(awk -F'\t' '!/^#/ && $4 ~ /^[0-9]+$/' "$REAL" | wc -l | tr -d ' ')"
B_N=$(expected "$B_T")
check_eq "bash: step 24's summary lists the genes the catalog query gives" \
  "$(listed "${GENOME_DIR}/${SAMPLE}/summary.json")" "ok ${B_N}"
check "bash: step 24's HTML report has the section with that count" grep -qF \
  "Not read, and no other step calls them (${B_N}): " "${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"

N_T="${GENOME_DIR}/nf-results/${SAMPLE}/coverage/${SAMPLE}_callability.bed.gz"
check "Nextflow: the pipeline's MOSDEPTH published the per-gene table" test -s "$N_T"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
nf_run report "    HTML_REPORT(channel.of([[id: '${SAMPLE}', sex: 'male'], file('${VCF}'),
        [file('${N_T}'), file('${REAL}'), file('${REPO}/assets/acmg_sf.tsv'), file('${REPO}/assets/carrier_genes.tsv')]]))"
N_N=$(expected "$N_T")
check_eq "Nextflow: HTML_REPORT's summary lists the genes the catalog query gives" \
  "$(listed "${CASE_TMP}/report/${SAMPLE}/summary.json")" "ok ${N_N}"
check "Nextflow: its HTML report has the section with that count" grep -qF \
  "Not read, and no other step calls them (${N_N}): " "${CASE_TMP}/report/${SAMPLE}/${SAMPLE}_report.html"

finish
