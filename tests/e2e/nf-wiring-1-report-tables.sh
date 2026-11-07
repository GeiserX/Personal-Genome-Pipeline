#!/usr/bin/env bash
# The Nextflow report gets every table the bash report reads, and CPSR's cDNA
# read check runs in Nextflow too.
#   1. CLINVAR_CDNA_CHECK with source 'cpsr', run for real (PYPGX_IMAGE), on a
#      CPSR table holding the planted hit as a P/LP row and a benign row, with
#      the spliced reads clinvar-clip-check-1 left: it publishes
#      cpsr/<id>_cpsr_cdna_check.tsv, the P/LP row is flagged with the counts
#      bin/clip_check.py --cpsr gives step 17, and the benign row is not
#      checked. Control: no exon file, and the row is not checked.
#   2. Sample NFW: step 24 and HTML_REPORT on the same tables: the ClinVar
#      hits, a VEP file with a gnomAD AF on each hit (3% on the planted one),
#      CPSR's table and its cDNA read check, step 02's read count and the
#      clinical filter's and slivar's summary tables. The read_count, clinvar,
#      cpsr, clinical and slivar sections of the two summary.json files must
#      be equal, with the AF, the CPSR classes and the flagged row in them.
#      Control: the module as it was (no place for these tables) stops.
#   3. main.nf passes them on: a stub run of every tool (-profile test_all)
#      stages each table in each HTML_REPORT task whose sample has the step's
#      main output. A misspelt option (--outdirr) stops the run and names it.
#      Control: a copy without the new joins and without the option check
#      stages none and takes --outdirr in silence.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

S=NFW
SD="${GENOME_DIR}/${S}"
C="${SAMPLE}cdna"
CHR=$(planted chrom)
POS=$(planted pos)
REF=$(planted ref)
ALT=$(planted alt)
HG_VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
HITS_VCF="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf"
HITS_TSV="${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv"
SPLICED="${GENOME_DIR}/${C}/aligned/${C}_sorted.bam"
for f in "$HG_VCF" "$HITS_VCF" "$HITS_TSV" "$SPLICED" "${SPLICED}.bai" \
         "${GENOME_DIR}/${SAMPLE}/qc/${SAMPLE}_read_count.tsv" \
         "${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv" \
         "${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" \
         "${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_slivar_summary.tsv" \
         "${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz"; do
  check "an earlier case left ${f#"$GENOME_DIR"/}" test -s "$f"
done

# A small pipeline in its own project folder, which gets the repository's bin/
NP="${CASE_TMP}/project"
rm -rf "$NP"
mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
# nf_run NAME: run NP/main.nf with the real containers, output under CASE_TMP/NAME
nf_run() {
  local name=$1 out="${CASE_TMP}/$1"
  rm -rf "$out"
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
  NF_RC=$?
}

# --- 1. CLINVAR_CDNA_CHECK on CPSR's rows -----------------------------------------
mkdir -p "${CASE_TMP}/in"
CLS="${CASE_TMP}/in/${S}.cpsr.grch38.classification.tsv.gz"
{
  printf 'VAR_ID\tGENOMIC_CHANGE\tSYMBOL\tCLINVAR_CLASSIFICATION\tCPSR_CLASSIFICATION\tCLASSIFICATION\n'
  printf '%s_%s_%s_%s\t%s:g.%s%s>%s\t%s\tPathogenic\tPathogenic\tPathogenic\n' "${CHR#chr}" "$POS" "$REF" "$ALT" \
    "${CHR#chr}" "$POS" "$REF" "$ALT" "$(planted gene)"
  printf '%s_%s_A_G\t%s:g.%sA>G\tOTHER\tBenign\tBenign\tBenign\n' "${CHR#chr}" "$((POS + 500))" "${CHR#chr}" "$((POS + 500))"
} | gzip -c > "$CLS"
# The exon file of clinvar-cpsr-quality-1: the planted hit's exon and the next one.
EXONS="${CASE_TMP}/in/exons.bed"
P0=$((POS - 1))
{
  printf '%s\t%s\t%s\n' "$CHR" "$((P0 - 30))" "$((P0 + 30))"
  printf '%s\t%s\t%s\n' "$CHR" "$((P0 + 3000))" "$((P0 + 3150))"
} > "$EXONS"
cdna_nf() {  # cdna_nf NAME EXONS-EXPRESSION
  cat > "${NP}/main.nf" <<NF
include { CLINVAR_CDNA_CHECK } from '${REPO}/modules/local/clinvar_screen/main'

workflow {
    CLINVAR_CDNA_CHECK(channel.of([[id: '${S}'], file('${CLS}'), file('${SPLICED}'), file('${SPLICED}.bai')]),
                       $2, [], 'cpsr')
}
NF
  nf_run "$1"
  check_eq "nextflow run ($1) exits 0" "$NF_RC" 0
  [ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/$1.log" | head -40
}
cdna_nf cdna "channel.value(file('${EXONS}'))"
NF_CDNA="${CASE_TMP}/cdna/results/${S}/cpsr/${S}_cpsr_cdna_check.tsv"
check "the table is published as cpsr/${S}_cpsr_cdna_check.tsv" test -s "$NF_CDNA"
check "and no ClinVar scan table beside it" test ! -e "${CASE_TMP}/cdna/results/${S}/cpsr/${S}_clinvar_cdna_scan.tsv"
# Step 17's call of bin/clip_check.py --cpsr on the same files
docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome:ro" -v "${CASE_TMP}:/t" -v "${REPO}/bin:/pgp-bin:ro" \
  "$PYPGX_IMAGE" python3 /pgp-bin/clip_check.py --cpsr "/t/in/${S}.cpsr.grch38.classification.tsv.gz" \
  --bam "/genome/${C}/aligned/${C}_sorted.bam" --exons /t/in/exons.bed --out /t/in/step17.tsv > /dev/null
ROW=$(awk -F'\t' -v p="$POS" 'NR > 1 && $2 == p' "$NF_CDNA" 2>/dev/null)
echo "Nextflow CPSR row: ${ROW}"
check_eq "only the P/LP row is in the table" "$(awk 'NR > 1' "$NF_CDNA" 2>/dev/null | wc -l | tr -d ' ')" 1
check_eq "the planted P/LP row is flagged" "$(cut -f9 <<< "$ROW")" flagged
check "the table is the one bin/clip_check.py --cpsr writes for step 17, byte for byte" cmp -s "${CASE_TMP}/in/step17.tsv" "$NF_CDNA"
cdna_nf cdna-no-exons "[]"
OFF=$(awk -F'\t' -v p="$POS" 'NR > 1 && $2 == p' "${CASE_TMP}/cdna-no-exons/results/${S}/cpsr/${S}_cpsr_cdna_check.tsv" 2>/dev/null)
check_eq "control: without exons the row is not checked, so the flagged check can fail" "$(cut -f9 <<< "$OFF")" not_checked

# --- 2. step 24 and HTML_REPORT on the same tables ---------------------------------
rm -rf "$SD"
mkdir -p "${SD}/vcf" "${SD}/clinvar" "${SD}/vep" "${SD}/cpsr" "${SD}/qc" "${SD}/clinical" "${SD}/slivar"
ln -f "$HG_VCF" "${SD}/vcf/${S}.vcf.gz"
ln -f "${HG_VCF}.tbi" "${SD}/vcf/${S}.vcf.gz.tbi"
cp "$HITS_VCF" "${SD}/clinvar/${S}_clinvar_hits.vcf"
cp "$CLS" "${SD}/cpsr/"
cp "$NF_CDNA" "${SD}/cpsr/"
cp "${GENOME_DIR}/${SAMPLE}/qc/${SAMPLE}_read_count.tsv" "${SD}/qc/${S}_read_count.tsv"
cp "${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv" "${SD}/clinical/${S}_clinical_summary.tsv"
cp "${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" "${SD}/clinical/${S}_clinical.vcf.gz"
cp "${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_slivar_summary.tsv" "${SD}/slivar/${S}_slivar_summary.tsv"
cp "${GENOME_DIR}/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" "${SD}/slivar/${S}_prioritized.vcf.gz"
# A VEP file with one record per hit: a gnomAD AF on each, 3% on the planted one.
python3 - "$HITS_TSV" "${SD}/vep/plain.vcf" "$CHR" "$POS" <<'PY'
import sys
hits, out, pchrom, ppos = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
rows = [l.rstrip("\n").split("\t") for l in open(hits)][1:]
recs = {}
for i, r in enumerate(rows):
    chrom, pos, ref, alt = r[0], int(r[1]), r[2], r[3]
    af = 0.03 if (chrom, pos) == (pchrom, ppos) else round(0.0001 * (i + 1), 6)
    snv = len(ref) == 1 and len(alt) == 1
    allele = alt if snv or alt[:1] != ref[:1] else (alt[1:] or "-")
    recs[(chrom, pos, ref, alt)] = f"{chrom}\t{pos}\trs{1000 + i}\t{ref}\t{alt}\t50\tPASS\tCSQ={allele}|x|LOW|G|{af}|{af / 2}|rs{1000 + i}\tGT\t0/1"
order = {c: i for i, c in enumerate(sorted({k[0] for k in recs}, key=lambda c: (len(c), c)))}
with open(out, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for c in order:
        f.write(f"##contig=<ID={c}>\n")
    f.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            'Format: Allele|Consequence|IMPACT|SYMBOL|gnomADe_AF|gnomADg_AF|Existing_variation">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
    for k in sorted(recs, key=lambda k: (order[k[0]], k[1])):
        f.write(recs[k] + "\n")
print(f"{len(recs)} VEP records")
PY
in_user() { docker run --rm -u "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -w /genome "$BCFTOOLS_IMAGE" "$@"; }
in_user bcftools view --no-version -Oz -o "${S}/vep/${S}_vep.vcf.gz" "${S}/vep/plain.vcf"
in_user bcftools index -t "${S}/vep/${S}_vep.vcf.gz"
rm -f "${SD}/vep/plain.vcf"
check "the VEP file is bgzipped and indexed" test -s "${SD}/vep/${S}_vep.vcf.gz.tbi"
run_step 24-html-report.sh "$S"
check_step_exit 24-html-report.sh

# The report assets main.nf passes HTML_REPORT (report_assets)
ASSET_FILES=$(for n in acmg_sf.tsv carrier_genes.tsv callability_genes.tsv mito_confirmed.tsv paralog_loci.tsv paraphase_genes.tsv; do
  printf "file('%s/assets/%s'), " "$REPO" "$n"; done | sed 's/, $//')
report_nf() {  # report_nf NAME MODULE
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '$2'

workflow {
    HTML_REPORT(channel.of([[id: '${S}', sex: 'male'], file('${SD}/vcf/${S}.vcf.gz'), [
        file('${SD}/clinvar/${S}_clinvar_hits.vcf'),
        file('${SD}/vep/${S}_vep.vcf.gz'), file('${SD}/vep/${S}_vep.vcf.gz.tbi'),
        file('${SD}/cpsr/${S}.cpsr.grch38.classification.tsv.gz'), file('${SD}/cpsr/${S}_cpsr_cdna_check.tsv'),
        file('${SD}/qc/${S}_read_count.tsv'),
        file('${SD}/clinical/${S}_clinical.vcf.gz'), file('${SD}/clinical/${S}_clinical_summary.tsv'),
        file('${SD}/slivar/${S}_prioritized.vcf.gz'), file('${SD}/slivar/${S}_slivar_summary.tsv'),
        ${ASSET_FILES}]]))
}
NF
  nf_run "$1"
}
report_nf report "${REPO}/modules/local/html_report/main"
check_eq "nextflow run (report) exits 0" "$NF_RC" 0
[ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/report.log" | head -40
NF_SUM="${CASE_TMP}/report/results/${S}/summary.json"
check "HTML_REPORT wrote the report" test -s "${CASE_TMP}/report/results/${S}/${S}_report.html"
# same_sections BASH NF: the five sections of the two summaries, equal but for the file times
same_sections() {
  python3 - "$1" "$2" <<'PY'
import json, sys
b, n = (json.load(open(p))["sections"] for p in sys.argv[1:3])
bad = 0
for k in ("read_count", "clinvar", "cpsr", "clinical", "slivar"):
    x, y = ({i: v for i, v in s.get(k, {"state": "absent"}).items() if i != "file_date"} for s in (b, n))
    ok = x == y and x["state"] == "ok"
    print(f"  {k}: bash {x['state']}, nextflow {y['state']}{'' if ok else ', differ'}")
    if not ok:
        for i in sorted(set(x) | set(y)):
            if x.get(i) != y.get(i):
                print(f"    {i}: bash {json.dumps(x.get(i))[:400]}")
                print(f"    {i}: nextflow {json.dumps(y.get(i))[:400]}")
    bad += not ok
sys.exit(bad)
PY
}
check "the five sections of the Nextflow summary equal the bash ones" same_sections "${SD}/summary.json" "$NF_SUM"
nf_py() { python3 -c "import json, sys; d = json.load(open(sys.argv[1]))['sections']; print($1)" "$NF_SUM" 2>/dev/null; }
check_eq "the ClinVar hits' frequencies come from the VEP file" \
  "$(nf_py "d['clinvar']['data']['frequency']['source']")" "vep/${S}_vep.vcf.gz"
check_eq "the planted hit has its 3% AF and is marked common" \
  "$(nf_py "[(h['af'], h['common']) for h in d['clinvar']['data']['hits'] if h['chrom'] == '${CHR}' and h['pos'] == ${POS}]")" \
  "[(0.03, True)]"
check_eq "the hits carry VEP's rsIDs" \
  "$(nf_py "sum(1 for h in d['clinvar']['data']['hits'] if h['rsids'])")" "$(awk 'NR > 1' "$HITS_TSV" | wc -l | tr -d ' ')"
check_eq "the CPSR card counts the classes" "$(nf_py "d['cpsr']['data']['classification']")" "{'Pathogenic': 1, 'Benign': 1}"
check_eq "and flags the P/LP row the read check flagged" "$(nf_py "d['cpsr']['data']['cdna_flagged']")" 1
check_eq "the read count section is there" "$(nf_py "d['read_count']['data']['read_count_check']")" \
  "$(awk -F'\t' '$1 == "read_count_check" {print $2}' "${SD}/qc/${S}_read_count.tsv")"
check_eq "the clinical card counts the summary table's rows" "$(nf_py "d['clinical']['data']['variants']")" \
  "$(awk 'NR > 1' "${SD}/clinical/${S}_clinical_summary.tsv" | wc -l | tr -d ' ')"
check "and has its impact breakdown" test -n "$(nf_py "d['clinical']['data']['by_impact']")"
check_eq "the slivar card counts the summary table's rows" "$(nf_py "d['slivar']['data']['prioritized']")" \
  "$(awk 'NR > 1' "${SD}/slivar/${S}_slivar_summary.tsv" | wc -l | tr -d ' ')"
# Control: the module as it was, without a place for the new tables.
OLD="${CASE_TMP}/old-module/main.nf"
mkdir -p "$(dirname "$OLD")"
sed -e '/_vep\.vcf\.gz|/d' -e 's/|\${id}_clinical_summary\.tsv//' -e 's/|\${id}_slivar_summary\.tsv//' \
    -e 's/|\${id}\.cpsr\.grch38\.classification\.tsv\.gz|\${id}_cpsr_cdna_check\.tsv//' \
    -e 's/|\${id}_read_count\.tsv//' "${REPO}/modules/local/html_report/main.nf" > "$OLD"
check_eq "control copy: the module names none of the new tables" \
  "$(grep -cE 'vep\.vcf|clinical_summary|slivar_summary|classification\.tsv|cpsr_cdna_check|read_count' "$OLD")" 0
report_nf report-old "${OLD%.nf}"
check "control: the module without the places stops on the first new table" test "$NF_RC" -ne 0
check "control: and says it has no place for it" grep -q 'HTML_REPORT got an input it has no place for' -r "${CASE_TMP}/work-report-old" --include=.command.err

# --- 3. main.nf passes every table on (stub run of every tool) ---------------------
# PAIRS: each new table and the step output (or published file, with "out:")
# that says the step ran for the sample
PAIRS="_vep.vcf.gz:_clinical.vcf.gz _vep.vcf.gz.tbi:_clinical.vcf.gz _clinical_summary.tsv:_clinical.vcf.gz \
_slivar_summary.tsv:_prioritized.vcf.gz .cpsr.grch38.classification.tsv.gz:.cpsr.grch38.html \
_cpsr_cdna_check.tsv:.cpsr.grch38.html _read_count.tsv:out:qc"
# stub_staged DIR: one line per pair, "<new> <HTML_REPORT tasks with the step output> <of them with the table>"
stub_staged() {
  local d id p new old hits got
  for p in $PAIRS; do
    new=${p%%:*} old=${p#*:} hits=0 got=0
    while IFS= read -r d; do
      d=$(dirname "$d")
      id=$(sed -n "s/^### name: [\"'].*HTML_REPORT (\(.*\))[\"']\$/\1/p" "${d}/.command.run")
      case "$old" in
        out:*) [ -e "$1/out/${id}/${old#out:}/${id}${new}" ] || continue ;;
        *) [ -e "${d}/input/${id}${old}" ] || continue ;;
      esac
      hits=$((hits + 1))
      [ ! -e "${d}/input/${id}${new}" ] || got=$((got + 1))
    done < <(grep -l "^### name: [\"'].*HTML_REPORT (" -r "$1/work" --include=.command.run 2>/dev/null)
    echo "$new $hits $got"
  done
}
stub_run() {  # stub_run NAME [OPTION...]
  local d="${CASE_TMP}/$1"; shift
  echo "+ nextflow run main.nf -profile test_all -stub ($(basename "$d") $*)"
  (cd "$d" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${d}/out" -work-dir "${d}/work" "$@") \
    > "${d}.log" 2>&1
}
for n in stub stub-control; do
  rm -rf "${CASE_TMP:?}/${n}"
  mkdir -p "${CASE_TMP}/${n}"
  cp -R "${REPO}/main.nf" "${REPO}/nextflow.config" "${REPO}/nextflow_schema.json" "${REPO}/versions.env" \
    "${REPO}/conf" "${REPO}/modules" "${REPO}/workflows" "${REPO}/bin" "${REPO}/assets" "${REPO}/scripts" "${CASE_TMP}/${n}/"
  mkdir -p "${CASE_TMP}/${n}/tests/fixtures"  # test_all's consumer_sites
  cp -R "${REPO}/tests/fixtures/consumer_export" "${CASE_TMP}/${n}/tests/fixtures/"
done
JOINS='ANNOTATION.out.vep_vcf.map|ANNOTATION.out.vep_vcf_index.map|ANNOTATION.out.clinical_summary.map|ANNOTATION.out.slivar_summary.map|CLINICAL.out.cpsr_tsv.map|CLINICAL.out.cpsr_cdna.map|UPSTREAM.out.read_count.map'
check_eq "main.nf joins the seven tables into the report inputs" "$(grep -cE "$JOINS" "${CASE_TMP}/stub/main.nf")" 7
sed -i -E "/${JOINS}/d" "${CASE_TMP}/stub-control/main.nf"
sed -i 's/^    if (unknown_params) {$/    if (false) {/' "${CASE_TMP}/stub-control/main.nf"
check_eq "control copy: the joins are gone" "$(grep -cE "$JOINS" "${CASE_TMP}/stub-control/main.nf")" 0
check_eq "control copy: the option check is off" "$(grep -c '^    if (false) {$' "${CASE_TMP}/stub-control/main.nf")" 1
stub_run stub
check_eq "stub run exits 0" "$?" 0
tail -n 5 "${CASE_TMP}/stub.log"
stub_run stub-control --outdirr x
check_eq "control: the stub run without the option check takes --outdirr and exits 0" "$?" 0
while read -r new hits got; do
  check_ge "stub run: tasks with the step before ${new}" "$hits" 1
  check_eq "stub run: each of them gets ${new}" "$got" "$hits"
done < <(stub_staged "${CASE_TMP}/stub")
while read -r new hits got; do
  check_ge "control: tasks with the step before ${new}" "$hits" 1
  check_eq "control: none gets ${new}" "$got" 0
done < <(stub_staged "${CASE_TMP}/stub-control")
check_ge "stub run: CPSR's read check ran for a sample" \
  "$(find "${CASE_TMP}/stub/out" -path '*/cpsr/*_cpsr_cdna_check.tsv' | wc -l | tr -d ' ')" 1
rm -rf "${CASE_TMP:?}/stub/work" "${CASE_TMP:?}/stub/out"
stub_run stub --outdirr x --max-cpus 2
check "a misspelt option stops the run" test "$?" -ne 0
check "and the message names --outdirr" grep -q 'unknown parameters .*--outdirr' "${CASE_TMP}/stub.log"
check "and the kebab-case --max-cpus" grep -Eq 'unknown parameters .*--(max-cpus|maxCpus)' "${CASE_TMP}/stub.log"
grep -m1 'unknown parameter' "${CASE_TMP}/stub.log"

rm -rf "$SD" "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-control" "${CASE_TMP:?}"/work-*
finish
