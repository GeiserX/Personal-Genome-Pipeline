#!/usr/bin/env bash
# The disease-gene contraindications of step 27 and CPIC_LOOKUP
# (assets/pgx_contraindications.tsv), in the python image. Everything is
# synthetic (no one's variant or medications):
#   - step 27 on a synthetic PharmCAT report and a synthetic step 06 hits
#     table with a pathogenic POLG hit: the report gives the valproic acid
#     WARNING, marked when valproate is on MEDICATIONS; fluconazole on the list
#     gives CYP2C9 a note in the fifth column;
#   - the same sample without the hits table lists the genes as not checked
#     and gives no valproic acid line (the bash control);
#   - Nextflow with --tools pharmcat,clinvar,cpic on the HG002 VCF and a
#     ClinVar file of one record, the fixture's planted SNV (HG002 is
#     homozygous for it) labelled POLG: CPIC_LOOKUP joins CLINVAR_SCREEN's
#     hits and gives the valproic acid WARNING;
#   - control: a copy of the pipeline whose main.nf passes [] for the
#     contraindication table gives no such line, so the check above depends
#     on main.nf's ch_contraindications line.
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
S2="${SAMPLE}contra"
rm -rf "${G:?}/${S2}"
mkdir -p "${G}/${S2}/pharmcat" "${G}/${S2}/clinvar"
cat > "${G}/${S2}/pharmcat/${S2}.report.json" <<'JSON'
{"pharmcatVersion": "3.4.0",
 "genes": {
  "CYP2C9": {"geneSymbol": "CYP2C9", "callSource": "MATCHER",
             "sourceDiplotypes": [{"label": "*1/*1", "allele1": {"name": "*1"}, "allele2": {"name": "*1"},
                                   "phenotypes": ["Normal Metabolizer"]}]}},
 "drugs": {}}
JSON
HITS="${G}/${S2}/clinvar/${S2}_clinvar_hits.tsv"
printf 'chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n' > "$HITS"
printf 'chr15\t89300000\tC\tT\t0/1\t9900001\tPOLG:5428\tPathogenic\tcriteria_provided,_single_submitter\n' >> "$HITS"
MED="${CASE_TMP}/medications.txt"
printf 'valproate 500 mg\nfluconazole\n' > "$MED"
REC="${G}/${S2}/cpic/${S2}_cpic_recommendations.txt"
PHENO="${G}/${S2}/cpic/${S2}_phenotypes.tsv"

# --- bash step 27 --------------------------------------------------------------
MEDICATIONS="$MED" run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
OUT=$(cat "$REC" 2>/dev/null)
check "step 27: a pathogenic POLG hit gives the valproic acid WARNING, marked as on the list" \
  has 'WARNING: valproic acid -- POLG chr15:89300000 C>T, genotype 0/1, ClinVar Pathogenic; ON THE MEDICATION LIST\.' "$OUT"
check "step 27: the line carries the label's note and its CPIC level" \
  has 'mitochondrial disorders caused by POLG mutations.*\(CPIC level A/B\)' "$OUT"
check "step 27: the step's log names the hits table" grep -q 'ClinVar hits (step 06): .*_clinvar_hits\.tsv' "$STEP_LOG"
check "step 27: fluconazole gives CYP2C9 a note in the fifth column" \
  grep -q $'^CYP2C9\t\\*1/\\*1\tNormal Metabolizer\tnormal\tnot computed: fluconazole inhibits CYP2C9' "$PHENO"

rm -f "$HITS"
run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
OUT=$(cat "$REC" 2>/dev/null)
check "control, step 27 without the hits table: the genes are listed as not checked" \
  has 'Not checked: no ClinVar hits table' "$OUT"
check "control, step 27 without the hits table: no valproic acid line" lacks 'WARNING: valproic acid' "$OUT"

# --- Nextflow: CLINVAR_SCREEN's hits reach CPIC_LOOKUP -----------------------------
if ! command -v nextflow >/dev/null; then
  fail "nextflow is not on PATH"
  finish
fi
CV="${G}/pgx-contra-clinvar"
rm -rf "${CV:?}"
mkdir -p "$CV"
read -r PC PP PR PA < <(awk -F'\t' 'NR == 2 {print $1, $2, $3, $4}' "${FIXTURE_DIR}/planted.tsv")
check "the fixture's planted record is readable (${PC:-?}:${PP:-?})" test -n "${PA:-}"
{
  printf '##fileformat=VCFv4.2\n'
  awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${G}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
  printf '##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="synthetic">\n'
  printf '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="synthetic">\n'
  printf '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="synthetic">\n'
  printf '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="synthetic">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf '%s\t%s\t9900001\t%s\t%s\t.\t.\tALLELEID=9900001;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;GENEINFO=POLG:5428\n' \
    "$PC" "$PP" "$PR" "$PA"
} > "${CV}/clinvar_polg.vcf"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -Oz -o pgx-contra-clinvar/clinvar_polg.vcf.gz pgx-contra-clinvar/clinvar_polg.vcf
  bcftools index -f -t pgx-contra-clinvar/clinvar_polg.vcf.gz"

SHEET="${CASE_TMP}/samplesheet.csv"
printf 'sample,sex,vcf,vcf_index\n%s,male,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET"
# nf_run <pipeline dir> <outdir> <log> [-resume]
nf_run() {
  ( cd "$CASE_TMP" && nextflow run "$1/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-contra" --input "$SHEET" \
      --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools pharmcat,clinvar,cpic \
      --clinvar "${CV}/clinvar_polg.vcf.gz" --clinvar_index "${CV}/clinvar_polg.vcf.gz.tbi" \
      --outdir "$2" --max_cpus 4 --max_memory 14.GB "${@:4}" ) > "$3" 2>&1
  local rc=$?
  sed 's/^/    | /' "$3" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
  return "$rc"
}
NF_OUT="${G}/nf-contra"
rm -rf "${NF_OUT:?}"
echo "+ nextflow run main.nf --tools pharmcat,clinvar,cpic (ClinVar: the planted record as POLG)"
nf_run "$REPO" "$NF_OUT" "${CASE_TMP}/nf.log"
check_eq "nextflow run exits 0" "$?" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
check_ge "CLINVAR_SCREEN: the POLG record is a hit" \
  "$(grep -c 'POLG:5428' "${NF_OUT}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv" 2>/dev/null || true)" 1
NF_REC=$(cat "${NF_OUT}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt" 2>/dev/null)
check "CPIC_LOOKUP: the POLG hit gives the valproic acid WARNING" \
  has "WARNING: valproic acid -- POLG ${PC}:${PP} ${PR}>${PA}, genotype [^,]+, ClinVar Pathogenic\." "$NF_REC"

# Control: the same run from a copy of the pipeline whose main.nf passes []
# for the contraindication table.
CTRL="${CASE_TMP}/pipeline-control"
rm -rf "${CTRL:?}"
mkdir -p "$CTRL"
for f in main.nf nextflow.config nextflow_schema.json versions.env assets bin conf modules workflows scripts; do
  cp -R "${REPO}/${f}" "${CTRL}/"
done
# shellcheck disable=SC2016  # ${projectDir} is Nextflow's, matched as written
LINE='Channel.value(file("${projectDir}/assets/pgx_contraindications.tsv", checkIfExists: true))'
check_eq "control: main.nf has the ch_contraindications line once" "$(grep -cF "$LINE" "${CTRL}/main.nf")" 1
python3 - "${CTRL}/main.nf" "$LINE" <<'PY'
import sys
p, line = sys.argv[1], sys.argv[2]
s = open(p).read()
open(p, "w").write(s.replace(line, "Channel.value([])"))
PY
check_eq "control: the copy's main.nf no longer reads the table" "$(grep -cF "$LINE" "${CTRL}/main.nf")" 0
CTRL_OUT="${G}/nf-contra-control"
rm -rf "${CTRL_OUT:?}"
echo "+ nextflow run <copy>/main.nf with ch_contraindications = [] (control)"
nf_run "$CTRL" "$CTRL_OUT" "${CASE_TMP}/nf-control.log" -resume
check_eq "control: nextflow run exits 0" "$?" 0
CTRL_REC=$(cat "${CTRL_OUT}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt" 2>/dev/null)
check "control: CPIC_LOOKUP wrote its report" test -n "$CTRL_REC"
check "control: without the table there is no valproic acid line" lacks 'WARNING: valproic acid' "$CTRL_REC"

rm -rf "${CTRL:?}" "${CTRL_OUT:?}"
finish
