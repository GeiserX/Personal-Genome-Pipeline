#!/usr/bin/env bash
# The opt-in medication list of step 27 and CPIC_LOOKUP, in the python image.
# The list is synthetic (no one's medications) and sits outside GENOME_DIR, so
# the step's own read-only mount of it is what is tested:
#   - step 27 on a synthetic report (CYP2D6 *1/*2 normal metabolizer, activity
#     score 2.0; CYP2C19 *1/*1 normal metabolizer): paroxetine (strong CYP2D6
#     inhibitor) gives 'Poor Metabolizer (phenoconverted by paroxetine)' in a
#     fifth column beside the unchanged genotype phenotype, fluconazole a
#     CYP2C19 note, and a name not in the table is listed as not checked;
#   - the same report without MEDICATIONS keeps the four-column table;
#   - a MEDICATIONS file that does not exist stops the step;
#   - Nextflow with --tools pharmcat,cpic --medications on the HG002 VCF: the
#     CPIC_LOOKUP task reads the list (the section and the fifth column are
#     there; HG002's CYP2D6 has no result on the fixture, so it is not
#     converted, and the line says why).
. "$(dirname "$0")/lib.sh"

G="$GENOME_DIR"
MED="${CASE_TMP}/medications.txt"
printf 'Paroxetine 20 mg\nfluconazole\nnotadrug\n# a comment\n' > "$MED"

S2="${SAMPLE}meds"
mkdir -p "${G}/${S2}/pharmcat"
cat > "${G}/${S2}/pharmcat/${S2}.report.json" <<'JSON'
{"pharmcatVersion": "3.4.0",
 "genes": {
  "CYP2D6": {"geneSymbol": "CYP2D6", "callSource": "MATCHER",
             "sourceDiplotypes": [{"label": "*1/*2", "allele1": {"name": "*1"}, "allele2": {"name": "*2"},
                                   "phenotypes": ["Normal Metabolizer"], "activityScore": "2.0"}]},
  "CYP2C19": {"geneSymbol": "CYP2C19", "callSource": "MATCHER",
              "sourceDiplotypes": [{"label": "*1/*1", "allele1": {"name": "*1"}, "allele2": {"name": "*1"},
                                    "phenotypes": ["Normal Metabolizer"]}]}},
 "drugs": {}}
JSON
PHENO="${G}/${S2}/cpic/${S2}_phenotypes.tsv"
REC="${G}/${S2}/cpic/${S2}_cpic_recommendations.txt"

echo "+ MEDICATIONS=${MED}"
MEDICATIONS="$MED" run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
check "the phenotypes table has the fifth column" \
  grep -q $'^Gene\tDiplotype\tPhenotype\tStatus\tPhenoconverted_phenotype$' "$PHENO"
check "CYP2D6: genotype phenotype kept, phenoconverted by paroxetine beside it" \
  grep -q $'^CYP2D6\t\\*1/\\*2\tNormal Metabolizer\tnormal\tPoor Metabolizer (phenoconverted by paroxetine)$' "$PHENO"
check "CYP2C19: a note, no converted phenotype" \
  grep -q $'^CYP2C19\t.*\tnot computed: fluconazole inhibits CYP2C19' "$PHENO"
OUT=$(cat "$REC" 2>/dev/null)
check "the report's WARNING line names both phenotypes" \
  has 'WARNING: CYP2D6 genotype phenotype Normal Metabolizer; phenoconverted phenotype Poor Metabolizer \(phenoconverted by paroxetine\)' "$OUT"
check "a name not in the inhibitor table is listed" has 'Not in the inhibitor table: notadrug\.' "$OUT"
check "the section says it is an aid for the clinician" has 'an aid for your clinician, not advice' "$OUT"

run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
check "without MEDICATIONS the table keeps four columns" \
  grep -q $'^Gene\tDiplotype\tPhenotype\tStatus$' "$PHENO"
check "without MEDICATIONS there is no phenoconversion section" \
  lacks 'Medications and Phenoconversion' "$(cat "$REC" 2>/dev/null)"

MEDICATIONS="${CASE_TMP}/missing.txt" run_step 27-cpic-lookup.sh "$S2"
check "a MEDICATIONS file that does not exist stops the step (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "and says which file" has 'MEDICATIONS=.*missing\.txt is not a readable file' "$(cat "$STEP_LOG")"

# --- Nextflow: CPIC_LOOKUP with --medications -------------------------------------
if ! command -v nextflow >/dev/null; then
  fail "nextflow is not on PATH"
  finish
fi
NF_OUT="${G}/nf-medications"
WORK="${E2E_WORK}/nf-work-medications"
SHEET="${CASE_TMP}/samplesheet.csv"
rm -rf "$NF_OUT" "$WORK"
printf 'sample,sex,vcf,vcf_index\n%s,male,%s,%s\n' "$SAMPLE" \
  "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" "${G}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "$SHEET"
echo "+ nextflow run main.nf -profile docker --tools pharmcat,cpic --medications"
( cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false -work-dir "$WORK" \
    --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --tools pharmcat,cpic \
    --medications "$MED" --outdir "$NF_OUT" --max_cpus 4 --max_memory 14.GB ) 2>&1 | tee "$STEP_LOG"
NF_RC=${PIPESTATUS[0]}
check_eq "nextflow run exits 0" "$NF_RC" 0
cp "${CASE_TMP}/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
NF_PHENO="${NF_OUT}/${SAMPLE}/cpic/${SAMPLE}_phenotypes.tsv"
NF_REC=$(cat "${NF_OUT}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt" 2>/dev/null)
check "CPIC_LOOKUP: the phenotypes table has the fifth column" \
  grep -q $'^Gene\tDiplotype\tPhenotype\tStatus\tPhenoconverted_phenotype$' "$NF_PHENO"
check "CPIC_LOOKUP: the phenoconversion section reads the list" \
  has 'Medication list: medications\.txt \(3 lines\)\.' "$NF_REC"
check "CPIC_LOOKUP: HG002's CYP2D6 on paroxetine is not converted, with the reason" \
  has 'CYP2D6 -- on paroxetine \(strong CYP2D6 inhibitor\): not converted, the gene is' "$NF_REC"
check "CPIC_LOOKUP: a name not in the inhibitor table is listed" has 'Not in the inhibitor table: notadrug\.' "$NF_REC"

finish
