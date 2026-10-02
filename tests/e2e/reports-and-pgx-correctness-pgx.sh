#!/usr/bin/env bash
# Step 27 with bin/pgx_parse.py, on the real PharmCAT report of case 31 and on
# two synthetic reports:
#   - HG002: every gene with a non-normal phenotype lists its drugs, and the
#     PharmCAT/pypgx comparison (moved here from step 32) is written;
#   - a report where PharmCAT has no CYP2D6 result while pypgx called it: the
#     CPIC report prints the warning with the drugs CYP2D6 affects, and a gene
#     the report names no drug for falls back to the static table;
#   - a report that parses to zero genes: the step exits non-zero and says so.
# The HG002 report.json is copied to the job's logs (the e2e-logs artifact):
# it is the committed parser fixture tests/fixtures/pharmcat/report-3.2.0.json.
. "$(dirname "$0")/lib.sh"

JSON="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.report.json"
check "case 31 left a PharmCAT report.json" test -s "$JSON"
cp "$JSON" "${E2E_WORK}/logs/${SAMPLE}.pharmcat-report.json" 2>/dev/null || true

# --- HG002 (step 27 ran in case 39, after pypgx in case 38) -------------------
PHENO="${GENOME_DIR}/${SAMPLE}/cpic/${SAMPLE}_phenotypes.tsv"
REC="${GENOME_DIR}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt"
check "the phenotypes table has the Status column" has '^Gene	Diplotype	Phenotype	Status$' "$(head -1 "$PHENO" 2>/dev/null)"
cat "$REC" 2>/dev/null
check "the recommendations do not say PARSING FAILED" lacks 'PARSING FAILED' "$(cat "$REC" 2>/dev/null)"
mapfile -t NONNORMAL < <(awk -F'\t' 'NR > 1 && $4 == "non-normal" {print $1}' "$PHENO" 2>/dev/null)
echo "genes with a non-normal phenotype: ${NONNORMAL[*]:-none}"
for g in "${NONNORMAL[@]}"; do
  # The gene's block runs from "  GENE -- " to its "Action:" line.
  BLOCK=$(awk -v g="$g" 'index($0, "  " g " -- ") == 1 {on = 1} on {print} on && /Action:/ {exit}' "$REC")
  check "${g}: its medications block lists drugs" has '^    (Drugs with guidance|Drugs PharmCAT links|Drugs \(pipeline fallback)' "$BLOCK"
  check "${g}: not left without a drug list" lacks 'is not in the drug table' "$BLOCK"
done
COMP="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pharmcat_comparison.tsv"
check "step 27 wrote the PharmCAT/pypgx comparison" has '^Gene	PharmCAT_diplotype	pypgx_diplotype	Match	Called_by$' "$(head -1 "$COMP" 2>/dev/null)"
check_ge "comparison rows" "$(awk 'NR > 1' "$COMP" 2>/dev/null | grep -c . || true)" 1
check "the comparison has a CYP2D6 row (pypgx called it in case 38)" grep -q '^CYP2D6	' "$COMP"

# --- synthetic: no CYP2D6 result from PharmCAT, a pypgx call -----------------
S2="${SAMPLE}pgx"
mkdir -p "${GENOME_DIR}/${S2}/pharmcat" "${GENOME_DIR}/${S2}/pypgx"
cat > "${GENOME_DIR}/${S2}/pharmcat/${S2}.report.json" <<'JSON'
{"pharmcatVersion": "3.2.0",
 "genes": {
  "CYP2D6": {"sourceDiplotypes": [{"allele1": {"name": "Unknown"}, "allele2": {"name": "Unknown"},
                                   "label": "Unknown/Unknown", "phenotypes": ["No Result"]}]},
  "CYP2C19": {"sourceDiplotypes": [{"allele1": {"name": "*2"}, "allele2": {"name": "*2"},
                                    "label": "*2/*2", "phenotypes": ["Poor Metabolizer"]}]},
  "SLCO1B1": {"sourceDiplotypes": [{"allele1": {"name": "*1"}, "allele2": {"name": "*1"},
                                    "label": "*1/*1", "phenotypes": ["Normal Function"]}]}
 },
 "drugs": {}}
JSON
printf 'Gene\tDiplotype\tPhenotype\tCNV_call\tSource\nCYP2D6\t*1/*4\tIntermediate Metabolizer\t.\tbam\nCYP2C19\t*2/*2\tPoor Metabolizer\t.\tvcf\n' \
  > "${GENOME_DIR}/${S2}/pypgx/${S2}_pypgx_summary.tsv"
run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
REC2="${GENOME_DIR}/${S2}/cpic/${S2}_cpic_recommendations.txt"
OUT2=$(cat "$REC2" 2>/dev/null)
check "the pypgx warning names CYP2D6 and its call" \
  has 'WARNING: PharmCAT has no result for CYP2D6, but pypgx \(step 32\) called \*1/\*4' "$OUT2"
check "the warning lists the drugs CYP2D6 affects" has 'Drugs affected by CYP2D6: .*codeine' "$OUT2"
check "CYP2C19 (no drug in the report) falls back to the static table" \
  has 'Drugs \(pipeline fallback table; the report names none\): clopidogrel' "$OUT2"
check "CYP2D6 is listed as not callable" has 'CYP2D6 -- No Result \(not callable' "$OUT2"
check "the comparison marks CYP2D6 as pypgx only" \
  grep -q $'^CYP2D6\tUnknown/Unknown\t\\*1/\\*4\tpypgx only' "${GENOME_DIR}/${S2}/pypgx/${S2}_pharmcat_comparison.tsv"

# --- synthetic: a report that yields no gene ---------------------------------
S3="${SAMPLE}pgxempty"
mkdir -p "${GENOME_DIR}/${S3}/pharmcat"
echo '{"pharmcatVersion": "3.2.0", "genes": {"CYP2C19": {"geneSymbol": "CYP2C19"}}}' > "${GENOME_DIR}/${S3}/pharmcat/${S3}.report.json"
run_step 27-cpic-lookup.sh "$S3"
check "a report with no readable gene makes step 27 exit non-zero" test "$STEP_RC" -ne 0
check "and the recommendations say PARSING FAILED" grep -q 'PARSING FAILED' "${GENOME_DIR}/${S3}/cpic/${S3}_cpic_recommendations.txt"

finish
