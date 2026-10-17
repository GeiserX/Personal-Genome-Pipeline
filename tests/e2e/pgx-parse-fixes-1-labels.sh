#!/usr/bin/env bash
# The PGx labels on the real tools, after cases 31, 38, 39 and
# nextflow-hardening-2-nf:
#   - step 32 and the PYPGX module write 'no pypgx phenotype rule' for a gene
#     pypgx 0.26.0 has no phenotype rule for (CYP4F2, VKORC1, COMT, MTHFR in
#     its gene table), where they wrote 'Indeterminate' before; a gene with a
#     rule (CYP2C19) keeps pypgx's phenotype;
#   - step 27's comparison carries both tools' phenotypes;
#   - step 27 reads the genotype at CYP4F2's *3 position from the real
#     PharmCAT 3.4.0 report.json (the slice does not cover it: not callable),
#     and on a report with C/T there it writes the labelled single-position
#     fallback with the CPIC warfarin rule.
. "$(dirname "$0")/lib.sh"

NO_RULE='^(CYP4F2|VKORC1|COMT|MTHFR)$'
# check_labels NAME SUMMARY: every called no-rule gene says so, and the gene
# with a rule does not.
check_labels() {
  local name=$1 tsv=$2 n bad
  n=$(awk -F'\t' -v re="$NO_RULE" 'NR > 1 && $1 ~ re && $2 != "FAILED"' "$tsv" 2>/dev/null | grep -c . || true)
  bad=$(awk -F'\t' -v re="$NO_RULE" 'NR > 1 && $1 ~ re && $2 != "FAILED" && $3 != "no pypgx phenotype rule" {print $1 ": " $3}' \
    "$tsv" 2>/dev/null)
  check_ge "${name}: called genes pypgx has no phenotype rule for" "$n" 1
  check "${name}: each says 'no pypgx phenotype rule' (${bad:-all do})" test -z "$bad"
  check "${name}: CYP2C19, which has a rule, keeps pypgx's phenotype" \
    lacks 'no pypgx phenotype rule' "$(awk -F'\t' '$1 == "CYP2C19"' "$tsv" 2>/dev/null)"
}
check_labels "step 32" "${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv"
check_labels "PYPGX module" "${GENOME_DIR}/nf-hardening/${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv"

COMP="${GENOME_DIR}/${SAMPLE}/pypgx/${SAMPLE}_pharmcat_comparison.tsv"
check "the comparison has both phenotype columns" \
  has '^Gene	PharmCAT_diplotype	pypgx_diplotype	Match	Called_by	PharmCAT_phenotype	pypgx_phenotype$' \
  "$(head -1 "$COMP" 2>/dev/null)"
check "the comparison's CYP4F2 row carries pypgx's 'no pypgx phenotype rule'" \
  grep -q $'^CYP4F2\t.*\tno pypgx phenotype rule$' "$COMP"

REC="${GENOME_DIR}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt"
check "HG002: the real report's rs2108622 entry is read, with no genotype in the slice" \
  grep -q 'CYP4F2 -- PharmCAT has no result; rs2108622 has no genotype in the VCF PharmCAT read: not callable' "$REC"

# --- synthetic: CYP4F2 No Result, C/T at rs2108622 ----------------------------
S2="${SAMPLE}cyp4f2"
mkdir -p "${GENOME_DIR}/${S2}/pharmcat"
cat > "${GENOME_DIR}/${S2}/pharmcat/${S2}.report.json" <<'JSON'
{"pharmcatVersion": "3.4.0",
 "genes": {
  "CYP4F2": {"geneSymbol": "CYP4F2", "callSource": "MATCHER",
             "sourceDiplotypes": [{"allele1": {"name": "Unknown"}, "allele2": {"name": "Unknown"},
                                   "label": "Unknown/Unknown", "phenotypes": []}],
             "variants": [{"chromosome": "chr19", "position": 15879621, "dbSnpId": "rs2108622",
                           "call": "C/T", "referenceAllele": "C"}]}},
 "drugs": {}}
JSON
run_step 27-cpic-lookup.sh "$S2"
check_step_exit 27-cpic-lookup.sh
OUT2=$(cat "${GENOME_DIR}/${S2}/cpic/${S2}_cpic_recommendations.txt" 2>/dev/null)
check "CYP4F2 C/T: the labelled fallback with the genotype" \
  has 'WARNING: PharmCAT has no result for CYP4F2; single-position fallback, not a PharmCAT call: rs2108622 reads C/T' "$OUT2"
check "CYP4F2 C/T: the CPIC warfarin rule under it" has 'Warfarin \(CPIC\): .*5-10%' "$OUT2"
check "CYP4F2 stays not called in the phenotypes table" \
  grep -q $'^CYP4F2\tUnknown/Unknown\tNo Result\tnot called$' "${GENOME_DIR}/${S2}/cpic/${S2}_phenotypes.tsv"

finish
