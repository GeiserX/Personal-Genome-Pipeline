#!/usr/bin/env bash
# Step 27 reads the PharmCAT report into gene rows; a sentinel line means it parsed nothing.
. "$(dirname "$0")/lib.sh"

run_step 27-cpic-lookup.sh "$SAMPLE"
check_step_exit 27-cpic-lookup.sh

PHENO="${GENOME_DIR}/${SAMPLE}/cpic/${SAMPLE}_phenotypes.tsv"
FIRST=$(awk 'NR == 1 {print $1}' "$PHENO" 2>/dev/null)
check "phenotypes TSV does not start with a parse-failure sentinel (${FIRST:-empty})" \
  lacks '^(NO_JSON_FOUND|UNKNOWN_FORMAT|PARSE_EMPTY|)$' "${FIRST:-}"
check_ge "gene rows in the phenotypes TSV" "$(awk -F'\t' 'NF >= 3' "$PHENO" 2>/dev/null | wc -l | tr -d ' ')" 1
check "recommendations file exists" nonempty "${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt"

finish
