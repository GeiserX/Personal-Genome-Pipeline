#!/usr/bin/env bash
# PharmCAT gets T1K's HLA types as outside calls, and the CPIC report has an
# HLA section sourced from the outside-call file.
#
# Bash: step 36 on HG002 (HLA from case bash-step-08, pypgx from case 38, no
# Cyrius), then steps 07 and 27. pypgx alone is one caller, so CYP2D6 stays
# indeterminate and does not reach PharmCAT.
# Nextflow: the leg of case nextflow-from-fastq-2 (pharmcat, cpic and
# hla_typing) ran PGX_CONSENSUS before PHARMCAT; its outputs are checked here.
. "$(dirname "$0")/lib.sh"

run_step 36-pgx-consensus.sh "$SAMPLE"
check_step_exit 36-pgx-consensus.sh
D="${GENOME_DIR}/${SAMPLE}/pgx_consensus"
CALLS="${D}/${SAMPLE}_outside_calls.tsv"
CONS="${D}/${SAMPLE}_pgx_consensus.tsv"
cat "$CALLS" "$CONS" 2>/dev/null
for g in HLA-A HLA-B; do
  check "${g} is an outside call (two-field alleles)" grep -Eq "^${g}"$'\t''\*[0-9]+:[0-9]+/\*[0-9]+:[0-9]+$' "$CALLS"
done
check "CYP2D6 is not passed on (pypgx alone)" lacks '^CYP2D6' "$(cat "$CALLS" 2>/dev/null)"
check "the table says CYP2D6 is indeterminate, from one caller" \
  grep -q $'^CYP2D6\tindeterminate\tno\tone caller only' "$CONS"

run_step 07-pharmacogenomics.sh "$SAMPLE"
check_step_exit 07-pharmacogenomics.sh
check "step 07 passed the outside calls" has '^Outside calls \(step 36\)' "$(cat "$STEP_LOG")"
# source REPORT GENE: PharmCAT's callSource for GENE
source_of() {
  python3 -c 'import json, sys; g = json.load(open(sys.argv[1]))["genes"][sys.argv[2]]; print(g.get("callSource", ""))' "$1" "$2" 2>/dev/null
}
JSON="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.report.json"
for g in HLA-A HLA-B; do
  check_eq "PharmCAT reports ${g} as an outside call" "$(source_of "$JSON" "$g")" OUTSIDE
done
check_eq "PharmCAT has no CYP2D6 outside call" "$(source_of "$JSON" CYP2D6)" NONE

run_step 27-cpic-lookup.sh "$SAMPLE"
check_step_exit 27-cpic-lookup.sh
REC="${GENOME_DIR}/${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt"
sed -n '/^Calls From Other Tools/,/^$/p' "$REC" 2>/dev/null
check "the CPIC report has the outside-call section" grep -q '^Calls From Other Tools (outside calls, step 36):' "$REC"
for g in HLA-A HLA-B; do
  check "the CPIC report lists ${g} as passed from T1K and seen by PharmCAT" \
    grep -Eq "^  ${g} +\*[0-9:/*]+ +passed to PharmCAT from T1K \(step 08\); PharmCAT reports it as an outside call\.$" "$REC"
done
check "the CPIC report says CYP2D6 was held back" grep -Eq '^  CYP2D6 +indeterminate +not passed to PharmCAT: one caller only' "$REC"

# --- the Nextflow leg ------------------------------------------------------------
P=HG002P
R="${GENOME_DIR}/nf-fastq/${P}"
check "PGX_CONSENSUS published the outside calls" grep -Eq $'^HLA-A\t' "${R}/pgx_consensus/${P}_outside_calls.tsv"
check_eq "PHARMCAT (Nextflow) reports HLA-B as an outside call" "$(source_of "${R}/pharmcat/${P}.report.json" HLA-B)" OUTSIDE
check "CPIC_LOOKUP lists HLA-B as passed from T1K" \
  grep -Eq '^  HLA-B +\*[0-9:/*]+ +passed to PharmCAT from T1K' "${R}/cpic/${P}_cpic_recommendations.txt"

finish
