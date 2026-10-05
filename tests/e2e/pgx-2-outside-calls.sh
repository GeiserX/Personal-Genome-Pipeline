#!/usr/bin/env bash
# PharmCAT gets T1K's HLA types as outside calls, and the CPIC report has an
# HLA section sourced from the outside-call file.
#
# Bash: steps 36, 07 and 27 on HG002pgx, a copy of HG002's VCF, gVCF, HLA
# types (case bash-step-08) and pypgx output (case 38), no Cyrius; a copy, so
# the later cases still read HG002's own PharmCAT report of case 31. pypgx
# alone is one caller, so CYP2D6 stays indeterminate and does not reach
# PharmCAT.
# Nextflow: the leg of case nextflow-from-fastq-2 (pharmcat, cpic and
# hla_typing) ran PGX_CONSENSUS before PHARMCAT; its outputs are checked here.
. "$(dirname "$0")/lib.sh"

SAMPLE_SRC=$SAMPLE
SAMPLE=${SAMPLE_SRC}pgx
rm -rf "${GENOME_DIR:?}/${SAMPLE}"
mkdir -p "${GENOME_DIR}/${SAMPLE}"
for d in vcf hla_t1k pypgx; do
  cp -r "${GENOME_DIR}/${SAMPLE_SRC}/${d}" "${GENOME_DIR}/${SAMPLE}/${d}"
  for f in "${GENOME_DIR}/${SAMPLE}/${d}/${SAMPLE_SRC}"[._]*; do
    [ -e "$f" ] && mv "$f" "${f%/*}/${SAMPLE}${f##*/"${SAMPLE_SRC}"}"
  done
done
rm -f "${GENOME_DIR}/${SAMPLE}/vcf/"*report* "${GENOME_DIR}/${SAMPLE}/vcf/"*.json
find "${GENOME_DIR}/${SAMPLE}" -maxdepth 2 -type f | head -40

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
# pypgx alone (the fixture's slice may give it no call): never passed on
check "the table says CYP2D6 is indeterminate: pypgx alone" \
  grep -Eq $'^CYP2D6\tindeterminate\tno\t(one caller only|no caller made a call)' "$CONS"
echo "- PGx consensus on the fixture: $(awk -F'\t' '$1 == "CYP2D6" {print $2 ", " $4 " (" $5 ")"}' "$CONS" 2>/dev/null)" >> "$E2E_NOTES"

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
check "the CPIC report says CYP2D6 was held back" grep -Eq '^  CYP2D6 +indeterminate +not passed to PharmCAT: (one caller only|no caller made a call)' "$REC"

# --- the Nextflow leg ------------------------------------------------------------
P=HG002P
R="${GENOME_DIR}/nf-fastq/${P}"
check "PGX_CONSENSUS published the outside calls" grep -Eq $'^HLA-A\t' "${R}/pgx_consensus/${P}_outside_calls.tsv"
check_eq "PHARMCAT (Nextflow) reports HLA-B as an outside call" "$(source_of "${R}/pharmcat/${P}.report.json" HLA-B)" OUTSIDE
check "CPIC_LOOKUP lists HLA-B as passed from T1K" \
  grep -Eq '^  HLA-B +\*[0-9:/*]+ +passed to PharmCAT from T1K' "${R}/cpic/${P}_cpic_recommendations.txt"

finish
