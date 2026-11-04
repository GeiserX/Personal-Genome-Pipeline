#!/usr/bin/env bash
# The 1.0 output contract (assets/output_contract.tsv) on the HG002 base outputs: every stable file
# the bash steps and the Nextflow run (case 60) wrote is there under its contracted name, with its
# required columns, VCF fields and JSON keys. Red control: a copy of the Nextflow results with one
# column dropped and one file renamed must fail, naming both.
. "$(dirname "$0")/lib.sh"

CHECK="${REPO}/scripts/ci/check-output-contract.py"
NF_OUT="${GENOME_DIR}/nf-results"

contract() {  # contract <dir> <leg> <log>: the run check, its output in <log>; returns its exit code
  python3 "$CHECK" --outputs "$1" --leg "$2" --sample "$SAMPLE" --observe > "$3" 2>&1
}

contract "$GENOME_DIR" bash "${CASE_TMP}/bash.log"; RC_BASH=$?
grep -v '^OBSERVED' "${CASE_TMP}/bash.log"
check_eq "the bash steps' outputs keep the contract" "$RC_BASH" 0
check_ge "stable bash files checked" "$(sed -n 's/^outputs (bash): \([0-9]*\) stable.*/\1/p' "${CASE_TMP}/bash.log")" 20

contract "$NF_OUT" nextflow "${CASE_TMP}/nf.log"; RC_NF=$?
grep -v '^OBSERVED' "${CASE_TMP}/nf.log"
check_eq "the Nextflow outputs keep the contract" "$RC_NF" 0
check_ge "stable Nextflow files checked" "$(sed -n 's/^outputs (nextflow): \([0-9]*\) stable.*/\1/p' "${CASE_TMP}/nf.log")" 10

# What each leg wrote, with its fields: the lines to read when the contract is written again.
echo "--- observed files (path, format, fields)"
grep -h '^OBSERVED' "${CASE_TMP}/bash.log" "${CASE_TMP}/nf.log"

# Red control on a copy: the ClinVar hits table loses its clnsig column, the ROH summary is renamed.
PLANT="${CASE_TMP}/planted"
rm -rf "$PLANT"; mkdir -p "$PLANT"
cp -a "${NF_OUT}/${SAMPLE}" "${PLANT}/"
HITS="${PLANT}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv"
ROH="${PLANT}/${SAMPLE}/roh/${SAMPLE}_roh_summary.txt"
if [ -f "$HITS" ] && [ -f "$ROH" ]; then
  awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "clnsig") c = i}
    {out = ""; for (i = 1; i <= NF; i++) if (i != c) out = out (out == "" ? "" : OFS) $i; print out}' \
    "$HITS" > "${HITS}.new" && mv "${HITS}.new" "$HITS"
  mv "$ROH" "${ROH%_roh_summary.txt}_roh_segments.txt"
fi
contract "$PLANT" nextflow "${CASE_TMP}/planted.log"; RC_PLANT=$?
grep '^FAIL' "${CASE_TMP}/planted.log"
check_eq "the planted copy fails the contract" "$RC_PLANT" 1
check "the dropped clnsig column is named" grep -q "BREAKING: ${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv .* has no clnsig" "${CASE_TMP}/planted.log"
check "the renamed ROH summary is named" grep -q "BREAKING: the stable file ${SAMPLE}/roh/${SAMPLE}_roh_summary.txt .* is missing" "${CASE_TMP}/planted.log"

finish
