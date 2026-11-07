#!/usr/bin/env bash
# Step 35's tables that no Parascopy model makes reach both reports and the
# couple screen, with the real images of steps 24, generate-report.sh and
# 39, on made-up tables in a GENOME_DIR of their own:
#   1. steps 24 and generate-report.sh on a sample with an SMN table, a
#      CYP21A2 table (a chimeric deletion: 1 copy, module total 3 of 4), the
#      gene-dosage table (HP 3, LCE3BC 1, CCL3L1 not assessed) and an interval
#      table with one 22q11.2 deletion: both reports show each of them, and
#      neither says CYP21A2 is not assessed. Red control: with the CYP21A2
#      table gone the text report says it is not in the output, with no
#      finding.
#   2. step 39 on a made-up couple, CNR_A (female) and CNR_B (male), each
#      with an empty ClinVar screen, an interval table with CFTRdele2,3
#      present and a CYP21A2 table (CNR_A the deletion, CNR_B none): step 39
#      passes both CYP21A2 tables, the couple table has CNR_A's CYP21A2
#      deletion, and the shared-risk heading names CFTRdele2,3 from read depth.
#      Red control: with CNR_A's CYP21A2 table gone there is no CYP21A2 row
#      and the report names the missing table.
. "$(dirname "$0")/lib.sh"

G="${CASE_TMP}/genome"
T="${SAMPLE}cnr"
rm -rf "$G"
mkdir -p "${G}/${T}/paralogs"
printf 'run\tscript\te2e\n' > "${G}/${T}/run_manifest.tsv"
P="${G}/${T}/paralogs"
printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n%s\n' \
  "chr5	70950920	70953320	SMN1	PASS	4	60	PASS	2,2	33,30	chr5:70075496-70077900:+" > "${P}/${T}_smn_copy_number.tsv"
cyp_del() {
  printf 'metric\tvalue\nstatus\tcalled\ncyp21a2_copies\t1\nmodule_total\t3\nexpected_total\t4\nestimate\t3.02\n'
  printf 'finding\tdeletion\nconversion_sites\t.\ntarget_depth\t45.30\ncontrol_depth\t30.10\nreason\t.\n'
  printf 'reading\tthe pattern of the 30 kb CYP21A1P/CYP21A2 chimeric deletion\n'
}
cyp_none() {
  printf 'metric\tvalue\nstatus\tcalled\ncyp21a2_copies\t2\nmodule_total\t4\nexpected_total\t4\nestimate\t4.01\n'
  printf 'finding\tnone\nconversion_sites\t.\ntarget_depth\t60.20\ncontrol_depth\t30.10\nreason\t.\nreading\t.\n'
}
cyp_del > "${P}/${T}_cyp21a2_psv.tsv"
{
  printf 'locus\tstatus\tcopies\testimate\treason\treading\n'
  printf 'HP\tcalled\t3\t3.04\t.\tthree copies of the duplicated HP exons: Hp2-1\n'
  printf 'LCE3BC\tcalled\t1\t0.98\t.\tLCE3B and LCE3C deleted on one chromosome\n'
  printf 'CCL3L1\tnot assessed\t.\t.\tcheck_failed\tthe depth check failed; see the log of step 35\n'
} > "${P}/${T}_cn_dosage.tsv"
ivh() { printf 'locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\testimate\treason\n'; }
{
  ivh
  printf '22q11.2_LCR_A-LCR_D\t22q11.2 recurrent region (A-D)\tchr22\t18924717\t21111383\tcalled\t1\t2\tdeletion\t1.02\t.\n'
  printf 'CFTR_dele2_3\tCFTRdele2,3\tchr7\t117498312\t117519392\tcalled\t2\t2\tno change\t1.97\t.\n'
} > "${P}/${T}_cn_intervals.tsv"

# --- 1. both reports ----------------------------------------------------------------
GENOME_DIR="$G" run_step 24-html-report.sh "$T"
check_step_exit 24-html-report.sh
HTML=$(cat "${G}/${T}/${T}_report.html" 2>/dev/null)
check "HTML: the CYP21A2 deletion" \
  has '>CYP21A2 \(RCCX module, from depth and six PSV sites\)<.*1 CYP21A2 copies, module total 3 of 4, finding deletion' "$HTML"
check "HTML: the HP count" has '>HP \(haptoglobin\) copies \(from depth\)<.*3 \(estimate 3.04\)' "$HTML"
check "HTML: the 22q11.2 deletion row" \
  has '<td>22q11.2_LCR_A-LCR_D</td><td>chr22:18924718-21111383</td><td>1 of 2</td><td>deletion</td>' "$HTML"
check "HTML: no claim that CYP21A2 is not assessed" lacks 'CYP21A2 has no Parascopy model and is not assessed' "$HTML"

GENOME_DIR="$G" run_step generate-report.sh "$T"
check_step_exit generate-report.sh
TXT=$(cat "${G}/${T}/${T}_report.txt" 2>/dev/null)
sed -n '/SMN1\/SMN2 Copy Number (Parascopy/,/^$/p' <<< "$TXT"
check "text: the CYP21A2 deletion" \
  has '^  CYP21A2 \(RCCX module, from depth and six PSV sites\): 1 CYP21A2 copies, module total 3 of 4, finding deletion' "$TXT"
check "text: the three dosage loci" \
  bash -c "grep -q '^  HP (haptoglobin) copies (from depth): 3' <<< \"\$1\" && grep -q '^  LCE3B/LCE3C copies (from depth): 1' <<< \"\$1\" \
    && grep -q '^  CCL3L1/CCL3L3 copies (from depth): not assessed (check_failed)' <<< \"\$1\"" _ "$TXT"
check "text: the interval counts and the deletion" \
  bash -c "grep -q 'from depth): 2 rows (2 called); 1 with a deletion or a duplication' <<< \"\$1\" \
    && grep -q '^    Deletion: 22q11.2 recurrent region (A-D) (22q11.2_LCR_A-LCR_D, chr22:18924718-21111383), 1 of 2 copies' <<< \"\$1\"" _ "$TXT"

mv "${P}/${T}_cyp21a2_psv.tsv" "${CASE_TMP}/cyp.tsv"
GENOME_DIR="$G" run_step generate-report.sh "$T"
TXT0=$(cat "${G}/${T}/${T}_report.txt" 2>/dev/null)
check "red: without the CYP21A2 table the text report says it is not in the output" \
  has "CYP21A2 \(RCCX module, from depth and six PSV sites\): not assessed \(step 35's <id>_cyp21a2_psv.tsv is not in this run's output\)" "$TXT0"
check "red: and gives no CYP21A2 finding" lacks 'finding deletion' "$TXT0"

# --- 2. step 39 -------------------------------------------------------------------
C="${G}/couple-case"
A=CNR_A; B=CNR_B
for s in "$A" "$B"; do
  mkdir -p "${C}/${s}/clinvar" "${C}/${s}/paralogs"
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$s" \
    > "${C}/${s}/clinvar/${s}_clinvar_hits.vcf"
  {
    ivh
    printf 'CFTR_dele2_3\tCFTRdele2,3\tchr7\t117498312\t117519392\tcalled\t1\t2\tdeletion\t1.01\t.\n'
  } > "${C}/${s}/paralogs/${s}_cn_intervals.tsv"
done
cyp_del > "${C}/${A}/paralogs/${A}_cyp21a2_psv.tsv"
cyp_none > "${C}/${B}/paralogs/${B}_cyp21a2_psv.tsv"
printf 'sample\tsex\tfather\tmother\tpartner\n%s\tfemale\t.\t.\t%s\n%s\tmale\t.\t.\t%s\n' "$A" "$B" "$B" "$A" > "${C}/pedigree.tsv"

GENOME_DIR="$C" run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
check "step 39 passes both CYP21A2 tables" \
  bash -c "grep -q '${A}: paralogs/${A}_cyp21a2_psv.tsv (step 35)' '$STEP_LOG' && grep -q '${B}: paralogs/${B}_cyp21a2_psv.tsv (step 35)' '$STEP_LOG'"
TSV=$(cat "${C}/couple/couple_${A}_${B}.tsv" 2>/dev/null)
MD=$(cat "${C}/couple/couple_${A}_${B}.md" 2>/dev/null)
check "step 39: CNR_A's CYP21A2 deletion is a row of the couple table" \
  has "^CYP21A2	.*one_partner	.*CYP21A2 deletion from depth and PSV sites \(CYP21A2 1 copies, RCCX module total 3 of 4\)" "$TSV"
check "step 39: the shared-risk heading names read depth" \
  has "No gene of the list where both partners carry a change in the gene table; read depth \(step 35\) adds a possible shared risk, below\." "$MD"
check "step 39: CFTRdele2,3 present in both is named" \
  has "^- \*\*CFTR\*\*: ${A}: CFTRdele2,3 present \(read depth, step 35: 1 of 2 copies\); ${B}: CFTRdele2,3 present" "$MD"

mv "${C}/${A}/paralogs/${A}_cyp21a2_psv.tsv" "${CASE_TMP}/a_cyp.tsv"
GENOME_DIR="$C" run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
check "red: without CNR_A's CYP21A2 table no CYP21A2 row" lacks "^CYP21A2	" "$(cat "${C}/couple/couple_${A}_${B}.tsv" 2>/dev/null)"
check "red: and the report names the missing table" \
  has "CYP21A2 deletion and conversions of ${A}: step 35's CYP21A2 table was not passed to the couple screen \(--cyp21a2\)\." \
  "$(cat "${C}/couple/couple_${A}_${B}.md" 2>/dev/null)"
docker run --rm -v "${CASE_TMP}:/w" "$PYTHON_IMAGE" rm -rf /w/genome

finish
