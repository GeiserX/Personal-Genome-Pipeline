#!/usr/bin/env bash
# Step 35's interval table (<id>_cn_intervals.tsv) and step 18's CNVpytor
# calls reach the couple screen and the Nextflow report.
#   1. step 39 on a made-up couple, FCN_A (male) and FCN_B (female), each
#      with an empty ClinVar screen, an interval table and a CNVpytor VCF in
#      their folders: the couple report lists FCN_A's del(GJB6-D13S1830) and
#      AZFc deletion as present, FCN_B's CFTRdele2,3 as absent with the
#      CNVpytor call beside it, and FCN_B's AZF line as not applicable.
#      Red control: with FCN_A's interval table gone the report has no line
#      of FCN_A's alleles and names the missing table.
#   2. COUPLE_CARRIER, run for real in its image through a two-line pipeline
#      on the same files: the same couple table and report as step 39.
#      Red control: a copy of the module without the two input kinds writes
#      no allele line.
#   3. HTML_REPORT, run for real through a two-line pipeline with the SMN
#      table and the interval table: the task places the interval table and
#      summary.json carries its rows with a change. Red control: a copy of
#      the module without that place line stops on the input.
#   4. main.nf passes both files on: a stub run of the stub pedigree (two
#      partners) with clinvar, parascopy, cnvpytor and html_report stages
#      each partner's interval table and CNVpytor VCF in COUPLE_CARRIER and
#      the interval table in each HTML_REPORT task. Red control: a copy of
#      the pipeline without the two lines stages none.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
check "an earlier case left ${VCF#"$GENOME_DIR"/}" nonempty "${VCF#"$GENOME_DIR"/}"

C="${E2E_WORK}/genome-followups-cn"
docker run --rm -v "${E2E_WORK}:/w" "$PYTHON_IMAGE" rm -rf /w/genome-followups-cn
A=FCN_A; B=FCN_B
for s in "$A" "$B"; do mkdir -p "${C}/${s}/clinvar" "${C}/${s}/paralogs" "${C}/${s}/cnvpytor"; done
G0="$GENOME_DIR"
export GENOME_DIR="$C"

ivh() { printf 'locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\testimate\treason\n'; }
iv() {  # iv LOCUS CHROM START END STATUS COPIES EXPECTED CHANGE REASON
  printf '%s\t%s name\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%s\n' "$1" "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9"
}
for s in "$A" "$B"; do
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$s" \
    > "${C}/${s}/clinvar/${s}_clinvar_hits.vcf"
done
{
  ivh
  iv GJB6_D13S1830 chr13 20223037 20531806 called 1 2 deletion .
  iv CFTR_dele2_3 chr7 117498312 117519392 called 2 2 "no change" .
  iv HBB_exon3_del chr11 5225254 5225875 called 2 2 "no change" .
  iv AZFa chrY 12537649 12920478 called 1 1 "no change" .
  iv AZFb chrY 19702297 20593154 called 1 1 "no change" .
  iv AZFc_DAZ chrY 23129354 24907040 called 2 4 deletion .
} > "${C}/${A}/paralogs/${A}_cn_intervals.tsv"
{
  ivh
  iv GJB6_D13S1830 chr13 20223037 20531806 called 2 2 "no change" .
  iv CFTR_dele2_3 chr7 117498312 117519392 called 2 2 "no change" .
  iv HBB_exon3_del chr11 5225254 5225875 called 2 2 "no change" .
  iv AZFa chrY 12537649 12920478 "not applicable" . . . no_y
} > "${C}/${B}/paralogs/${B}_cn_intervals.tsv"
# CNVpytor: FCN_A has no call, FCN_B a deletion over CFTRdele2,3 (absent by depth)
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -n > "${C}/${A}/cnvpytor/${A}_cnvs.vcf.gz"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr7\t117498001\tCNV1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=117520000\n' \
  | gzip -n > "${C}/${B}/cnvpytor/${B}_cnvs.vcf.gz"
printf 'sample\tsex\tfather\tmother\tpartner\n%s\tmale\t.\t.\t%s\n%s\tfemale\t.\t.\t%s\n' "$A" "$B" "$B" "$A" > "${C}/pedigree.tsv"

GJB6_A="- ${A}, del\(GJB6-D13S1830\): \*\*present\*\* \(read depth, step 35: 1 of 2 copies\); CNVpytor \(step 18\) calls no deletion here"
AZFC_A="- ${A}, AZFc deletion \(DAZ copies; gr/gr and b2/b3 leave two of four\): \*\*present\*\* \(read depth, step 35: 2 of 4 copies\)"
CFTR_B="- ${B}, CFTRdele2,3: \*\*absent\*\* \(read depth, step 35: 2 of 2 copies\); CNVpytor \(step 18\) calls a deletion here"
AZF_B="- ${B}, AZF regions: not applicable \(no Y chromosome\)"

# --- 1. step 39 -------------------------------------------------------------------
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
check "step 39 names both inputs of each partner" \
  sh -c "grep -q '${A}: paralogs/${A}_cn_intervals.tsv (step 35)' '$STEP_LOG' && grep -q '${B}: cnvpytor/${B}_cnvs.vcf.gz (step 18)' '$STEP_LOG'"
MD="${C}/couple/couple_${A}_${B}.md"
TXT=$(cat "$MD" 2>/dev/null)
grep -E "^- (${A}|${B}), " <<< "$TXT"
check "step 39: ${A}'s del(GJB6-D13S1830) present, CNVpytor disagreeing" has "$GJB6_A" "$TXT"
check "step 39: ${A}'s AZFc deletion present" has "$AZFC_A" "$TXT"
check "step 39: ${B}'s CFTRdele2,3 absent, the CNVpytor deletion beside it" has "$CFTR_B" "$TXT"
check "step 39: ${B}'s AZF line not applicable" has "$AZF_B" "$TXT"
cp "$MD" "${CASE_TMP}/step.md"
cp "${C}/couple/couple_${A}_${B}.tsv" "${CASE_TMP}/step.tsv"

mv "${C}/${A}/paralogs/${A}_cn_intervals.tsv" "${CASE_TMP}/a_iv.tsv"
run_step 39-couple.sh couple "${C}/pedigree.tsv"
check_step_exit 39-couple.sh
TXT0=$(cat "$MD" 2>/dev/null)
check "red: without ${A}'s interval table no line of ${A}'s alleles" lacks "^- ${A}, (del\(GJB6|AZFc)" "$TXT0"
check "red: and the report names the missing table" has "Carrier deletion alleles of ${A}: no interval table \(step 35\)" "$TXT0"
mv "${CASE_TMP}/a_iv.tsv" "${C}/${A}/paralogs/${A}_cn_intervals.tsv"

# --- 2. COUPLE_CARRIER --------------------------------------------------------------
NP="${CASE_TMP}/project"
rm -rf "$NP"; mkdir -p "$NP"
cp -R "${REPO}/bin" "${NP}/bin"
mkdir -p "${NP}/couple-nokinds"
sed -e "/'cn-intervals': '_cn_intervals.tsv', cnvpytor: '_cnvs.vcf.gz'/d" \
  -e "s|smn: '_smn_copy_number.tsv',\$|smn: '_smn_copy_number.tsv']|" "${REPO}/modules/local/couple/main.nf" \
  > "${NP}/couple-nokinds/main.nf"
check_eq "control copy of the couple module: the two kinds are gone" \
  "$(grep -c "_cn_intervals.tsv'" "${NP}/couple-nokinds/main.nf")" 0
nfconf() {  # nfconf OUT
  printf "params.outdir = '%s'\nparams.publish_dir_mode = 'copy'\nparams.family_name = 'couple'\ndocker.enabled = true\ndocker.runOptions = '--user %s:%s --network none'\n" \
    "$1" "$(id -u)" "$(id -g)" > "${NP}/nextflow.config"
}
nfrun() {  # nfrun NAME: the pipeline in NP/main.nf, log CASE_TMP/NAME.log
  echo "+ nextflow run ($1)"
  ( cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false -c "${REPO}/conf/containers.config" \
      -work-dir "${CASE_TMP}/work-$1" ) > "${CASE_TMP}/$1.log" 2>&1
  NF_RC=$?
  [ "$NF_RC" = 0 ] || grep -E 'ERROR|Error executing process|Command error' -A 8 "${CASE_TMP}/$1.log" | head -40
}
couple_run() {  # couple_run NAME MODULE
  rm -rf "${CASE_TMP:?}/$1"
  nfconf "${CASE_TMP}/$1"
  cat > "${NP}/main.nf" <<NF
include { COUPLE_CARRIER } from '$2'

workflow {
    COUPLE_CARRIER(file('${C}/pedigree.tsv'), ['${A}', '${B}'], ['male', 'female'],
        [file('${C}/${A}/clinvar/${A}_clinvar_hits.vcf'), file('${C}/${B}/clinvar/${B}_clinvar_hits.vcf'),
         file('${C}/${A}/paralogs/${A}_cn_intervals.tsv'), file('${C}/${B}/paralogs/${B}_cn_intervals.tsv'),
         file('${C}/${A}/cnvpytor/${A}_cnvs.vcf.gz'), file('${C}/${B}/cnvpytor/${B}_cnvs.vcf.gz')],
        file('${REPO}/assets/carrier_genes.tsv'), [], [], file('${REPO}/assets/par_grch38.bed'))
}
NF
  nfrun "$1"
  check_eq "COUPLE_CARRIER ($1) exits 0" "$NF_RC" 0
}
couple_run couple "${REPO}/modules/local/couple/main"
NMD="${CASE_TMP}/couple/couple/couple_${A}_${B}.md"
check "COUPLE_CARRIER writes step 39's report, byte for byte" cmp -s "$NMD" "${CASE_TMP}/step.md"
check "COUPLE_CARRIER writes step 39's table, byte for byte" \
  cmp -s "${CASE_TMP}/couple/couple/couple_${A}_${B}.tsv" "${CASE_TMP}/step.tsv"
check "COUPLE_CARRIER: ${A}'s del(GJB6-D13S1830) present" has "$GJB6_A" "$(cat "$NMD" 2>/dev/null)"
couple_run couple-nokinds "${NP}/couple-nokinds/main"
check "red: the module without the two kinds writes no allele line" \
  lacks "^- ${A}, del\(GJB6" "$(cat "${CASE_TMP}/couple-nokinds/couple/couple_${A}_${B}.md" 2>/dev/null)"

# --- 3. HTML_REPORT -----------------------------------------------------------------
T="${CASE_TMP}/tables"
rm -rf "$T"; mkdir -p "$T"
printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n' \
  > "${T}/${SAMPLE}_smn_copy_number.tsv"
cp "${C}/${A}/paralogs/${A}_cn_intervals.tsv" "${T}/${SAMPLE}_cn_intervals.tsv"
mkdir -p "${NP}/report-noplace"
sed '/_cn_intervals.tsv)  *place paralogs/d' "${REPO}/modules/local/html_report/main.nf" > "${NP}/report-noplace/main.nf"
check_eq "control copy of the report module: the place line is gone" \
  "$(grep -c '_cn_intervals.tsv)' "${NP}/report-noplace/main.nf")" 0
report_run() {  # report_run NAME MODULE
  rm -rf "${CASE_TMP:?}/$1"
  nfconf "${CASE_TMP}/$1"
  cat > "${NP}/main.nf" <<NF
include { HTML_REPORT } from '$2'

workflow {
    HTML_REPORT(channel.of(
        [[id: '${SAMPLE}', sex: 'male'], file('${G0}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz'),
         [file('${T}/${SAMPLE}_smn_copy_number.tsv'), file('${T}/${SAMPLE}_cn_intervals.tsv')]]))
}
NF
  nfrun "$1"
}
report_run report "${REPO}/modules/local/html_report/main"
check_eq "HTML_REPORT with the interval table exits 0" "$NF_RC" 0
changes=$(python3 - "${CASE_TMP}/report/${SAMPLE}/summary.json" <<'PY'
import json, sys
try:
    d = json.load(open(sys.argv[1]))
except (OSError, ValueError):
    print("no summary"); sys.exit(0)
ci = ((d["sections"].get("parascopy") or {}).get("data") or {}).get("cn_intervals") or {}
print(",".join(c["locus"] for c in ci.get("changes", [])) + "/" + str(ci.get("rows")))
PY
)
check_eq "summary.json: the interval rows with a change, and the row count" "$changes" "GJB6_D13S1830,AZFc_DAZ/6"
report_run report-noplace "${NP}/report-noplace/main"
check "red: the module without the place line stops (exit ${NF_RC})" test "$NF_RC" -ne 0
check "and says it has no place for the table" grep -q "no place for: ${SAMPLE}_cn_intervals.tsv" "${CASE_TMP}/report-noplace.log"

# --- 4. main.nf: stub run of the stub pedigree ---------------------------------------
staged_in() {  # staged_in WORK PROCESS PATTERN: tasks of PROCESS in WORK that stage a file matching PATTERN (ERE)
  find "$1" -name .command.run -exec grep -l "^### name: .*$2" {} + 2>/dev/null \
    | while read -r f; do grep -lE -- "$3" "$f"; done | wc -l | tr -d ' '
}
for n in stub stub-nowire; do
  rm -rf "${CASE_TMP:?}/${n}"; mkdir -p "${CASE_TMP}/${n}"
  git -C "$REPO" ls-files -z | (cd "$REPO" && tar --null -T - -cf -) | tar -xf - -C "${CASE_TMP}/${n}"
done
W1='BAM_ANALYSIS.out.cn_intervals.map'; W2='BAM_ANALYSIS.out.cn_intervals, SV.out.cnvpytor_vcf'
check_eq "main.nf joins the interval table into the report inputs once" "$(grep -c "$W1" "${CASE_TMP}/stub/main.nf")" 1
check_eq "main.nf hands both files to the couple screen once" "$(grep -c "$W2" "${CASE_TMP}/stub/main.nf")" 1
sed -i -e "/${W1}/d" -e "s/BAM_ANALYSIS.out.smn_copy_number,\$/BAM_ANALYSIS.out.smn_copy_number/" -e "/${W2}/d" \
  "${CASE_TMP}/stub-nowire/main.nf"
check_eq "control copy: both lines are gone" "$(grep -c -e "$W1" -e "$W2" "${CASE_TMP}/stub-nowire/main.nf")" 0
echo '{"kir": false}' > "${CASE_TMP}/no-kir.json"
for n in stub stub-nowire; do
  echo "+ nextflow run main.nf -profile test_all -stub --pedigree assets/stub/pedigree.tsv (${n})"
  ( cd "${CASE_TMP}/${n}" && nextflow run main.nf -profile test_all -stub -ansi-log false --outdir "${CASE_TMP}/${n}/out" \
      -work-dir "${CASE_TMP}/${n}/work" --pedigree assets/stub/pedigree.tsv --tools clinvar,parascopy,cnvpytor,html_report \
      -params-file "${CASE_TMP}/no-kir.json" ) > "${CASE_TMP}/${n}.log" 2>&1
  rc=$?
  tail -n 6 "${CASE_TMP}/${n}.log"
  check_eq "stub run (${n}) exits 0" "$rc" 0
done
S="${CASE_TMP}/stub/work"; S0="${CASE_TMP}/stub-nowire/work"
check_eq "COUPLE_CARRIER stages the interval table" "$(staged_in "$S" COUPLE_CARRIER '_cn_intervals.tsv')" 1
check_eq "COUPLE_CARRIER stages the CNVpytor VCF" "$(staged_in "$S" COUPLE_CARRIER '_cnvs.vcf.gz')" 1
check_eq "both partners' tables reach it" \
  "$(find "$S" -name .command.run -exec grep -l '^### name: .*COUPLE_CARRIER' {} + 2>/dev/null | while read -r f; do grep -ho 'stub_[a-z]*_cn_intervals.tsv' "$f"; done | sort -u | wc -l | tr -d ' ')" 2
check_eq "every HTML_REPORT task stages the interval table" \
  "$(staged_in "$S" HTML_REPORT '_cn_intervals.tsv')" "$(staged_in "$S" HTML_REPORT '_smn_copy_number.tsv')"
check_ge "HTML_REPORT tasks with the SMN table" "$(staged_in "$S" HTML_REPORT '_smn_copy_number.tsv')" 1
check_ge "control: the stub run without the lines still ran COUPLE_CARRIER" "$(staged_in "$S0" COUPLE_CARRIER '_clinvar_hits.vcf')" 1
check_eq "red: without the lines COUPLE_CARRIER stages neither file" \
  "$(staged_in "$S0" COUPLE_CARRIER '_cn_intervals[.]tsv|_cnvs[.]vcf[.]gz')" 0
check_eq "red: and no HTML_REPORT task stages the interval table" "$(staged_in "$S0" HTML_REPORT '_cn_intervals.tsv')" 0

rm -rf "${CASE_TMP:?}/stub" "${CASE_TMP:?}/stub-nowire" "${CASE_TMP:?}"/work-*
docker run --rm -v "${E2E_WORK}:/w" "$PYTHON_IMAGE" rm -rf /w/genome-followups-cn
finish
