#!/usr/bin/env bash
# shellcheck disable=SC2016  # the $ fields are awk's, inside its single-quoted programs
# Step 48 (scripts/48-exomiser.sh) and the EXOMISER_* modules, with the real
# Exomiser image, on Exomiser's own worked example: the Pfeiffer syndrome exome
# its repository ships (a planted FGFR2 variant, 10:123256215 T>G) with the
# four HPO terms its documentation gives for it, where FGFR2 ranks first.
# The data is Exomiser's cut-down test release (hg19, 2.8 GB to download, 4.7
# GB unpacked), which covers the variants of that exome; the hg38 release the
# pipeline installs is about 37 GB, too large for this runner.
#   - the single-sample run: FGFR2 at rank 1 with that variant, a dominant
#     mode, an ACMG class, the term names from the release's hp.obo, the note;
#   - red controls: no terms, a malformed term, an unknown term and a missing
#     data folder stop the step before Exomiser runs; unrelated terms lower
#     FGFR2's symptom match; a variants file with a column cut stops the report;
#   - the family run (EXOMISER_FAMILY) on the Pfeiffer quartet laid out as step
#     38 leaves it, parents and sibling unaffected: FGFR2 (de novo) at rank 1;
#     the sibling, who lacks the variant, named affected: FGFR2 drops out of
#     the dominant ranks;
#   - the three modules in a small Nextflow workflow: the same files as step 48,
#     byte for byte;
#   - main.nf -stub with --hpo: the tool turns on, its files reach HTML_REPORT;
#     --tools exomiser without --hpo, and a malformed --hpo, stop the run.
. "$(dirname "$0")/lib.sh"

# Exomiser's test release (exomiser.readthedocs.io, "Installation"; size and
# sha256 measured 2026-10-25, Exomiser publishes no checksum).
XURL=https://g-879a9f.f5dc97.75bc.dn.glob.us/data/exomiser-test-data.zip
XSHA=869e358b9a3e73207d133e1a378d4d0cb6448a04eb4131ec9f75f3e382b0bd01
# The Pfeiffer files of the Exomiser 15.1.0 tag.
PURL=https://raw.githubusercontent.com/exomiser/Exomiser/15.1.0/exomiser-cli/src/main/resources/vcf
PF_SHA=ac497e6d941b87bb96f0b4f2b0b612ec7e4d1d8d8bd2ba68d111961518f18f71
QUARTET_SHA=db5fa80abb740672405bd57a5fc3daf1748f373e3f2589310fd3a5d9053b17b7
PED_SHA=4e357bcdab7e79dd9e643a2a0cba5b739030fa13d39b21340c0a54cd1ce0ae64
# Exomiser's docs: Pfeiffer syndrome as HPO terms.
HPO_PF=HP:0001156,HP:0001363,HP:0011304,HP:0010055
# Hepatomegaly, Proteinuria: nothing FGFR2's diseases are known for.
HPO_OTHER=HP:0002240,HP:0000093

DREL=reference/exomiser-test
D="${GENOME_DIR}/${DREL}"
export EXOMISER_DATA=$DREL EXOMISER_ASSEMBLY=hg19 EXOMISER_MEM=6

# get URL SHA FILE: download, resuming after a dropped connection (the data
# server closes long ones), until the file has the pinned sha256.
get() {
  local i
  for i in 1 2 3 4 5 6 7 8 9 10; do
    if [ -s "$3" ] && echo "$2  $3" | sha256sum -c --status; then return 0; fi
    [ "$i" = 1 ] || echo "download of ${1##*/} broke, resuming (try ${i})"
    curl -fsSL --retry 3 --retry-delay 10 -C - -o "$3" "$1" || true
  done
  echo "$2  $3" | sha256sum -c --status
}

if [ ! -d "${D}/test_hg19" ] || [ ! -d "${D}/test_phenotype" ]; then
  echo "+ downloading Exomiser's test release"
  Z="${E2E_WORK}/exomiser-test-data.zip"
  if ! get "$XURL" "$XSHA" "$Z"; then
    fail "Exomiser's test release did not download with sha256 ${XSHA}"
    finish
  fi
  rm -rf "${D}" "${D}.part" && mkdir -p "${D}.part"
  unzip -q "$Z" -d "${D}.part" && mv "${D}.part" "$D"
  rm -f "$Z"
fi
check "the test release has test_hg19 and test_phenotype" test -d "${D}/test_hg19" -a -d "${D}/test_phenotype"
check "and the HPO ontology" test -s "${D}/test_phenotype/hp.obo"

F="${CASE_TMP}/pfeiffer-files"
mkdir -p "$F"
if ! { get "${PURL}/Pfeiffer.vcf.gz" "$PF_SHA" "${F}/Pfeiffer.vcf.gz" \
        && get "${PURL}/Pfeiffer-quartet.vcf.gz" "$QUARTET_SHA" "${F}/Pfeiffer-quartet.vcf.gz" \
        && get "${PURL}/Pfeiffer-quartet.ped" "$PED_SHA" "${F}/Pfeiffer-quartet.ped"; }; then
  fail "the Pfeiffer files did not download with their pinned sha256"
  finish
fi

# One person: <GENOME_DIR>/pfeiffer/vcf/pfeiffer.vcf.gz, as step 03 leaves it.
S=pfeiffer
rm -rf "${GENOME_DIR:?}/${S}"
mkdir -p "${GENOME_DIR}/${S}/vcf"
cp "${F}/Pfeiffer.vcf.gz" "${GENOME_DIR}/${S}/vcf/${S}.vcf.gz"

X="${GENOME_DIR}/${S}/exomiser"
# cell FILE ROW COLUMN: a value of the top table (ROW 1 is the first data row).
cell() {
  awk -F'\t' -v r="$2" -v c="$3" 'NR == 1 {for (i = 1; i <= NF; i++) k[$i] = i; next} NR == r + 1 {print $k[c]; exit}' "$1" 2>/dev/null
}
# gene_rows FILE GENE [MOI]: rows of the top table for GENE (under MOI).
gene_rows() {
  awk -F'\t' -v g="$2" -v m="${3:-}" 'NR == 1 {for (i = 1; i <= NF; i++) k[$i] = i; next}
    $k["gene"] == g && (m == "" || $k["moi"] == m) {n++} END {print n + 0}' "$1" 2>/dev/null
}
# pheno_score FILE GENE: FGFR2's symptom match in the top table, or in Exomiser's genes file.
pheno_of() {
  awk -F'\t' -v g="$2" 'NR == 1 {for (i = 1; i <= NF; i++) {h = $i; sub(/^#/, "", h); k[h] = i}; next}
    $k["GENE_SYMBOL"] == g {print $k["EXOMISER_GENE_PHENO_SCORE"]; exit}' "$1" 2>/dev/null
}

# --- 1. the single-sample run ----------------------------------------------------------
HPO=$HPO_PF run_step 48-exomiser.sh "$S"
check_step_exit 48-exomiser.sh
TOP="${X}/${S}_exomiser_top.tsv"
check_eq "rank 1 is FGFR2" "$(cell "$TOP" 1 gene)" FGFR2
check_eq "with the planted variant" "$(cell "$TOP" 1 variant)" "10:123256215 T>G"
check_eq "heterozygous" "$(cell "$TOP" 1 genotype)" "0/1"
check_eq "under a dominant mode" "$(cell "$TOP" 1 moi)" AD
check "an ACMG class is given" has '^(PATHOGENIC|LIKELY_PATHOGENIC|UNCERTAIN_SIGNIFICANCE|LIKELY_BENIGN|BENIGN)$' "$(cell "$TOP" 1 acmg)"
check "an Exomiser score above 0.9" awk -v s="$(cell "$TOP" 1 exomiser_score)" 'BEGIN {exit !(s + 0 > 0.9)}'
check "at most 10 ranks" awk -F'\t' 'NR > 1 {r[$1]} END {exit !(length(r) >= 1 && length(r) <= 10)}' "$TOP"
TXT=$(cat "${X}/${S}_exomiser.txt" 2>/dev/null)
check "the words name the terms with their hp.obo labels" has '^Symptoms searched: Brachydactyly \(HP:0001156\), .*Broad thumb \(HP:0011304\)' "$(head -n 1 <<< "$TXT")"
check "and rank 1 in a sentence" has '^1\. FGFR2 \(autosomal dominant\): Exomiser score' "$TXT"
check "and end with the not-a-diagnosis note" has '^A ranking to bring to a clinical geneticist, not a diagnosis\.' "$(tail -n 1 <<< "$TXT")"
check_eq "the terms file holds the four terms" "$(awk 'NR > 1' "${X}/${S}_exomiser_terms.tsv" 2>/dev/null | wc -l | tr -d ' ')" 4
check "Exomiser's genes file and HTML are there" test -s "${X}/${S}_exomiser.genes.tsv" -a -s "${X}/${S}_exomiser.html"
PHENO_PF=$(pheno_of "${X}/${S}_exomiser.genes.tsv" FGFR2)
cp -R "$X" "${CASE_TMP}/bash-out"

# --- 2. red controls ----------------------------------------------------------------------
# Each red run starts without the variants file of the run before, so
# step_stops can tell whether Exomiser ran.
VARS="${X}/${S}_exomiser.variants.tsv"
# step_stops DESC MESSAGE: the step exited non-zero, said MESSAGE, and Exomiser did not run.
step_stops() {
  check "$1: the step stops (exit ${STEP_RC})" test "$STEP_RC" -ne 0
  check "$1: it says why" grep -Eq -- "$2" "$STEP_LOG"
  check "$1: Exomiser did not run" test ! -e "$VARS"
}
rm -f "$VARS"; HPO="" run_step 48-exomiser.sh "$S"
step_stops "no HPO terms" 'no HPO terms'
rm -f "$VARS"; HPO=HP:12345 run_step 48-exomiser.sh "$S"
step_stops "a malformed term" 'HP:12345'
rm -f "$VARS"; HPO=HP:0001156,HP:9999999 run_step 48-exomiser.sh "$S"
step_stops "a term the ontology does not have" 'HP:9999999'
rm -f "$VARS"; HPO=$HPO_PF EXOMISER_DATA=reference/exomiser-none run_step 48-exomiser.sh "$S"
step_stops "no data folder" 'setup.sh --exomiser-data'

HPO=$HPO_OTHER run_step 48-exomiser.sh "$S"
check_step_exit 48-exomiser.sh
PHENO_OTHER=$(pheno_of "${X}/${S}_exomiser.genes.tsv" FGFR2)
check "unrelated terms give FGFR2 a lower symptom match (${PHENO_OTHER:-none} < ${PHENO_PF:-none})" \
  awk -v a="${PHENO_OTHER:-0}" -v b="${PHENO_PF:-0}" 'BEGIN {exit !(b + 0 > 0 && a + 0 < b + 0)}'

# A variants file with its MOI column cut: the report stops and names it.
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "MOI") m = i}
  {o = ""; for (i = 1; i <= NF; i++) if (i != m) o = o (o == "" ? "" : "\t") $i; print o}' \
  "${CASE_TMP}/bash-out/${S}_exomiser.variants.tsv" > "${CASE_TMP}/cut.variants.tsv"
check "the cut file still has rows" test "$(wc -l < "${CASE_TMP}/cut.variants.tsv")" -gt 1
OUT=$(python3 -I "${REPO}/bin/exomiser_report.py" report --sample "$S" --variants "${CASE_TMP}/cut.variants.tsv" \
  --terms "${CASE_TMP}/bash-out/${S}_exomiser_terms.tsv" --out "${CASE_TMP}/cut.tsv" --text "${CASE_TMP}/cut.txt" 2>&1)
RC=$?
check "a variants file without MOI: the report stops (exit ${RC})" test "$RC" -ne 0
check "and names the column" has 'lacks the column MOI' "$OUT"

# --- 3. the family run --------------------------------------------------------------------
# Step 38's folder: family.vcf.gz (bgzip, indexed) and pedigree.tsv.
FAM=pfeiffer-family
P=ISDBM322017
rm -rf "${GENOME_DIR:?}/${FAM}" "${GENOME_DIR:?}/${P}"
mkdir -p "${GENOME_DIR}/${FAM}"
cp "${F}/Pfeiffer-quartet.vcf.gz" "${GENOME_DIR}/${FAM}/raw.vcf.gz"
bcf view -Oz -o "${FAM}/family.vcf.gz" "${FAM}/raw.vcf.gz" && bcf index -t "${FAM}/family.vcf.gz"
awk -v OFS='\t' 'BEGIN {print "sample", "sex", "father", "mother", "partner"}
  !/^#/ && NF >= 6 {print $2, ($5 == 1 ? "male" : $5 == 2 ? "female" : "unknown"), ($3 == "0" ? "" : $3), ($4 == "0" ? "" : $4), ""}' \
  "${F}/Pfeiffer-quartet.ped" > "${GENOME_DIR}/${FAM}/pedigree.tsv"
check_eq "the pedigree names the four members" "$(awk 'NR > 1' "${GENOME_DIR}/${FAM}/pedigree.tsv" | wc -l | tr -d ' ')" 4
XF="${GENOME_DIR}/${P}/exomiser"

HPO=$HPO_PF EXOMISER_FAMILY=$FAM EXOMISER_UNAFFECTED=ISDBM322015,ISDBM322016,ISDBM322018 run_step 48-exomiser.sh "$P"
check_step_exit 48-exomiser.sh
check_eq "family, parents and sibling unaffected: rank 1 is FGFR2" "$(cell "${XF}/${P}_exomiser_top.tsv" 1 gene)" FGFR2
check_eq "with the de novo variant" "$(cell "${XF}/${P}_exomiser_top.tsv" 1 variant)" "10:123256215 T>G"
check "the PED gives the proband status 2 and the others 1" \
  awk -F'\t' -v p="$P" '{s[$2] = $6} END {exit !(s[p] == 2 && s["ISDBM322015"] == 1 && s["ISDBM322016"] == 1 && s["ISDBM322018"] == 1)}' "${XF}/${P}_exomiser.ped"
check "and the parents of the proband" awk -F'\t' -v p="$P" '$2 == p && $3 == "ISDBM322016" && $4 == "ISDBM322018" {f = 1} END {exit !f}' "${XF}/${P}_exomiser.ped"

HPO=$HPO_PF EXOMISER_FAMILY=$FAM EXOMISER_AFFECTED=ISDBM322015 EXOMISER_UNAFFECTED=ISDBM322016,ISDBM322018 run_step 48-exomiser.sh "$P"
check_step_exit 48-exomiser.sh
check_eq "the sibling without the variant named affected: FGFR2 leaves the dominant ranks" \
  "$(gene_rows "${XF}/${P}_exomiser_top.tsv" FGFR2 AD)" 0
check "and rank 1 is another gene" test "$(cell "${XF}/${P}_exomiser_top.tsv" 1 gene)" != FGFR2

HPO=$HPO_PF EXOMISER_FAMILY=$FAM EXOMISER_UNAFFECTED=$P run_step 48-exomiser.sh "$P"
check "the proband named unaffected: the step stops (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "and says why" grep -q 'it cannot be unaffected' "$STEP_LOG"

# --- 4. the modules, in a small workflow --------------------------------------------------
if command -v nextflow >/dev/null; then
  NP="${CASE_TMP}/project"
  rm -rf "$NP" && mkdir -p "$NP"
  cp -R "${REPO}/bin" "${NP}/bin"
  cat > "${NP}/main.nf" <<NF
include { EXOMISER_PREPARE; EXOMISER; EXOMISER_REPORT } from '${REPO}/modules/local/exomiser/main'

workflow {
    data = channel.value(file('${D}'))
    EXOMISER_PREPARE(channel.of([[id: '${S}'], file('${GENOME_DIR}/${S}/vcf/${S}.vcf.gz')]), '${HPO_PF}', data)
    EXOMISER(EXOMISER_PREPARE.out.analysis, data)
    EXOMISER_REPORT(EXOMISER.out.variants.join(EXOMISER_PREPARE.out.terms))
}
NF
  cat > "${NP}/nextflow.config" <<NFCONF
params.outdir = '${CASE_TMP}/nf-out'
params.publish_dir_mode = 'copy'
params.exomiser_assembly = 'hg19'
process.resourceLimits = [ cpus: 4, memory: '8.GB' ]
docker.enabled = true
docker.runOptions = '--user $(id -u):$(id -g) --network none'
NFCONF
  echo "+ nextflow run (EXOMISER_PREPARE, EXOMISER, EXOMISER_REPORT on ${S})"
  ( cd "$CASE_TMP" && nextflow run "${NP}/main.nf" -ansi-log false \
      -c "${REPO}/conf/containers.config" -work-dir "${CASE_TMP}/work" ) 2>&1 | tee "$STEP_LOG"
  NF_RC=${PIPESTATUS[0]}
  check_eq "nextflow run exits 0" "$NF_RC" 0
  [ "$NF_RC" = 0 ] || grep -E 'Error executing process|Command exit status|Command error' -A 8 "${CASE_TMP}/.nextflow.log" | head -60
  for f in "${S}_exomiser_top.tsv" "${S}_exomiser.txt" "${S}_exomiser_terms.tsv" "${S}_exomiser.genes.tsv"; do
    check "the module's ${f} is step 48's, byte for byte" cmp -s "${CASE_TMP}/bash-out/${f}" "${CASE_TMP}/nf-out/${S}/exomiser/${f}"
  done
  check "the module publishes Exomiser's HTML" test -s "${CASE_TMP}/nf-out/${S}/exomiser/${S}_exomiser.html"
else
  fail "nextflow is not on PATH"
fi

# --- 5. main.nf, stub ---------------------------------------------------------------------
ST="${REPO}/assets/stub"
# stub_run NAME ARGS...: main.nf -profile test_all -stub on a one-sample sheet; sets RC, LOG and DIR
# (--kir false: the profile turns KIR on, which needs hla_typing in --tools)
stub_run() {
  DIR="${CASE_TMP}/stub-$1"
  shift
  rm -rf "$DIR" && mkdir -p "$DIR"
  printf 'sample,bam,bam_index,vcf,vcf_index,sex\nS1,%s,%s,%s,%s,male\n' \
    "${ST}/sample.bam" "${ST}/sample.bam.bai" "${ST}/sample.vcf.gz" "${ST}/sample.vcf.gz.tbi" > "${DIR}/samplesheet.csv"
  echo "+ nextflow run main.nf -profile test_all -stub $*"
  ( cd "$DIR" && nextflow run "${REPO}/main.nf" -profile test_all -stub -ansi-log false \
      -work-dir "${DIR}/work" --input "${DIR}/samplesheet.csv" --outdir "${DIR}/out" --kir false "$@" ) > "${DIR}/run.log" 2>&1
  RC=$?
  LOG=$(cat "${DIR}/run.log")
  echo "$LOG" | tail -n 30
}
if command -v nextflow >/dev/null; then
  stub_run on --tools html_report --hpo "$HPO_PF" --exomiser_data "${ST}/db"
  check_eq "--hpo with --tools html_report: the stub run exits 0" "$RC" 0
  check "the EXOMISER task ran" has 'EXOMISER \(S1\)' "$LOG"
  check "the top table is published" test -e "${DIR}/out/S1/exomiser/S1_exomiser_top.tsv"
  check "and the terms" test -e "${DIR}/out/S1/exomiser/S1_exomiser_terms.tsv"
  REP=$(grep -l "^### name: [\"']HTML_REPORT (" "${DIR}"/work/*/*/.command.run 2>/dev/null | head -n 1)
  check "HTML_REPORT received the top table" test -n "$REP" -a -e "$(dirname "${REP:-/x}")/S1_exomiser_top.tsv"

  stub_run off --tools exomiser,html_report --exomiser_data "${ST}/db"
  check "--tools exomiser without --hpo: the run stops (exit ${RC})" test "$RC" -ne 0
  check "and says to give the terms with --hpo" has 'give them with --hpo' "$LOG"
  stub_run bad --tools html_report --hpo HP:12 --exomiser_data "${ST}/db"
  check "a malformed --hpo: the run stops (exit ${RC})" test "$RC" -ne 0
  check "and names the term" has 'HP:12 is not an HPO id' "$LOG"
fi

# 4.7 GB of data and the work folders: no later case reads them.
rm -rf "${GENOME_DIR:?}/${FAM}" "${GENOME_DIR:?}/${P}" "${D:?}" "${CASE_TMP}/work" "${CASE_TMP}"/stub-*/work
finish
