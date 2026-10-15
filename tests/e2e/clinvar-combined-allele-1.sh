#!/usr/bin/env bash
# Step 06 and the CLINVAR_SCREEN module on a delins the caller wrote as two
# records 3 bp apart (tests/fixtures/vcf/combined_allele_*.vcf, at real chr22
# positions of the fixture reference). ClinVar holds a 2-star record for the
# combined allele and 1-star records for its two parts. The combined record
# must be a hit, both part hits must say "part of a combined allele", and a
# single SNV hit must be unchanged; both reports show the note, and the
# module writes the same hits as the bash step.
. "$(dirname "$0")/lib.sh"

FX="${REPO}/tests/fixtures/vcf"
S2=cca
G2="${E2E_WORK}/cca-genome"   # a GENOME_DIR of its own: its ClinVar file is the fixture's
rm -rf "${G2:?}"
mkdir -p "${G2}/reference" "${G2}/clinvar" "${G2}/${S2}/vcf"
for f in GRCh38_no_alt_analysis_set.fasta GRCh38_no_alt_analysis_set.fasta.fai GRCh38_no_alt_analysis_set.dict; do
  ln "${GENOME_DIR}/reference/${f}" "${G2}/reference/${f}" 2>/dev/null || cp "${GENOME_DIR}/reference/${f}" "${G2}/reference/${f}"
done
cp "${FX}/combined_allele_sample.vcf" "${G2}/${S2}/vcf/${S2}.vcf"
cp "${FX}/combined_allele_clinvar.vcf" "${G2}/clinvar/clinvar_pathogenic_chr.vcf"
export GENOME_DIR="$G2"   # lib.sh's bcf and sam, and the steps, read this one now
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools view -Oz -o ${S2}/vcf/${S2}.vcf.gz ${S2}/vcf/${S2}.vcf
  bcftools index -f -t ${S2}/vcf/${S2}.vcf.gz
  bcftools view -Oz -o clinvar/clinvar_pathogenic_chr.vcf.gz clinvar/clinvar_pathogenic_chr.vcf
  bcftools index -f -t clinvar/clinvar_pathogenic_chr.vcf.gz"
# The fixture records' REF bases are the reference's: the two parts and the
# bases between them, and the SNV.
check_eq "the reference under the two parts (chr22:30000037-30000042)" \
  "$(sam faidx reference/GRCh38_no_alt_analysis_set.fasta chr22:30000037-30000042 | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]')" AGCCAT
check_eq "the reference at the SNV (chr22:30000210)" \
  "$(sam faidx reference/GRCh38_no_alt_analysis_set.fasta chr22:30000210-30000210 | sed 1d | tr -d '\n' | tr '[:lower:]' '[:upper:]')" C

# --- bash step 06 -------------------------------------------------------------
run_step 06-clinvar-screen.sh "$S2"
check_step_exit 06-clinvar-screen.sh
HITS="${G2}/${S2}/clinvar/${S2}_clinvar_hits.vcf"
TSV="${G2}/${S2}/clinvar/${S2}_clinvar_hits.tsv"
# hit POS REF ALT: the hit's ID and INFO, tab-separated; empty when there is none.
hit() { awk -F'\t' -v p="$1" -v r="$2" -v a="$3" '!/^#/ && $2 == p && $4 == r && $5 == a {print $3 "\t" $8}' "$HITS" 2>/dev/null; }
COMB=$(hit 30000038 GCCAT CATGGG)
check "the combined allele is a hit: ClinVar's 2-star record 9900003" \
  has '^9900003[[:space:]].*CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts' "${COMB:-none}"
check "the combined hit names its two parts" \
  has 'PGP_JOINED_FROM=chr22:30000037:AGC:A\|chr22:30000042:T:TGGG' "${COMB:-none}"
for p in "30000037 AGC A 9900001" "30000042 T TGGG 9900002"; do
  read -r pos ref alt id <<< "$p"
  H=$(hit "$pos" "$ref" "$alt")
  check "the 1-star part ${pos} ${ref}>${alt} is a hit (${id}), marked part of the combined allele" \
    has "^${id}[[:space:]].*PGP_PART_OF=chr22:30000038:GCCAT:CATGGG" "${H:-none}"
done
SNV=$(hit 30000210 C T)
check "the single SNV hit is unchanged: ClinVar's record 9900004" has '^9900004[[:space:]]' "${SNV:-none}"
check "the single SNV hit carries no combined-allele tag" lacks 'PGP_' "${SNV:-none}"
check_eq "hits in the VCF" "$(grep -vc '^#' "$HITS" 2>/dev/null || true)" 4
check_eq "rows in the hits TSV" "$(awk 'NR > 1' "$TSV" 2>/dev/null | grep -c . || true)" 4
check "the step's log lists the parts as parts of the combined allele" \
  grep -q 'part of a combined allele: chr22:30000038 GCCAT>CATGGG' "$STEP_LOG"

# --- both reports ---------------------------------------------------------------
run_step 24-html-report.sh "$S2"
check_step_exit 24-html-report.sh
HTML="${G2}/${S2}/${S2}_report.html"
check "the HTML hits table has the combined-allele column" grep -q '<th>Stars</th><th>Combined allele</th>' "$HTML"
# The fixture's parts are unphased het calls, so every note says the parts may
# sit on different haplotypes.
UNPH='; unphased, so the parts may sit on different haplotypes'
check "the HTML report marks a part hit" grep -qF "<td>part of a combined allele: chr22:30000038 GCCAT&gt;CATGGG${UNPH}</td>" "$HTML"
check "the HTML report lists the combined hit with 2 stars" \
  grep -q "<td>30000038</td>.*<td>2</td><td>combined allele of chr22:30000037 AGC&gt;A + chr22:30000042 T&gt;TGGG${UNPH}</td>" "$HTML"
# The combined-allele cell follows Stars; step 06's cDNA read check column may follow it.
check "the HTML row of the single SNV has an empty combined-allele cell" grep -q '<td>30000210</td>.*<td>[0-4]</td><td></td>' "$HTML"
run_step generate-report.sh "$S2"
check_step_exit generate-report.sh
check "the text report marks a part hit" \
  grep -qF "(part of a combined allele: chr22:30000038 GCCAT>CATGGG${UNPH})" "${G2}/${S2}/${S2}_report.txt"

# --- the Nextflow module --------------------------------------------------------
if command -v nextflow >/dev/null; then
  SHEET="${CASE_TMP}/sheet.csv"
  printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$S2" "${G2}/${S2}/vcf/${S2}.vcf.gz" "${G2}/${S2}/vcf/${S2}.vcf.gz.tbi" > "$SHEET"
  echo "+ nextflow run main.nf --tools clinvar"
  (
    cd "$CASE_TMP" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
      -work-dir "${E2E_WORK}/nf-work-cca" --input "$SHEET" \
      --reference "${G2}/reference/GRCh38_no_alt_analysis_set.fasta" --tools clinvar \
      --clinvar "${G2}/clinvar/clinvar_pathogenic_chr.vcf.gz" \
      --clinvar_index "${G2}/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi" \
      --outdir "${G2}/nf" --max_cpus 4 --max_memory 14.GB
  ) > "${CASE_TMP}/nf.log" 2>&1
  NF_RC=$?
  sed 's/^/    | /' "${CASE_TMP}/nf.log" | grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists'
  check_eq "nextflow run exits 0" "$NF_RC" 0
  NF_HITS="${G2}/nf/${S2}/clinvar/${S2}_clinvar_hits.vcf"
  check_eq "CLINVAR_SCREEN writes the same hits as step 06 (CHROM to INFO, GT)" \
    "$(grep -v '^#' "$NF_HITS" 2>/dev/null | cut -f1-8,10 | md5sum)" "$(grep -v '^#' "$HITS" | cut -f1-8,10 | md5sum)"
  check_eq "CLINVAR_SCREEN hits" "$(grep -vc '^#' "$NF_HITS" 2>/dev/null || true)" 4
else
  fail "nextflow is not on PATH"
fi

finish
