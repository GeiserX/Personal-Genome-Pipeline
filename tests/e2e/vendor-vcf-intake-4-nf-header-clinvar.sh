#!/usr/bin/env bash
# Two inputs that stopped a Nextflow run inside a task:
#   - a VCF whose header has a ##FILTER line with escaped quotes (valid VCF,
#     what `bcftools filter -s LowDP -e '... GT!="0/0"'` writes): PharmCAT's
#     Java step stopped on it. PHARMCAT now rewrites its own copy and calls
#     exactly what case 2 called without that line;
#   - the full chr-renamed ClinVar file as --clinvar, with a record on
#     NT_113889.1 as the real file has: bcftools norm stopped with exit 255 on
#     the contig the reference lacks. CLINVAR_SCREEN now leaves it out.
# The run selects html_report without roh, mito_haplogroup and cpic, so those
# three cards must say "Not run".
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

# The chr slice with the escaped ##FILTER line (case 1 wrote the header line).
mkdir -p "${INTAKE}/lowdp"
bcf annotate --no-version -h intake/lowdp.hdr -Oz -o "intake/lowdp/${SAMPLE}.vcf.gz" "${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
bcf index -f -t "intake/lowdp/${SAMPLE}.vcf.gz"
check_eq "the input header carries the escaped quotes" \
  "$(bcf view -h "intake/lowdp/${SAMPLE}.vcf.gz" | grep -c 'GT!=\\"0/0\\"' || true)" 1

# The fixture's full chr-renamed ClinVar plus one record on NT_113889.1, a
# contig the reference does not have.
mkdir -p "${INTAKE}/clinvar"
{
  bcf view -h clinvar/clinvar_chr.vcf.gz | grep '^##'
  echo '##contig=<ID=NT_113889.1,length=161147>'
  bcf view -h clinvar/clinvar_chr.vcf.gz | grep '^#CHROM'
  bcf view -H clinvar/clinvar_chr.vcf.gz
  printf 'NT_113889.1\t100\t9999999\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=NONE:0;CLNREVSTAT=no_assertion_criteria_provided\n'
} > "${INTAKE}/clinvar/unsorted.vcf"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools sort -Oz -o intake/clinvar/clinvar_chr_nt.vcf.gz intake/clinvar/unsorted.vcf
  bcftools index -f -t intake/clinvar/clinvar_chr_nt.vcf.gz"
check_eq "the ClinVar copy has its NT_113889.1 record" \
  "$(bcf view -H intake/clinvar/clinvar_chr_nt.vcf.gz NT_113889.1 2>/dev/null | wc -l | tr -d ' ')" 1

sheet "${CASE_TMP}/lowdp.csv" "${INTAKE}/lowdp/${SAMPLE}.vcf.gz"
nf_run lowdp "${CASE_TMP}/lowdp.csv" pharmcat,clinvar,html_report \
  --clinvar "${INTAKE}/clinvar/clinvar_chr_nt.vcf.gz" \
  --clinvar_index "${INTAKE}/clinvar/clinvar_chr_nt.vcf.gz.tbi"
check_eq "nextflow run exits 0" "$NF_RC" 0
check_eq "PHARMCAT completed" "$(trace_col ':PHARMCAT ' status)" COMPLETED
check_eq "CLINVAR_SCREEN completed" "$(trace_col ':CLINVAR_SCREEN ' status)" COMPLETED
R="${NF_OUT}/${SAMPLE}"

for f in match phenotype; do
  A=$(json_calls "${INTAKE}/chr-run/${SAMPLE}.${f}.json")
  B=$(json_calls "${R}/pharmcat/${SAMPLE}.${f}.json")
  check "PHARMCAT ${f}.json is readable" lacks '^unreadable' "$B"
  if [ "$A" = "$B" ]; then
    pass "PHARMCAT ${f}.json calls are the same as without the escaped line"
  else
    fail "PHARMCAT ${f}.json calls differ from case 2's"
    diff <(tr ',' '\n' <<< "$A") <(tr ',' '\n' <<< "$B") | head -20
  fi
done

check "CLINVAR_SCREEN wrote a readable hits VCF" vcf_ok "nf-lowdp/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf"
check "CLINVAR_SCREEN wrote the hits TSV" test -s "${R}/clinvar/${SAMPLE}_clinvar_hits.tsv"
check_eq "no hit on NT_113889.1" "$(grep -c '^NT_113889' "${R}/clinvar/${SAMPLE}_clinvar_hits.vcf" 2>/dev/null || true)" 0

HTML="${R}/${SAMPLE}_report.html"
check "the report exists" test -s "$HTML"
for h in 'Runs of Homozygosity' 'Mitochondrial Haplogroup' 'CPIC Drug Recommendations'; do
  check_eq "report card '${h}' says Not run" "$(html_stat "$HTML" "$h" Status)" 'Not run'
done
check_eq "report card 'Pharmacogenomics' is filled" "$(html_stat "$HTML" Pharmacogenomics 'PharmCAT report')" Complete

finish
