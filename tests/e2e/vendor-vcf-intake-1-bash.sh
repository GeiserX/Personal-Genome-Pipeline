#!/usr/bin/env bash
# Bash steps on vendor-style input:
#   - step 06 (run by case 30) writes clinvar/<sample>_clinvar_hits.tsv, one
#     row per hit with genotype, gene, significance and review status, and
#     only hits the person carries;
#   - step 07 runs on a VCF whose header has a ##FILTER line with escaped
#     quotes, as `bcftools filter -s LowDP -e '... GT!="0/0"'` writes it, and
#     PharmCAT calls exactly what it called on the same records without that
#     line (case 31).
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

# --- Step 06: the hits table ---------------------------------------------------
HITS="${G}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf"
TSV="${G}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.tsv"
check "step 06 wrote the hits TSV" test -s "$TSV"
check_eq "hits TSV header" "$(head -n 1 "$TSV" 2>/dev/null)" \
  "$(printf 'chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat')"
check_eq "hits TSV has one row per hit in the VCF" \
  "$(awk 'NR > 1' "$TSV" 2>/dev/null | grep -c . || true)" "$(grep -vc '^#' "$HITS" 2>/dev/null || true)"
check_ge "the planted hit is a TSV row with its gene ($(planted gene))" \
  "$(awk -F'\t' -v c="$(planted chrom)" -v p="$(planted pos)" -v g="$(planted gene)" \
      '$1 == c && $2 == p && index($7, g) && $5 ~ /1/' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 1
check_eq "every hit's genotype carries an ALT allele" \
  "$(awk -F'\t' 'NR > 1 && $5 !~ /[1-9]/' "$TSV" 2>/dev/null | wc -l | tr -d ' ')" 0

# --- Step 07: a backslash in a ## header line ----------------------------------
S2="${SAMPLE}_lowdp"
mkdir -p "${G}/${S2}/vcf"
printf '%s\n' '##FILTER=<ID=LowDP,Description="Set if true: FORMAT/DP<10 && GT!=\"0/0\"">' > "${INTAKE}/lowdp.hdr"
bcf annotate --no-version -h intake/lowdp.hdr -Oz -o "${S2}/vcf/${S2}.vcf.gz" "${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
bcf index -f -t "${S2}/vcf/${S2}.vcf.gz"
check_eq "the copy's header carries the escaped quotes" \
  "$(bcf view -h "${S2}/vcf/${S2}.vcf.gz" | grep -c 'GT!=\\"0/0\\"' || true)" 1
check_eq "the copy has the same records" \
  "$(bcf view -H "${S2}/vcf/${S2}.vcf.gz" | md5sum)" "$(bcf view -H "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | md5sum)"

run_step 07-pharmacogenomics.sh "$S2"
check_step_exit 07-pharmacogenomics.sh
for f in match phenotype; do
  A=$(json_calls "${G}/${SAMPLE}/vcf/${SAMPLE}.${f}.json")
  B=$(json_calls "${G}/${S2}/vcf/${S2}.${f}.json")
  check "step 07 ${f}.json is readable with the escaped line" lacks '^unreadable' "$B"
  if [ "$A" = "$B" ]; then
    pass "step 07 ${f}.json calls are the same as without the escaped line"
  else
    fail "step 07 ${f}.json calls differ from case 31's"
    diff <(tr ',' '\n' <<< "$A") <(tr ',' '\n' <<< "$B") | head -20
  fi
done
check "step 07 leaves no rewritten copy behind" test ! -e "${G}/${S2}/vcf/${S2}.pharmcat_input.vcf"

finish
