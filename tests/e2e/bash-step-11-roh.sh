#!/usr/bin/env bash
# Step 11 and the Nextflow ROH module (run by 60-nextflow on the same VCF) read
# the same records, PASS and unfiltered only, and write the same summary.
. "$(dirname "$0")/lib.sh"

run_step 11-roh-analysis.sh "$SAMPLE"
check_step_exit 11-roh-analysis.sh
B="${GENOME_DIR}/${SAMPLE}/vcf"
N="${GENOME_DIR}/nf-results/${SAMPLE}/roh"
check "the bash step wrote a summary" test -s "${B}/${SAMPLE}_roh_summary.txt"
check "the module wrote a summary" test -s "${N}/${SAMPLE}_roh_summary.txt"
check "bash and module summaries are equal" cmp -s "${B}/${SAMPLE}_roh_summary.txt" "${N}/${SAMPLE}_roh_summary.txt"
diff "${B}/${SAMPLE}_roh_summary.txt" "${N}/${SAMPLE}_roh_summary.txt" | head -n 10
grep '^RG' "${B}/${SAMPLE}_roh.txt" 2>/dev/null | sort > "${CASE_TMP}/rg-bash.txt"
grep '^RG' "${N}/${SAMPLE}_roh.txt" 2>/dev/null | sort > "${CASE_TMP}/rg-module.txt"
check_ge "ROH segments (RG lines) from the bash step" "$(grep -c . "${CASE_TMP}/rg-bash.txt" || true)" 1
check "bash and module find the same segments" cmp -s "${CASE_TMP}/rg-bash.txt" "${CASE_TMP}/rg-module.txt"

# No per-site line at a record DeepVariant did not PASS (RefCall and others).
bcf query -i 'FILTER!="PASS" && FILTER!="."' -f '%CHROM\t%POS\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
  | sort -u > "${CASE_TMP}/filtered.txt"
check_ge "filtered records in the VCF (the check needs some)" "$(grep -c . "${CASE_TMP}/filtered.txt" || true)" 1
awk -F'\t' '$1 == "ST" {print $3 "\t" $4}' "${B}/${SAMPLE}_roh.txt" | sort -u > "${CASE_TMP}/sites.txt"
check_eq "per-site lines at filtered records" "$(comm -12 "${CASE_TMP}/filtered.txt" "${CASE_TMP}/sites.txt" | grep -c . || true)" 0

finish
