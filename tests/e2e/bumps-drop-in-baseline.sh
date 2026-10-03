#!/usr/bin/env bash
# What a drop-in tool bump must not change: record counts of steps 06, 11 and
# 14, the chrM call set of step 20 and the HLA alleles of step 08, compared
# with the values measured on this fixture with the images before package
# bumps-drop-in (samtools 1.20, bcftools 1.21, GATK 4.6.2.0, T1K 1.0.9).
#
# The case reads what earlier cases left: 30-clinvar-screen-06 (step 06),
# bash-step-11-roh (step 11), bash-step-14-imputation (step 14),
# bash-step-20-numt (step 20) and bash-step-08-hla (step 08, IPD-IMGT/HLA
# 3.64.0). A bump that changes a value fails here and the diff names the line.
# When a change is expected (a tool fix, a new fixture), update BASELINE in
# the same pull request and say why in its description.
. "$(dirname "$0")/lib.sh"

BASELINE=$(cat <<'BASE'
step06 pass.vcf.gz records 27744
step06 clinvar_hits.vcf records 1
step11 ST lines 27021
step11 RG lines 38
step14 chr1 records 242
step14 chr2 records 239
step14 chr4 records 508
step14 chr5 records 284
step14 chr6 records 22264
step14 chr10 records 130
step14 chr12 records 226
step14 chr16 records 70
step14 chr19 records 577
step14 chr20 records 1061
step14 chr22 records 577
step14 chrX records 1161
step20 chrM records 31
step20 chrM PASS records 19
step20 chrM call set md5 6015b376c274d306b5dd66d3571c62ee
step08 HLA-A HLA-A*01:01 HLA-A*26:01
step08 HLA-B HLA-B*35:08 HLA-B*38:01
step08 HLA-C HLA-C*04:01 HLA-C*12:03
BASE
)

S=$SAMPLE
OUT="${CASE_TMP}/measured.txt"
{
  # The isec intermediate (_shared.vcf.gz) is deleted by the step itself.
  for f in pass.vcf.gz clinvar_hits.vcf; do
    echo "step06 ${f} records $(vcf_count "${S}/clinvar/${S}_${f}")"
  done
  for k in ST RG; do
    echo "step11 ${k} lines $(grep -c "^${k}" "${GENOME_DIR}/${S}/vcf/${S}_roh.txt" 2>/dev/null || true)"
  done
  for f in "${GENOME_DIR}/${S}/imputation/mis_ready/${S}"_chr*.vcf.gz; do
    [ -e "$f" ] || continue
    c=${f##*_}; c=${c%.vcf.gz}
    echo "step14 ${c} records $(vcf_count "${f#"${GENOME_DIR}/"}")"
  done | sort -V
  M="${S}/mito/${S}_chrM_filtered.vcf.gz"
  bcf query -f '%CHROM:%POS:%REF:%ALT:%FILTER\n' "$M" 2>/dev/null > "${CASE_TMP}/chrM.txt"
  echo "step20 chrM records $(grep -c . "${CASE_TMP}/chrM.txt" || true)"
  echo "step20 chrM PASS records $(grep -c ':PASS$' "${CASE_TMP}/chrM.txt" || true)"
  echo "step20 chrM call set md5 $(md5sum < "${CASE_TMP}/chrM.txt" | cut -c1-32)"
  # Alleles at two-field resolution, the two of a gene in sorted order.
  awk -F'\t' '$1 ~ /^HLA-[ABC]$/ {
      split($3, x, ":"); a = x[1] ":" x[2]; b = "-"
      if ($6 ~ /^HLA-/) { split($6, y, ":"); b = y[1] ":" y[2] }
      if (b != "-" && b < a) { t = a; a = b; b = t }
      print "step08 " $1 " " a " " b
    }' "${GENOME_DIR}/${S}/hla_t1k/${S}_hla_genotype.tsv" 2>/dev/null | sort
} > "$OUT"

echo "Measured (copy into BASELINE when a change is expected):"
cat "$OUT"
echo "chrM calls:"
cat "${CASE_TMP}/chrM.txt"

check_ge "values measured" "$(grep -c . "$OUT" || true)" 15
diff <(printf '%s\n' "$BASELINE") "$OUT" > "${CASE_TMP}/diff.txt"
cat "${CASE_TMP}/diff.txt"
check_eq "lines that differ from the baseline" "$(grep -c '^[<>]' "${CASE_TMP}/diff.txt" || true)" 0
check "HLA-A reports both alleles of HG002 (A*01:01 and A*26:01)" \
  grep -qx 'step08 HLA-A HLA-A\*01:01 HLA-A\*26:01' "$OUT"

finish
