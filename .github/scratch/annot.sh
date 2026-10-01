#!/usr/bin/env bash
# Steps 30, 23, 31 and 12 on tiny synthetic inputs (pgp-9ms.1, .4, .6). Temporary.
# shellcheck source=lib.sh
. "$(dirname "$0")/lib.sh"
prelude

G=/tmp/an/gnew
GO=/tmp/an/gold
mkdir -p "$G/s1/vep" "$G/s1/vcf" "$G/annotations"

# --- A VEP-style VCF: two HIGH variants in GENE1 (a compound-het pair), one MODERATE in GENE2 ---
CSQ_FMT='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|Existing_variation|gnomADe_AF|CLIN_SIG'
{
  echo '##fileformat=VCFv4.2'
  echo '##FILTER=<ID=PASS,Description="All filters passed">'
  echo '##contig=<ID=chr20,length=64444167>'
  echo '##contig=<ID=chrM,length=16569>'
  echo "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ${CSQ_FMT}\">"
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t1000\t.\tA\tG\t50\tPASS\tCSQ=G|stop_gained|HIGH|GENE1|ENSG01|Transcript|ENST01|protein_coding|1/5||||rs1|0.0001|pathogenic\tGT\t0/1\n'
  printf 'chr20\t2000\t.\tC\tT\t50\tPASS\tCSQ=T|stop_gained|HIGH|GENE1|ENSG01|Transcript|ENST01|protein_coding|2/5|||||0.0002|\tGT\t0/1\n'
  printf 'chr20\t3000\t.\tG\tA\t50\tPASS\tCSQ=A|missense_variant|MODERATE|GENE2|ENSG02|Transcript|ENST02|protein_coding|3/7|||||0.0003|\tGT\t0/1\n'
} > "$G/s1/vep/s1_vep.vcf"

# --- Tiny annotation tracks: CADD (bare chromosome names) and a masked SpliceAI file ---
{
  echo '## CADD GRCh38-v1.7 (synthetic)'
  printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n'
  printf '20\t1000\tA\tG\t5.1\t35.0\n'
  printf '20\t3000\tG\tA\t3.2\t24.5\n'
} | bgzip -c > "$G/annotations/whole_genome_SNVs.tsv.gz"
tabix -s1 -b2 -e2 -c'#' "$G/annotations/whole_genome_SNVs.tsv.gz"
{
  echo '##fileformat=VCFv4.0'
  echo '##contig=<ID=chr20,length=64444167>'
  echo '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr20\t3000\t.\tG\tA\t.\t.\tSpliceAI=A|GENE2|0.50|0.00|0.00|0.00|1|2|3|4\n'
} | bgzip -c > "$G/annotations/spliceai_scores.masked.snv.hg38.vcf.gz"
tabix -p vcf "$G/annotations/spliceai_scores.masked.snv.hg38.vcf.gz"

# --- A chrM VCF for step 12: rCRS positions of a few common non-H markers ---
{
  echo '##fileformat=VCFv4.2'
  echo '##FILTER=<ID=PASS,Description="All filters passed">'
  echo '##contig=<ID=chr20,length=64444167>'
  echo '##contig=<ID=chrM,length=16569>'
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n'
  printf 'chr20\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n'
  for v in 73:A:G 263:A:G 750:A:G 1438:A:G 2706:A:G 4769:A:G 7028:C:T 8860:A:G 11719:G:A 14766:C:T 15326:A:G; do
    IFS=: read -r p r a <<< "$v"
    printf 'chrM\t%s\t.\t%s\t%s\t50\tPASS\t.\tGT\t1\n' "$p" "$r" "$a"
  done
} | bgzip -c > "$G/s1/vcf/s1.vcf.gz"
tabix -p vcf "$G/s1/vcf/s1.vcf.gz"
cp -r "$G" "$GO"

# --- Step 30 ---
expect_ok "30 new: vcfanno with CADD + masked SpliceAI" env GENOME_DIR="$G" bash "$NEW/scripts/30-vcfanno.sh" s1
OUTF="$G/s1/vep/s1_annotated.vcf.gz"
N=$(bcftools view -H "$OUTF" 2>/dev/null | wc -l)
check "30 new: annotated VCF records" "$N" "3" "$([ "$N" = 3 ] && echo 1 || echo 0)"
check "30 new: annotated VCF index" "$(ls "$OUTF.tbi" 2>&1)" "present" "$([ -f "$OUTF.tbi" ] && echo 1 || echo 0)"
H=$(bcftools view -h "$OUTF" | grep -oE 'ID=(CADD_PHRED|SpliceAI),' | tr '\n' ' ')
check "30 new: CADD_PHRED and SpliceAI in header" "$H" "both" "$(grep -q CADD_PHRED <<< "$H" && grep -q SpliceAI <<< "$H" && echo 1 || echo 0)"
grep -E 'masked' "$LOGS/30_new:_vcfanno_with_CADD_+_masked_SpliceAI.log" || true
expect_ok "30 new: second run skips the complete output" env GENOME_DIR="$G" bash "$NEW/scripts/30-vcfanno.sh" s1
check "30 new: second run says skipping" "$(grep -c 'Skipping' "$LOGS/30_new:_second_run_skips_the_complete_output.log")" "1" "$(grep -q 'Skipping' "$LOGS/30_new:_second_run_skips_the_complete_output.log" && echo 1 || echo 0)"
expect_fail "30 old (origin-main, red-first): vcfanno" env GENOME_DIR="$GO" bash "$OLD/scripts/30-vcfanno.sh" s1
SZ=$(stat -c %s "$GO/s1/vep/s1_vep.vcf.gz" 2>/dev/null || echo missing)
check "30 old (red-first): _vep.vcf.gz size in bytes" "$SZ" "0" "$([ "$SZ" = 0 ] && echo 1 || echo 0)"

# --- Step 23 on step 30's output ---
expect_ok "23 new: clinical filter" env GENOME_DIR="$G" bash "$NEW/scripts/23-clinical-filter.sh" s1
N=$(bcftools view -H "$G/s1/clinical/s1_spliceai_high.vcf.gz" 2>/dev/null | wc -l)
check "23 new: SpliceAI-high records" "$N" "1" "$([ "$N" = 1 ] && echo 1 || echo 0)"
N=$(bcftools view -H "$G/s1/clinical/s1_clinical.vcf.gz" 2>/dev/null | wc -l)
check "23 new: clinical VCF records" "$N" ">0" "$([ "$N" -gt 0 ] && echo 1 || echo 0)"
# Red-first: origin/main step 23 on the same (good) input still calls bgzip
cp "$G/s1/vep/s1_annotated.vcf.gz" "$G/s1/vep/s1_annotated.vcf.gz.tbi" "$GO/s1/vep/"
expect_fail "23 old (origin-main, red-first): clinical filter" env GENOME_DIR="$GO" bash "$OLD/scripts/23-clinical-filter.sh" s1
grep -m2 -E 'bgzip|tabix' "$LOGS/23_old_(origin-main,_red-first):_clinical_filter.log" || true

# --- Step 31 ---
expect_ok "31 new: slivar" env GENOME_DIR="$G" bash "$NEW/scripts/31-slivar.sh" s1
grep -E 'Found: .* compound het|no compound heterozygote' "$LOGS/31_new:_slivar.log" || true
check "31 new: summary TSV has rows" "$(wc -l < "$G/s1/slivar/s1_slivar_summary.tsv")" ">1 line" "$([ "$(wc -l < "$G/s1/slivar/s1_slivar_summary.tsv")" -gt 1 ] && echo 1 || echo 0)"
# Wrong slivar image: must fail, with docker's own error
rm -rf /tmp/badslivar && cp -r "$NEW" /tmp/badslivar
sed -i 's#^SLIVAR_IMAGE=.*#SLIVAR_IMAGE="quay.io/biocontainers/slivar:0.0.0--doesnotexist"#' /tmp/badslivar/versions.env
expect_fail "31 new with a wrong slivar image" env GENOME_DIR="$G" bash /tmp/badslivar/scripts/31-slivar.sh s1
grep -m3 -E 'ERROR|manifest|not found|denied' "$LOGS/31_new_with_a_wrong_slivar_image.log" || true
# Red-first: origin/main step 31 with the same wrong image exits 0 and says "no candidates"
cp "$G/s1/vep/s1_annotated.vcf.gz" "$G/s1/vep/s1_annotated.vcf.gz.tbi" "$GO/s1/vep/"
rm -rf /tmp/badold && cp -r "$OLD" /tmp/badold
echo 'SLIVAR_IMAGE="quay.io/biocontainers/slivar:0.0.0--doesnotexist"' >> /tmp/badold/versions.env
observe "31 old (origin-main, red-first) with a wrong slivar image" env GENOME_DIR="$GO" bash /tmp/badold/scripts/31-slivar.sh s1
L="$LOGS/31_old_(origin-main,_red-first)_with_a_wrong_slivar_image.log"
grep -E 'No compound heterozygote candidates found|WARNING: slivar' "$L" || true
check "31 old (red-first): wrong image reported as no candidates" "$(grep -c 'No compound heterozygote candidates found' "$L")" "1" "$(grep -q 'No compound heterozygote candidates found' "$L" && echo 1 || echo 0)"
expect_fail "31 old (origin-main, red-first) with origin-main versions.env" env GENOME_DIR="$GO" bash "$OLD/scripts/31-slivar.sh" s1
grep -m1 'unbound variable' "$LOGS/31_old_(origin-main,_red-first)_with_origin-main_versions.env.log" || true

# --- Step 12 ---
expect_ok "12 new: haplogrep3 classify" env GENOME_DIR="$G" bash "$NEW/scripts/12-mito-haplogroup.sh" s1
HG="$G/s1/mito/s1_haplogroup.txt"
cat "$HG" 2>/dev/null | head -3
COL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){h=$i; gsub(/"/,"",h); if(h=="Haplogroup") c=i}} NR==2 && c {v=$c; gsub(/"/,"",v); print v}' "$HG" 2>/dev/null)
check "12 new: haplogroup column" "${COL:-<empty>}" "non-empty" "$([ -n "$COL" ] && echo 1 || echo 0)"
expect_ok "12 new: second run" env GENOME_DIR="$G" bash "$NEW/scripts/12-mito-haplogroup.sh" s1
expect_fail "12 old (origin-main, red-first): haplogrep3" env GENOME_DIR="$GO" bash "$OLD/scripts/12-mito-haplogroup.sh" s1
grep -m2 -E 'executable file not found|classify' "$LOGS/12_old_(origin-main,_red-first):_haplogrep3.log" || true

finish
