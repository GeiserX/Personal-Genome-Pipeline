#!/usr/bin/env bash
# The commands of docs/vcf-first.md and the header recipe of docs/nextflow.md
# ("Before you share outputs") run as pasted on a vendor-style copy of the
# fixture VCF: Ensembl contig names, gVCF reference blocks, and header lines
# that name the sample the way callers and providers write them. Only the
# page's Setup block is replaced, by the same four variables pointing here.
. "$(dirname "$0")/lib.sh"
. "$(dirname "$0")/vendor-vcf-intake.inc"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }

# bash_blocks FILE: every ```bash block of FILE, in order, separated by a line
# "#--- block N".
bash_blocks() {
  awk '/^```bash$/ {n++; print "#--- block " n; on = 1; next} /^```/ {on = 0} on' "$1"
}

# --- A vendor-style file ----------------------------------------------------------
# Ensembl names, a ##GVCFBlock line, two reference blocks, and four header lines
# that carry the sample id the way different tools write it.
V="${G}/vendor"
mkdir -p "$V"
{
  bcf view -h "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^##' | grep -v '^##INFO=<ID=END,'
  echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the reference block">'
  echo '##GVCFBlock0-20=minGQ=0(inclusive),maxGQ=20(exclusive)'
  echo "##bcftoolsCommand=mpileup -f ref.fa ${SAMPLE}_vendor.bam"
  echo "##commandline=\"caller --sample ${SAMPLE} --out ${SAMPLE}.genome.vcf\""
  echo "##cmdline=caller ${SAMPLE}"
  echo "##GATKCommandLine=<ID=HaplotypeCaller,CommandLine=\"HaplotypeCaller -I ${SAMPLE}.bam -O ${SAMPLE}.g.vcf.gz\">"
  bcf view -h "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^#CHROM'
  bcf view -H "${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  printf 'chr20\t10000001\t.\tN\t<*>\t0\t.\tEND=10000050\tGT\t0/0\n'
  printf 'chr10\t94700001\t.\tN\t<NON_REF>\t.\t.\tEND=94700100\tGT\t0/0\n'
} > "${V}/unsorted.vcf"
in_genome "$BCFTOOLS_IMAGE" sh -c "set -e
  bcftools sort -Ou vendor/unsorted.vcf | bcftools annotate --no-version --rename-chrs intake/to_ensembl.txt \
    -Oz -o vendor/${SAMPLE}.vendor.g.vcf.gz
  bcftools index -f -t vendor/${SAMPLE}.vendor.g.vcf.gz"
rm -f "${V}/unsorted.vcf" 2>/dev/null
VENDOR="${V}/${SAMPLE}.vendor.g.vcf.gz"
check "the vendor-style file has Ensembl names" lacks '^chr' "$(bcf index -s "vendor/${SAMPLE}.vendor.g.vcf.gz" | cut -f1)"
check_ge "the vendor-style file has header lines naming the sample" \
  "$(gzip -dc "$VENDOR" | grep '^##' | grep -c "$SAMPLE" || true)" 4

# --- docs/vcf-first.md, as pasted ----------------------------------------------------
PAGE="${REPO}/docs/vcf-first.md"
bash_blocks "$PAGE" > "${CASE_TMP}/page.sh"
check "the page's first bash block is its Setup block" \
  has '^# Setup:' "$(awk '/^#--- block 1$/ {on = 1; next} /^#--- block/ {on = 0} on' "${CASE_TMP}/page.sh")"
{
  echo 'set -euo pipefail'
  printf 'PGP=%q\nREF=%q\nVCF=%q\nLABEL=%q\n' "$REPO" "$REF" "$VENDOR" sample1
  awk '/^#--- block 1$/ {skip = 1; next} /^#--- block/ {skip = 0} !skip' "${CASE_TMP}/page.sh"
} > "${CASE_TMP}/run-page.sh"
echo "+ the page's commands:"; sed 's/^/    /' "${CASE_TMP}/run-page.sh"
(cd "$CASE_TMP" && bash "${CASE_TMP}/run-page.sh") > "${CASE_TMP}/run-page.log" 2>&1
PAGE_RC=$?
grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "${CASE_TMP}/run-page.log" | sed 's/^/    | /'
check_eq "the page's commands run as pasted and exit 0" "$PAGE_RC" 0

OUT="${V}/sample1.vcf.gz"
check "the page wrote sample1.vcf.gz and its index" test -s "$OUT" -a -s "${OUT}.tbi"
check "every contig holding records is chr-named" \
  lacks '^[^c]' "$(bcf index -s vendor/sample1.vcf.gz | cut -f1)"
check_eq "no reference block left" "$(vcf_count -i 'INFO/END!="."' vendor/sample1.vcf.gz)" 0
check_eq "same records as the fixture VCF" \
  "$(bcf query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n' vendor/sample1.vcf.gz | md5sum)" \
  "$(bcf query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n' "${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | md5sum)"
check_eq "the sample column is the label" "$(bcf query -l vendor/sample1.vcf.gz)" sample1
HDR=$(gzip -dc "$OUT" | grep '^#')
check_eq "header lines outside the kept list" \
  "$(grep -vcE '^(##(fileformat|FILTER|INFO|FORMAT|ALT|contig)=|#CHROM)' <<< "$HDR" || true)" 0
check_eq "header lines naming the sample" "$(grep -c "$SAMPLE" <<< "$HDR" || true)" 0
for k in bcftoolsCommand commandline cmdline GATKCommandLine; do
  check "no ##${k} line" lacks "^##${k}=" "$HDR"
done
R="${V}/results/sample1"
check "the page's Nextflow run wrote the report" test -s "${R}/sample1_report.html"
check_eq "its ROH card is filled" "$(html_stat "${R}/sample1_report.html" 'Runs of Homozygosity' Status)" Complete
check_eq "its haplogroup card is filled" "$(html_stat "${R}/sample1_report.html" 'Mitochondrial Haplogroup' Status)" Complete
JSON="${R}/pharmcat/sample1.report.json"
check_ge "its PharmCAT calls genes" "$( [ -s "$JSON" ] && pharmcat_called "$JSON" | grep -c . || echo 0)" 1

# --- docs/nextflow.md "Before you share outputs", as pasted ------------------------------
RECIPE=$(awk '/^## Before you share outputs/ {on = 1} on && /^```bash$/ {b = 1; next} b && /^```/ {exit} b' \
  "${REPO}/docs/nextflow.md")
check "the recipe block was found" test -n "$RECIPE"
S="${CASE_TMP}/share"
mkdir -p "$S"
cp "$VENDOR" "${S}/in.vcf.gz"
cp "${VENDOR}.tbi" "${S}/in.vcf.gz.tbi"
(
  cd "$S" || exit 1
  eval "$bcftools_fn"
  set -euo pipefail
  eval "$RECIPE"
) > "${CASE_TMP}/share.log" 2>&1
check_eq "the recipe runs as pasted and exits 0" "$?" 0
cat "${CASE_TMP}/share.log"
HDR=$(gzip -dc "${S}/SAMPLE.vcf.gz" 2>/dev/null | grep '^#')
for k in bcftoolsCommand commandline cmdline GATKCommandLine; do
  check "recipe: no ##${k} line" lacks "^##${k}=" "$HDR"
done
check_eq "recipe: header lines naming the old sample" "$(grep -c "$SAMPLE" <<< "$HDR" || true)" 0
check "recipe: the sample column is renamed" has $'^#CHROM\t.*\tSAMPLE$' "$HDR"
check_eq "recipe: records unchanged" \
  "$(gzip -dc "${S}/SAMPLE.vcf.gz" | grep -v '^#' | md5sum)" "$(gzip -dc "${S}/in.vcf.gz" | grep -v '^#' | md5sum)"

finish
