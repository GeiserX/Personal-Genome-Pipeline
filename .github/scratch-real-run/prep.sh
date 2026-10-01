#!/usr/bin/env bash
# Scratch helper (deleted before merge): build a small real HG002 input set
# on chr20 from public GIAB, UCSC, NCBI and Illumina data.
set -euo pipefail

D=${1:?usage: prep.sh <data_dir>}
mkdir -p "$D"
cd "$D"

REGION=chr20:2000000-3200000
BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
VCF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
CLINVAR_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
CATALOG_URL=https://raw.githubusercontent.com/Illumina/ExpansionHunter/v5.0.0/variant_catalog/grch38/variant_catalog.json

echo "== reference: chr20 only (UCSC hg38, same sequence as the analysis set)"
curl -fsSL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz | gunzip > ref.fa
samtools faidx ref.fa
cut -f1,2 ref.fa.fai

echo "== BAM slice ${REGION}, header cut to chr20, mates elsewhere dropped"
samtools view -h "$BAM_URL" "$REGION" \
  | awk -F'\t' '/^@SQ/ { if ($2 == "SN:chr20") print; next } /^@/ { print; next } $7 == "=" || $7 == "*" || $7 == "chr20" { print }' \
  | samtools view -b -o HG002.bam -
samtools index HG002.bam
samtools idxstats HG002.bam | awk '$3 > 0'

echo "== VCF slice (GIAB v4.2.1 benchmark)"
bcftools view -r "$REGION" -m2 -M2 "$VCF_URL" -Oz -o HG002.vcf.gz
bcftools index -t HG002.vcf.gz
echo "records: $(bcftools view -H HG002.vcf.gz | wc -l)"
echo "FILTER values:"
bcftools query -f '%FILTER\n' HG002.vcf.gz | sort | uniq -c

echo "== same VCF with every FILTER set to '.'"
bcftools annotate -x FILTER HG002.vcf.gz -Oz -o HG002_nofilter.vcf.gz
bcftools index -t HG002_nofilter.vcf.gz
bcftools query -f '%FILTER\n' HG002_nofilter.vcf.gz | sort | uniq -c

echo "== ClinVar P/LP on the slice, chr-prefixed"
echo "20 chr20" > chr_rename.txt
bcftools view -r 20:2000000-3200000 "$CLINVAR_URL" -Ou \
  | bcftools annotate --rename-chrs chr_rename.txt -Ou \
  | bcftools view -i 'CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic"' -Oz -o clinvar_pathogenic_chr.vcf.gz
bcftools index -t clinvar_pathogenic_chr.vcf.gz
echo "ClinVar P/LP records: $(bcftools view -H clinvar_pathogenic_chr.vcf.gz | wc -l)"

echo "== tiny CADD-format score file (bare chromosome names, like CADD)"
{
  echo "## CADD GRCh38-v1.7 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health at Charite 2013-2023. All rights reserved."
  printf '#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n'
  bcftools query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%REF\t%ALT\n' HG002.vcf.gz \
    | awk -v OFS='\t' 'NR % 3 == 1 { sub(/^chr/, "", $1); print $1, $2, $3, $4, "1.5", "25.5" }'
} | bgzip > cadd_tiny.tsv.gz
tabix -s 1 -b 2 -e 2 -c '#' cadd_tiny.tsv.gz
echo "CADD rows: $(zcat cadd_tiny.tsv.gz | grep -vc '^#')"

echo "== tiny SpliceAI-format VCF (chr names, the masked file name)"
{
  echo '##fileformat=VCFv4.0'
  echo '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">'
  echo '##contig=<ID=chr20,length=64444167>'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  bcftools query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%REF\t%ALT\n' HG002.vcf.gz \
    | awk -v OFS='\t' 'NR % 5 == 2 { print $1, $2, ".", $3, $4, ".", ".", "SpliceAI=" $4 "|TESTGENE|0.91|0.00|0.00|0.00|1|-2|3|-4" }'
} | bgzip > spliceai_scores.masked.snv.hg38.vcf.gz
tabix -p vcf spliceai_scores.masked.snv.hg38.vcf.gz
echo "SpliceAI rows: $(zcat spliceai_scores.masked.snv.hg38.vcf.gz | grep -vc '^#')"

echo "== ExpansionHunter catalog: NOP56 (chr20) only"
curl -fsSL "$CATALOG_URL" \
  | jq '[ .[] | select(.LocusId == "NOP56") | .ReferenceRegion |= ("chr" + .) ]' > eh_catalog_chr20.json
cat eh_catalog_chr20.json

echo "== Delly exclude map (chr20 p-arm telomere)"
printf 'chr20\t0\t10000\ttelomere\n' > delly_exclude.tsv

echo "== samplesheets"
printf 'sample,vcf,vcf_index,bam,bam_index,sex\nHG002,%s/HG002.vcf.gz,%s/HG002.vcf.gz.tbi,%s/HG002.bam,%s/HG002.bam.bai,male\n' "$D" "$D" "$D" "$D" > samplesheet.csv
printf 'sample,vcf,vcf_index,bam,bam_index\nHG002,%s/HG002.vcf.gz,%s/HG002.vcf.gz.tbi,%s/HG002.bam,%s/HG002.bam.bai\n' "$D" "$D" "$D" "$D" > samplesheet_nosex.csv
{ cat samplesheet.csv; tail -n 1 samplesheet.csv; } > samplesheet_dup.csv
printf 'sample,vcf,vcf_index,bam,bam_index,sex\nHG002_nofilter,%s/HG002_nofilter.vcf.gz,%s/HG002_nofilter.vcf.gz.tbi,%s/HG002.bam,%s/HG002.bam.bai,male\n' "$D" "$D" "$D" "$D" > samplesheet_nofilter.csv
cat samplesheet.csv samplesheet_dup.csv

cat > ci.config <<'CICONF'
trace.enabled = false
dag.enabled = false
timeline.enabled = false
report.enabled = false
docker.runOptions = '-u 0:0'
process.resourceLimits = [ cpus: 4, memory: '14.GB', time: '2.h' ]
CICONF

ls -la "$D"
