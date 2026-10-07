# shellcheck shell=sh
# VEP_IMAGE row for step 13's --custom ClinVar annotation: a one-record
# ClinVar-style VCF (the first record of the 50, CLNSIG Pathogenic) is
# annotated the way step 13 and the VEP module pass the ClinVar file, and
# its CLNSIG must come back as the CSQ field ClinVar_CLNSIG of that record.
set -e
{
  printf '##fileformat=VCFv4.2\n'
  printf '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n'
  printf '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">\n'
  printf '##INFO=<ID=CLNDN,Number=.,Type=String,Description="Disease name">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk -F'\t' 'BEGIN {OFS = "\t"} /^#/ {next} !n++ {print $1, $2, "900000001", $4, $5, ".", ".", "CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Smoke_test"}' /in/sample50.vcf
} > cv.vcf
bgzip -f cv.vcf && tabix -f -p vcf cv.vcf.gz
awk -F'\t' '/^#/ {next} {print $1 ":" $2; exit}' /in/sample50.vcf > first.txt
vep --input_file /in/sample50.vcf --output_file vep_cv.vcf --vcf --database --assembly GRCh38 --symbol \
  --force_overwrite --no_stats \
  --custom file=/out/cv.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
grep -m1 '^##INFO=<ID=CSQ' vep_cv.vcf
