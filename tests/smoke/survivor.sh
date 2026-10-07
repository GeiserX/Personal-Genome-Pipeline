# shellcheck shell=sh
# SURVIVOR_IMAGE row: two callers, the same 2 kb deletion called 2 bp apart
# across a 1 kb boundary, and two deletions of 600 bp and 5 kb that start in
# one 1 kb window. Step 22's parameters must keep the first, once, with
# support 2, and drop the second pair (tests/fixtures/sv/ has the three-caller set).
set -e
head='##fileformat=VCFv4.2
##contig=<ID=chr20,length=64444167>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
printf '%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\n' "$head" > a.vcf
printf 'chr20\t10100999\ta1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=10102999;SVLEN=-2000\tGT\t0/1\n' >> a.vcf
printf 'chr20\t10200100\ta2\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=10200700;SVLEN=-600\tGT\t0/1\n' >> a.vcf
printf '%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdelly\n' "$head" > b.vcf
printf 'chr20\t10101001\tb1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=10103001;SVLEN=-2000\tGT\t0/1\n' >> b.vcf
printf 'chr20\t10200300\tb2\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=10205300;SVLEN=-5000\tGT\t0/1\n' >> b.vcf
printf '%s\n' "$PWD/a.vcf" "$PWD/b.vcf" > list.txt
SURVIVOR merge list.txt 1000 2 1 1 0 50 merged.vcf
grep -v '^##' merged.vcf
