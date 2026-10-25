#!/usr/bin/env bash
# Step 05 (AnnotSV, the real image) keeps a long Delly call that Delly's own
# read depth (FORMAT/RDCN) contradicts out of AnnotSV's input, on the made-up
# AnnotSV data of tests/fixtures/annotsv/class45.sh, in a GENOME_DIR of its
# own. The calls describe no person:
#   Manta   DEL3 chr20:3000001-3450000 over H01-H40 (class 4);
#   Delly   ART  chr20:2000001-2450000 over G01-G40, PASS, RDCN 2: no loss
#                in read depth, so it is set aside;
#           KEEP chr20:4000001-4100000, PASS, RDCN 1: added;
#           SHORT chr20:5000001-5000500, PASS, RDCN 2: 500 bp, added.
# Checks: AnnotSV's input lacks ART and holds KEEP and SHORT with
# SVCALLER=delly; the log line and the ##svcallers_added header count one
# call set aside; the class 4 and 5 table holds DEL3 and not ART.
# Red control: the same Delly calls without the RDCN field (as a VCF not
# written by Delly 2.x) put ART in the input and in the class 4 and 5 table.
. "$(dirname "$0")/lib.sh"
# shellcheck source=../fixtures/annotsv/class45.sh
. "${REPO}/tests/fixtures/annotsv/class45.sh"

G="${CASE_TMP}/genome"
T="${SAMPLE}rdcn"
D="${G}/${T}"
rm -rf "$G"
class45_bundle "${G}/annotsv_annotations"
mkdir -p "${D}/manta/results/variants" "${D}/delly"

r=$'chr20\t3000001\tDEL3\tN\t<DEL>\t500\tPASS\tSVTYPE=DEL;END=3450000;SVLEN=-449999\tGT\t0/1\n'
_class45_vcf "$T" "$r" > "${D}/manta/results/variants/diploidSV.vcf"

# Delly's FORMAT, in the order Delly 2.x writes it.
dl() {  # dl POS ID END RDCN
  printf 'chr20\t%s\t%s\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=%s\tGT:GL:GQ:FT:RCL:RC:RCR:RDCN:DR:DV:RR:RV\t0/1:-10,0,-10:40:PASS:100:100:100:%s:10:10:0:0\n' \
    "$1" "$2" "$3" "$4"
}
{
  printf '##fileformat=VCFv4.2\n##contig=<ID=chr20,length=64444167>\n##FILTER=<ID=PASS,Description="All filters passed">\n'
  printf '##ALT=<ID=DEL,Description="Deletion">\n'
  printf '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
  printf '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods">\n'
  printf '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n'
  printf '##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">\n'
  printf '##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Raw high-quality read counts for the left control region">\n'
  printf '##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts for the SV">\n'
  printf '##FORMAT=<ID=RCR,Number=1,Type=Integer,Description="Raw high-quality read counts for the right control region">\n'
  printf '##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description="Read-depth based copy-number estimate for autosomal sites">\n'
  printf '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">\n'
  printf '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">\n'
  printf '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">\n'
  printf '##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$T"
  dl 2000001 ART 2450000 2
  dl 4000001 KEEP 4100000 1
  dl 5000001 SHORT 5000500 2
} > "${D}/delly/${T}_sv.vcf"
# shellcheck disable=SC2016  # $1 expands inside the container
docker run --rm --network none --user "$(id -u):$(id -g)" -v "${D}:/d" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
  for f in /d/manta/results/variants/diploidSV.vcf "/d/delly/$1_sv.vcf"; do
    bcftools view -Oz -o "$f.gz" "$f" && bcftools index -t "$f.gz" && rm "$f"
  done' _ "$T"

IN="${D}/annotsv/${T}_annotsv_input.vcf"
C45="${D}/annotsv/${T}_sv_class4_5.tsv"
delly_ids() { awk -F'\t' '!/^#/ && $8 ~ /SVCALLER=delly/ {print $3}' "$IN" 2>/dev/null | sort | tr '\n' ' '; }
class45_ids() {
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next } { print $h["ID"] }' "$C45" 2>/dev/null | sort -u | tr '\n' ' '
}

GENOME_DIR="$G" run_step 05-annotsv.sh "$T" male
check_step_exit 05-annotsv.sh
check_eq "AnnotSV's input: the Delly calls added" "$(delly_ids)" "KEEP SHORT "
check "the log counts the call set aside" has 'delly: 2 PASS DEL/DUP calls, 2 added, and 1 PASS calls longer than 1000 bp set aside' "$(cat "$STEP_LOG")"
check "the input's header counts it" has '^##svcallers_added=.*and 1 PASS calls longer than 1000 bp set aside' "$(cat "$IN" 2>/dev/null)"
check_eq "class 4 and 5: DEL3 and not ART" "$(class45_ids)" "DEL3 "

# Red control: the same calls with no RDCN field are read as before.
# shellcheck disable=SC2016  # $1 expands inside the container
docker run --rm --network none --user "$(id -u):$(id -g)" -v "${D}:/d" "$BCFTOOLS_IMAGE" bash -euo pipefail -c '
  f="/d/delly/$1_sv.vcf.gz"
  bcftools annotate -x FORMAT/RDCN -Oz -o "$f.tmp" "$f" && mv "$f.tmp" "$f" && bcftools index -f -t "$f"' _ "$T"
GENOME_DIR="$G" run_step 05-annotsv.sh "$T" male
check "red control: without RDCN ART is in AnnotSV's input" has 'ART' "$(delly_ids)"
check "red control: without RDCN ART reaches the class 4 and 5 table" has 'ART' "$(class45_ids)"
rm -rf "$G"

finish
