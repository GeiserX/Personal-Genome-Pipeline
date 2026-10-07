#!/usr/bin/env bash
# Step 30's REVEL checks, on the fixture's VEP output under a second sample
# name, with REVEL tables built here from the fixture's own missense SNVs:
#   1. the documented layout ('#chr pos ref alt REVEL'): the step exits 0, the
#      check sees the listed sites and each comes back with INFO/REVEL;
#   2. a header copied from the original 9-column REVEL file over 5 data
#      columns: the step stops before vcfanno, naming the header;
#   3. the same wrong header without '#', indexed with tabix -S 1 (vcfanno
#      reads it as the header, the header check does not see it): vcfanno
#      matches nothing, and the data check stops the step.
# The fixture's synthetic REVEL table is put back afterwards.
. "$(dirname "$0")/lib.sh"

if ! command -v tabix >/dev/null || ! command -v bgzip >/dev/null; then
  sudo apt-get update -q >/dev/null && sudo apt-get install -y -q tabix >/dev/null
fi
check "bgzip and tabix are available" command -v tabix

S2="${SAMPLE}revel"
ANN="${GENOME_DIR}/annotations"
REVEL="${ANN}/revel_grch38.tsv.gz"
mv "$REVEL" "${CASE_TMP}/revel.orig.tsv.gz"
mv "${REVEL}.tbi" "${CASE_TMP}/revel.orig.tsv.gz.tbi"
restore() {
  mv -f "${CASE_TMP}/revel.orig.tsv.gz" "$REVEL"
  mv -f "${CASE_TMP}/revel.orig.tsv.gz.tbi" "${REVEL}.tbi"
}
trap restore EXIT

mkdir -p "${GENOME_DIR}/${S2}/vep"
cp "${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf" "${GENOME_DIR}/${S2}/vep/${S2}_vep.vcf"

# Biallelic missense SNVs of the fixture, as "chrom pos ref alt"
SITES=$(in_genome "$BCFTOOLS_IMAGE" sh -c "bcftools view -Oz -o /tmp/in.vcf.gz ${S2}/vep/${S2}_vep.vcf && bcftools index -t /tmp/in.vcf.gz &&
  bcftools view -m2 -M2 -v snps -i 'INFO/CSQ~\"missense_variant\"' /tmp/in.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'" | head -5)
printf 'missense SNVs used:\n%s\n' "${SITES:-none}"
check_ge "missense SNVs in the fixture's VEP output" "$(grep -c . <<< "$SITES" || true)" 1

# table HEADER_LINE TABIX_ARGS...: write the REVEL table with these sites (score 0.7)
table() {
  local header=$1; shift
  rm -f "$REVEL" "${REVEL}.tbi"
  { printf '%s\n' "$header"; awk -F'\t' 'BEGIN {OFS = "\t"} {print $1, $2, $3, $4, "0.7"}' <<< "$SITES"; } \
    | bgzip -c > "$REVEL"
  tabix -f "$@" "$REVEL"
}
rerun() {
  rm -f "${GENOME_DIR}/${S2}/vep/${S2}_annotated.vcf.gz" "${GENOME_DIR}/${S2}/vep/${S2}_annotated.vcf.gz.tbi"
  run_step 30-vcfanno.sh "$S2"
}

# 1. the documented layout
table $'#chr\tpos\tref\talt\tREVEL' -s 1 -b 2 -e 2
rerun
check_step_exit 30-vcfanno.sh
N=$(grep -o 'Sites the REVEL table lists: [0-9]*' "$STEP_LOG" | grep -o '[0-9]*$')
check_ge "the check found the table's sites in the output" "${N:-0}" 1
check "and none came back without INFO/REVEL" grep -q 'without INFO/REVEL after vcfanno: 0$' "$STEP_LOG"
check_ge "annotated records carrying REVEL" "$(vcf_count -i 'INFO/REVEL!="."' "${S2}/vep/${S2}_annotated.vcf.gz")" 1

# 2. a 9-column header over 5 data columns, as a '#' line
WRONG=$'chr\thg19_pos\tgrch38_pos\tref\talt\taaref\taaalt\tREVEL\tEnsembl_transcriptid'
table "#${WRONG}" -s 1 -b 2 -e 2
rerun
check "a wrong '#' header makes step 30 exit non-zero" test "$STEP_RC" -ne 0
check "it names the header" grep -q 'has the header line' "$STEP_LOG"
check "it leaves no annotated VCF" test ! -e "${GENOME_DIR}/${S2}/vep/${S2}_annotated.vcf.gz"

# 3. the same header without '#', skipped by tabix -S 1
table "$WRONG" -s 1 -b 2 -e 2 -S 1
rerun
check "a wrong header vcfanno reads but the header check cannot see makes step 30 exit non-zero" test "$STEP_RC" -ne 0
check "the data check names the sites that came back without INFO/REVEL" grep -q 'came back without INFO/REVEL' "$STEP_LOG"

finish
