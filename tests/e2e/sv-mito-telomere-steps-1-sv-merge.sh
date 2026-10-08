#!/usr/bin/env bash
# Step 22 (SURVIVOR merge) on the synthetic three-caller set of
# tests/fixtures/sv/ (its README has the known answer): one 2 kb deletion
# called 2 bp apart across a 1 kb boundary must give one record with support
# 2; two deletions of 600 bp and 5 kb starting in the same 1 kb window must
# not merge; the call one caller made must not reach the consensus. Then a
# GRIDSS VCF beside them is left out of the consensus, and the log says so.
. "$(dirname "$0")/lib.sh"

T="${SAMPLE}sv"
FIX="${REPO}/tests/fixtures/sv"
rm -rf "${GENOME_DIR:?}/${T}"
mkdir -p "${GENOME_DIR}/${T}/manta/results/variants" "${GENOME_DIR}/${T}/delly" "${GENOME_DIR}/${T}/cnvpytor" \
  "${GENOME_DIR}/${T}/sv_gridss" "${CASE_TMP}/in"
cp "${FIX}/manta.vcf.in" "${FIX}/delly.vcf.in" "${FIX}/cnvpytor.vcf.in" "${CASE_TMP}/in/"
bgz() {  # bgz NAME DEST: the fixture VCF NAME, bgzipped and indexed at GENOME_DIR/DEST
  docker run --rm -i --user "$(id -u):$(id -g)" -v "${CASE_TMP}/in:/in:ro" -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
    sh -c "bcftools view -Oz -o /genome/$2 /in/$1 && bcftools index -t /genome/$2"
}
bgz manta.vcf.in "${T}/manta/results/variants/diploidSV.vcf.gz"
bgz delly.vcf.in "${T}/delly/${T}_sv.vcf.gz"
bgz cnvpytor.vcf.in "${T}/cnvpytor/${T}_cnvs.vcf.gz"
run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
VCF="${T}/sv_merged/${T}_sv_consensus.vcf.gz"
check "consensus VCF is readable" vcf_ok "$VCF"
echo "Consensus records (CHROM POS END SVTYPE, then SUPP and SUPP_VEC when the header has them):"
# The positions first, from fields every consensus VCF has, so a VCF without
# SUPP still lists its records; SUPP and SUPP_VEC are pasted on when declared.
bcf query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' "$VCF" 2>/dev/null > "${CASE_TMP}/pos.tsv"
if bcf view -h "$VCF" 2>/dev/null | grep -q '^##INFO=<ID=SUPP_VEC,'; then
  bcf query -f '%INFO/SUPP\t%INFO/SUPP_VEC\n' "$VCF" 2>/dev/null > "${CASE_TMP}/supp.tsv"
else
  awk '{print ".\t."}' "${CASE_TMP}/pos.tsv" > "${CASE_TMP}/supp.tsv"
fi
paste "${CASE_TMP}/pos.tsv" "${CASE_TMP}/supp.tsv" | tee "${CASE_TMP}/consensus.tsv"
SAME=$(awk -F'\t' '$2 >= 10100999 && $2 <= 10101001' "${CASE_TMP}/consensus.tsv")
check_eq "one record for the deletion called at 10,100,999 and 10,101,001" "$(grep -c . <<<"$SAME" || true)" 1
check_eq "its support (SUPP)" "$(cut -f5 <<<"$SAME")" 2
check_eq "its callers (SUPP_VEC: manta, delly, cnvpytor)" "$(cut -f6 <<<"$SAME")" 110
check_eq "records for the 600 bp and 5 kb deletions of one 1 kb window" \
  "$(awk -F'\t' '$2 >= 10200000 && $2 <= 10206000' "${CASE_TMP}/consensus.tsv" | grep -c . || true)" 0
check_eq "records for the call only CNVpytor made" \
  "$(awk -F'\t' '$2 >= 10299000 && $2 <= 10311000' "${CASE_TMP}/consensus.tsv" | grep -c . || true)" 0
check_eq "consensus records in all" "$(vcf_count "$VCF")" 1

# GRIDSS reports breakends; a copy of the Manta calls stands in for its file.
bgz manta.vcf.in "${T}/sv_gridss/${T}_gridss.vcf.gz"
run_step 22-survivor-merge.sh "$T"
check_step_exit 22-survivor-merge.sh
check "the log says GRIDSS is left out" has 'gridss: left out' "$(cat "$STEP_LOG")"
check_eq "consensus records with a GRIDSS VCF beside the others" "$(vcf_count "$VCF")" 1
# The GRIDSS copy matches the Manta call, so had it been merged the support would change.
WITH=$(bcf query -f '%POS\t%INFO/SUPP\t%INFO/SUPP_VEC\n' "$VCF" 2>/dev/null | awk -F'\t' '$1 >= 10100999 && $1 <= 10101001')
check_eq "support with a GRIDSS VCF beside the others" "$(cut -f2 <<<"$WITH")" 2
check_eq "callers with a GRIDSS VCF beside the others" "$(cut -f3 <<<"$WITH")" 110

finish
