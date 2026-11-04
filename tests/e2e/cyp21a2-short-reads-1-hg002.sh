#!/usr/bin/env bash
# The CYP21A2 row of bin/cn_loci.tsv (call psv) on real reads: HG002 from the
# 35x BWA-MEM GRCh38 no-ALT BAM of DeepVariant's case study, over the row's
# segments (CYP21A2, CYP21A1P, the two flanks). Three BAMs:
#   - unmodified: 2 CYP21A2 copies, module total 4, no finding (or not
#     assessed, never a deletion or a conversion);
#   - a planted In2G conversion: every second read that carries CYP21A2's
#     base at In2G (chr6:32039081, a read without clips or indels) gets
#     CYP21A1P's G there: a possible conversion at In2G, its ratio under 0.7;
#   - a planted 30 kb chimeric deletion: every second read is dropped that
#     starts on CYP21A2 before 32,039,600 (its In2G and I172N part) or on the
#     matching CYP21A1P stretch after 32,006,864 (the rest), so one
#     chromosome holds a chimera of CYP21A1P's 5' part and CYP21A2's 3' part:
#     module total 3, 1 CYP21A2 copy, finding deletion.
# The block of step 35 runs on each (cut from scripts/35-paralogs.sh, with
# scripts/lib/common.sh's run_in), and the PARASCOPY module's block, cut from
# its main.nf, writes the same table on the unmodified BAM, byte for byte.
. "$(dirname "$0")/lib.sh"

SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
P=HG002cyp
W="${GENOME_DIR}/${P}"
rm -rf "$W"; mkdir -p "${W}/aligned" "${W}/paralogs"
python3 "${REPO}/bin/copy_number.py" bed --locus CYP21A2 | sort -k1,1 -k2,2n > "${W}/fetch.bed"
check_ge "segments of the CYP21A2 row to fetch" "$(wc -l < "${W}/fetch.bed" | tr -d ' ')" 10

retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
fetch_reads() {
  rm -f "${W}/win.bam" "${W}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" "$SAMTOOLS_IMAGE" \
    samtools view -b -M -L "/genome/${P}/fetch.bed" -o "/genome/${P}/win.bam" "$SRC"
}
check "HG002's reads over the CYP21A2 row download" retry fetch_reads
rm -f "${W}/$(basename "$SRC").bai"
sam view -h "${P}/win.bam" > "${W}/orig.sam"
rm -f "${W}/win.bam"

# bam NAME SAM: a sorted, indexed BAM of a SAM file, at ${P}/<NAME>.bam
bam() {
  sam sort -o "${P}/$1.bam" "${P}/$2" && sam index "${P}/$1.bam"
}
bam unmodified orig.sam
check_ge "HG002 reads over the CYP21A2 row" "$(sam view -c "${P}/unmodified.bam" 2>/dev/null || echo 0)" 5000

# The planted In2G conversion.
awk -F'\t' -v OFS='\t' '/^@/ {print; next}
  $3 == "chr6" && $6 ~ /^[0-9]+M$/ && $4 <= 32039081 && 32039081 < $4 + length($10) {
    o = 32039081 - $4 + 1; b = substr($10, o, 1)
    if (b == "C" || b == "A") { n++; if (n % 2) { $10 = substr($10, 1, o - 1) "G" substr($10, o + 1); f++ } }
  }
  {print}
  END {print f " of " n " reads with the CYP21A2 base at In2G changed" > "/dev/stderr"}' \
  "${W}/orig.sam" > "${W}/conversion.sam" 2> "${CASE_TMP}/conversion.log"
cat "${CASE_TMP}/conversion.log"
check_ge "reads changed at In2G" "$(awk '{print $1}' "${CASE_TMP}/conversion.log")" 10
bam conversion conversion.sam

# The planted chimeric deletion.
awk -F'\t' '/^@/ {print; next}
  $3 == "chr6" && (($4 >= 32038800 && $4 < 32039600) || ($4 >= 32006864 && $4 < 32008909)) { if (n++ % 2) next }
  {print}' "${W}/orig.sam" > "${W}/deletion.sam"
bam deletion deletion.sam
rm -f "${W}"/*.sam

col() { awk -F'\t' -v k="$2" '$1 == k {print $2}' "$1" 2>/dev/null; }

# --- step 35's CYP21A2 block, with common.sh's run_in ---------------------------------
sed -n '/^CYP21_OUT=/,/^rm -rf "\$CYP21_WORK"$/p' "${REPO}/scripts/35-paralogs.sh" > "${CASE_TMP}/step_body.sh"
check_ge "step 35's CYP21A2 block was cut from the script" "$(wc -l < "${CASE_TMP}/step_body.sh" | tr -d ' ')" 20
step_block() {  # step_block NAME: the block on ${P}/<NAME>.bam, its table copied to CASE_TMP/<NAME>.tsv
  { echo 'set -euo pipefail'
    printf 'SAMPLE=%q GENOME_DIR=%q\n' "$P" "$GENOME_DIR"
    printf '. %q\n' "${REPO}/scripts/lib/common.sh"
    printf 'OUTDIR=%q BAM=%q\n' "${W}/paralogs" "${W}/$1.bam"
    cat "${CASE_TMP}/step_body.sh"; } > "${CASE_TMP}/block.sh"
  echo "+ step 35's CYP21A2 block ($1)"
  bash "${CASE_TMP}/block.sh" 2>&1 | tail -n 5
  cp "${W}/paralogs/${P}_cyp21a2_psv.tsv" "${CASE_TMP}/$1.tsv" 2>/dev/null || : > "${CASE_TMP}/$1.tsv"
  cat "${CASE_TMP}/$1.tsv"
}
for v in unmodified conversion deletion; do step_block "$v"; done

U="${CASE_TMP}/unmodified.tsv" C="${CASE_TMP}/conversion.tsv" D="${CASE_TMP}/deletion.tsv"
if [ "$(col "$U" status)" = called ]; then
  check_eq "unmodified: 2 CYP21A2 copies, module total 4, no finding" \
    "$(col "$U" cyp21a2_copies)/$(col "$U" module_total)/$(col "$U" finding)" "2/4/none"
else
  check_eq "unmodified: not called, so not assessed (never a finding)" "$(col "$U" status)/$(col "$U" finding)" \
    "not assessed/."
fi
check_eq "planted In2G conversion: called, a possible conversion" "$(col "$C" status)/$(col "$C" finding)" \
  "called/conversion"
check "planted In2G conversion: In2G among the conversion sites" has '(^|,)PSV_In2G(,|$)' "$(col "$C" conversion_sites)"
ratio=$(col "$C" psv_PSV_In2G | sed -n 's/.*ratio \([0-9.]*\)$/\1/p')
check "planted In2G conversion: the In2G ratio (${ratio:-none}) is under 0.7" \
  awk -v r="$ratio" 'BEGIN {exit !(r != "" && r + 0 < 0.7)}'
check_eq "planted chimeric deletion: 1 CYP21A2 copy, module total 3, finding deletion" \
  "$(col "$D" cyp21a2_copies)/$(col "$D" module_total)/$(col "$D" finding)" "1/3/deletion"
check "control: the unmodified and the conversion tables differ" bash -c '! cmp -s "$1" "$2"' _ "$U" "$C"
echo "cyp21a2-short-reads: HG002 $(col "$U" status), $(col "$U" cyp21a2_copies) copies, module total" \
  "$(col "$U" module_total) ($(col "$U" estimate)), finding $(col "$U" finding); In2G $(col "$U" psv_PSV_In2G);" \
  "planted conversion In2G ratio ${ratio:-none}; planted deletion $(col "$D" cyp21a2_copies) copies," \
  "module total $(col "$D" module_total) ($(col "$D" estimate))" >> "$E2E_NOTES"

# --- the PARASCOPY module's CYP21A2 block, in its image -------------------------------
MOD="${REPO}/modules/local/parascopy/main.nf"
sed -n '/# CYP21A2 from depth and PSVs (as scripts\/35-paralogs.sh)/,/rm -f cyp21a2_regions.bed/p' "$MOD" \
  | sed -e 's/\\\\/\\/g' -e 's/\\\$/$/g' > "${CASE_TMP}/module_body.sh"
check_ge "the module's CYP21A2 block was cut from main.nf" "$(wc -l < "${CASE_TMP}/module_body.sh" | tr -d ' ')" 15
M="${W}/module"
rm -rf "$M"; mkdir -p "$M"
sed -e "s|\${bam}|/genome/${P}/unmodified.bam|g" -e "s|\${meta.id}|${P}|g" "${CASE_TMP}/module_body.sh" > "${M}/block.sh"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${P}/module" \
  "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
check "the module's block writes step 35's table, byte for byte" cmp -s "$U" "${M}/${P}_cyp21a2_psv.tsv"

rm -rf "$W"
finish
