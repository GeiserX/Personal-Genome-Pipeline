#!/usr/bin/env bash
# The interval rows of bin/cn_loci.tsv on real reads: HG002 (male) from the
# 35x BWA-MEM GRCh38 no-ALT BAM of DeepVariant's case study, at the segments
# of four rows (AZFa, AZFb, AZFc_DAZ and 22q11.2 A-B, targets and controls).
#   - the interval block of step 35 (cut from scripts/35-paralogs.sh, run
#     with scripts/lib/common.sh's run_in) on the sample read as male: AZFa
#     and AZFb 1 copy of 1, 22q11.2 A-B 2 of 2, no change; every other row
#     has no read and is not assessed, never a call;
#   - the interval block of modules/local/parascopy/main.nf, cut from the
#     file, in PARASCOPY_IMAGE on the same BAM: the same table, byte for byte;
#   - the same reads read as female: the three AZF rows not applicable,
#     22q11.2 A-B still 2 copies.
# Control: the male table and the female table differ.
. "$(dirname "$0")/lib.sh"

SRC="https://storage.googleapis.com/deepvariant/case-study-testdata/HG002.novaseq.pcr-free.35x.dedup.grch38_no_alt.bam"
P=HG002iv
W="${GENOME_DIR}/${P}"
rm -rf "$W"; mkdir -p "${W}/aligned" "${W}/paralogs"
ROWS=(AZFa AZFb AZFc_DAZ 22q11.2_DGS_proximal_A-B)
args=(); for r in "${ROWS[@]}"; do args+=(--locus "$r"); done
python3 "${REPO}/bin/copy_number.py" bed "${args[@]}" | sort -k1,1 -k2,2n -u > "${W}/fetch.bed"
check_ge "segments of the four rows to fetch" "$(wc -l < "${W}/fetch.bed" | tr -d ' ')" 10

retry() { local i; for i in 1 2 3 4; do "$@" && return 0; [ "$i" -lt 4 ] && { echo "  try ${i} of 4 failed; again in $((i * 15)) s"; sleep $((i * 15)); }; done; return 1; }
fetch_reads() {
  rm -f "${W}/win.bam" "${W}/$(basename "$SRC").bai"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp -v /etc/ssl/certs:/etc/ssl/certs:ro \
    -e CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt -v "${GENOME_DIR}:/genome" -w "/genome/${P}" "$SAMTOOLS_IMAGE" \
    samtools view -b -M -L "/genome/${P}/fetch.bed" -o "/genome/${P}/win.bam" "$SRC"
}
check "HG002's reads at the four rows download" retry fetch_reads
rm -f "${W}/$(basename "$SRC").bai"
sam sort -o "${P}/aligned/${P}_sorted.bam" "${P}/win.bam"
sam index "${P}/aligned/${P}_sorted.bam"
rm -f "${W}/win.bam"
check_ge "HG002 reads over the four rows" "$(sam view -c "${P}/aligned/${P}_sorted.bam" 2>/dev/null || echo 0)" 50000

col() { awk -F'\t' -v l="$2" -v k="$3" 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} $1 == l {print $h[k]}' "$1" 2>/dev/null; }

# --- step 35's interval block, with common.sh's run_in --------------------------------
BLK="${CASE_TMP}/step_block.sh"
sed -n '/^INTERVAL_OUT=/,/^rm -rf "\$INTERVAL_WORK"$/p' "${REPO}/scripts/35-paralogs.sh" > "${CASE_TMP}/step_body.sh"
check_ge "step 35's interval block was cut from the script" "$(wc -l < "${CASE_TMP}/step_body.sh" | tr -d ' ')" 15
step_block() {  # step_block SEX: the block on HG002iv, its table in paralogs/
  { echo 'set -euo pipefail'
    printf 'SAMPLE=%q GENOME_DIR=%q\n' "$P" "$GENOME_DIR"
    printf '. %q\n' "${REPO}/scripts/lib/common.sh"
    printf 'OUTDIR=%q BAM=%q SEX=%q\n' "${W}/paralogs" "${W}/aligned/${P}_sorted.bam" "$1"
    cat "${CASE_TMP}/step_body.sh"; } > "$BLK"
  echo "+ step 35's interval block (sex $1)"
  bash "$BLK" 2>&1 | tail -n 5
}
step_block male
T="${W}/paralogs/${P}_cn_intervals.tsv"
cp "$T" "${CASE_TMP}/male.tsv"
awk -F'\t' '$1 ~ /^(AZF|22q11.2_DGS_proximal_A-B)/' "$T"
check_ge "the table lists every interval row" "$(($(wc -l < "$T") - 1))" 88
for r in AZFa AZFb; do
  check_eq "male: ${r} called" "$(col "$T" "$r" status)" called
  check_eq "male: ${r} 1 copy of 1, no change" "$(col "$T" "$r" copies)/$(col "$T" "$r" expected)/$(col "$T" "$r" change)" \
    "1/1/no change"
done
check_eq "male: 22q11.2 A-B 2 copies of 2, no change" \
  "$(col "$T" 22q11.2_DGS_proximal_A-B copies)/$(col "$T" 22q11.2_DGS_proximal_A-B expected)/$(col "$T" 22q11.2_DGS_proximal_A-B change)" \
  "2/2/no change"
check_eq "rows without a read are not assessed, never called" \
  "$(awk -F'\t' 'NR > 1 && $1 !~ /^(AZF|22q11.2_DGS_proximal_A-B$)/ && $6 == "called"' "$T" | wc -l | tr -d ' ')" 0
echo "cn-intervals: HG002 AZFa $(col "$T" AZFa estimate), AZFb $(col "$T" AZFb estimate)," \
  "AZFc DAZ $(col "$T" AZFc_DAZ estimate) ($(col "$T" AZFc_DAZ status)), 22q11.2 A-B $(col "$T" 22q11.2_DGS_proximal_A-B estimate)" >> "$E2E_NOTES"

# --- the PARASCOPY module's interval block, in its image ------------------------------
MOD="${REPO}/modules/local/parascopy/main.nf"
sed -n '/# Interval rows of bin\/cn_loci.tsv (as scripts\/35-paralogs.sh)/,/rm -f interval_regions.bed interval_bedcov.tsv/p' "$MOD" \
  | sed 's/\\\\/\\/g' > "${CASE_TMP}/module_body.sh"
check_ge "the module's interval block was cut from main.nf" "$(wc -l < "${CASE_TMP}/module_body.sh" | tr -d ' ')" 8
M="${W}/module"
rm -rf "$M"; mkdir -p "$M"
sed -e "s|\${bam}|/genome/${P}/aligned/${P}_sorted.bam|g" -e "s|\${meta.id}|${P}|g" -e "s|\${sex}|male|g" \
  "${CASE_TMP}/module_body.sh" > "${M}/block.sh"
docker run --rm --user "$(id -u):$(id -g)" -v "${GENOME_DIR}:/genome" -v "${REPO}/bin:/pgp-bin:ro" -w "/genome/${P}/module" \
  "$PARASCOPY_IMAGE" sh -c 'PATH="/pgp-bin:$PATH"; sh -e block.sh' 2>&1 | tail -n 3
check "the module's block writes step 35's table, byte for byte" cmp -s "$T" "${M}/${P}_cn_intervals.tsv"

# --- the same reads read as female ----------------------------------------------------
step_block female
for r in AZFa AZFb AZFc_DAZ; do
  check_eq "female: ${r} not applicable" "$(col "$T" "$r" status)" "not applicable"
done
check_eq "female: 22q11.2 A-B still 2 copies" "$(col "$T" 22q11.2_DGS_proximal_A-B copies)" 2
check "control: the male and the female tables differ" bash -c '! cmp -s "$1" "$2"' _ "$T" "${CASE_TMP}/male.tsv"

rm -rf "$W"
finish
