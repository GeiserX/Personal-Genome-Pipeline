#!/usr/bin/env bash
# The default reference has no ALT or HLA contigs, and validate-setup.sh
# holds a run to it, with real samtools reading the BAM step 02 wrote:
#   - the fixture reference (the no-ALT analysis set) passes, and the BAM
#     aligned to it matches its .fai sequence for sequence; its header passes
#     the checks docs/realignment.md gives;
#   - the same FASTA with an ALT contig in its .fai fails (named), and that
#     BAM then fails as aligned to another reference, naming the first
#     sequence that differs; ALLOW_ALT_REFERENCE=true turns the ALT failure
#     into a warning. REF_FASTA is given relative to GENOME_DIR, as the docs
#     write it.
. "$(dirname "$0")/lib.sh"

REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
# shellcheck disable=SC2016  # an awk program
check "the fixture reference has no ALT or HLA contig" \
  awk -F'\t' '$1 ~ /_alt$/ || $1 ~ /^HLA-/ {found = 1} END {exit found}' "${REF}.fai"
check "step 02 left a BAM (case 20)" test -s "$BAM"
N=$(grep -c . "${REF}.fai")

# The header checks docs/realignment.md gives for a realigned BAM.
HDR=$(sam view -H "${SAMPLE}/aligned/${SAMPLE}_sorted.bam")
check_eq "the BAM lists as many sequences as the .fai" "$(grep -c '^@SQ' <<< "$HDR")" "$N"
check_eq "the BAM lists no ALT or HLA contig" "$(grep -c -E 'SN:(chr[^[:space:]]*_alt|HLA-)' <<< "$HDR" || true)" 0
check "the minimap2 @PG line names the reference's index" \
  has 'GRCh38_no_alt_analysis_set\.sr\.mmi' "$(grep '^@PG' <<< "$HDR")"

# validate NAME [VAR=VALUE...]: validate-setup.sh for the sample, output in
# CASE_TMP/NAME.log. Its exit code is not checked: on the runner it also fails
# on the VEP cache and the images not pulled.
validate() {
  local name=$1
  shift
  echo "+ $* scripts/validate-setup.sh ${SAMPLE}"
  env "$@" "${REPO}/scripts/validate-setup.sh" "$SAMPLE" > "${CASE_TMP}/${name}.log" 2>&1
  echo "+ exit $?"
  grep -E 'ALT|BAM|Reference|reference' "${CASE_TMP}/${name}.log"
}

validate default
OUT=$(cat "${CASE_TMP}/default.log")
check "default: the reference passes the ALT check" has "\[OK\] +Reference has no ALT or HLA contigs \(${N} sequences\)" "$OUT"
check "default: the BAM matches the reference" has "\[OK\] +BAM header matches the reference: the same ${N} sequences in the same order" "$OUT"
check "default: the BAM is coordinate-sorted" has '\[OK\] +BAM is coordinate-sorted' "$OUT"
check "default: the BAM passes samtools quickcheck" has '\[OK\] +BAM passes samtools quickcheck' "$OUT"
check "default: no reference or BAM failure" lacks '\[FAIL\] +(Reference has|this BAM|BAM is sorted|BAM fails)' "$OUT"

# The same FASTA (a hard link, no copy) under another name, with one ALT
# contig added to its .fai. validate-setup.sh reads only the .fai here.
ALT_NAME=with_alt_e2e
ln -f "$REF" "${GENOME_DIR}/reference/${ALT_NAME}.fasta"
{ cat "${REF}.fai"; printf 'chr22_KI270879v1_alt\t304135\t0\t60\t61\n'; } > "${GENOME_DIR}/reference/${ALT_NAME}.fasta.fai"

validate with-alt REF_FASTA="reference/${ALT_NAME}.fasta"
OUT=$(cat "${CASE_TMP}/with-alt.log")
check "with ALT: the ALT contig fails the reference" has '\[FAIL\] +Reference has 1 ALT/HLA contigs \(first: chr22_KI270879v1_alt\)' "$OUT"
check "with ALT: the BAM fails as aligned to another reference" has '\[FAIL\] +this BAM was aligned to a different reference: realign \(docs/realignment\.md\)' "$OUT"
check "with ALT: the first difference is the extra sequence" \
  has "First difference: sequence $((N + 1)) is missing in the BAM and chr22_KI270879v1_alt \(304135 bp\) in the reference\." "$OUT"

validate with-alt-allowed REF_FASTA="reference/${ALT_NAME}.fasta" ALLOW_ALT_REFERENCE=true
OUT=$(cat "${CASE_TMP}/with-alt-allowed.log")
check "allowed: the ALT contig is a warning" has '\[WARN\] +Reference has 1 ALT/HLA contigs \(allowed by ALLOW_ALT_REFERENCE=true\)' "$OUT"
check "allowed: no failure for the ALT contig" lacks '\[FAIL\] +Reference has' "$OUT"

rm -f "${GENOME_DIR}/reference/${ALT_NAME}.fasta" "${GENOME_DIR}/reference/${ALT_NAME}.fasta.fai"
finish
