#!/usr/bin/env bash
# validate-setup.sh, check_bam_masked_reference: a BAM aligned to the default
# reference is refused when REF_FASTA is the masked one (setup.sh
# --masked-reference), and a BAM aligned to the masked one is refused when
# REF_FASTA is the default, for both aligners. Masking keeps every @SQ name and
# length, so only the aligner's @PG command line tells them apart:
#   minimap2 (step 02, the Nextflow run) names <reference base>.sr.mmi,
#   bwa-mem2 (step 02a) names the FASTA itself.
# The function is read out of validate-setup.sh and run on hand-made headers,
# without Docker. Each verdict is counted from the [OK]/[WARN]/[FAIL] lines,
# and the two refusals are the red controls of the two passes beside them.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0

awk '/^check_bam_masked_reference\(\) \{$/ {on = 1} on {print} on && /^}$/ {exit}' \
  "${REPO}/scripts/validate-setup.sh" > "${WORK}/fn.sh"
if grep -q '_maskedGRC' "${WORK}/fn.sh"; then echo "[PASS] check_bam_masked_reference read from validate-setup.sh"
else echo "[FAIL] check_bam_masked_reference not found in validate-setup.sh"; exit 1; fi

STOCK=/genome/reference/GRCh38_no_alt_analysis_set
MASKED=${STOCK}_maskedGRC
sq=$'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr21\tLN:46709983'
mm2() { printf '%s\n@PG\tID:minimap2\tPN:minimap2\tVN:2.28\tCL:minimap2 -t 16 -a -x sr %s.sr.mmi R1.fq.gz R2.fq.gz\n@PG\tID:samtools\tPN:samtools\tPP:minimap2\tCL:samtools sort -o x.bam' "$sq" "$1"; }
bwa() { printf '%s\n@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\tCL:bwa-mem2 mem -t 16 %s.fasta R1.fq.gz R2.fq.gz' "$sq" "$1"; }
none() { printf '%s\n@PG\tID:DRAGEN\tPN:DRAGEN\tCL:dragen -r /ref' "$sq"; }

# verdict REF HEADER: OK, WARN, FAIL or none, from the function's output.
verdict() {
  (
    # The function sourced below calls these three.
    # shellcheck disable=SC2329
    pass() { echo "OK $1"; }
    # shellcheck disable=SC2329
    warn() { echo "WARN $1"; }
    # shellcheck disable=SC2329
    fail() { echo "FAIL $1"; }
    # shellcheck disable=SC1091
    . "${WORK}/fn.sh"
    REF_FASTA=$1 check_bam_masked_reference "$2"
  ) | awk 'NR == 1 {print $1; n++} END {if (!n) print "none"}'
}

expect() {  # expect DESC REF HEADER WANT
  local got
  got=$(verdict "$2" "$3")
  if [ "$got" = "$4" ]; then echo "[PASS] $1: $got"
  else echo "[FAIL] $1: got $got, want $4"; FAILS=$((FAILS + 1)); fi
}

expect "masked REF_FASTA, minimap2 BAM on the masked index" "${MASKED}.fasta" "$(mm2 "$MASKED")" OK
expect "masked REF_FASTA, minimap2 BAM on the default index (refused)" "${MASKED}.fasta" "$(mm2 "$STOCK")" FAIL
expect "masked REF_FASTA, bwa-mem2 BAM on the masked FASTA" "${MASKED}.fasta" "$(bwa "$MASKED")" OK
expect "masked REF_FASTA, bwa-mem2 BAM on the default FASTA (refused)" "${MASKED}.fasta" "$(bwa "$STOCK")" FAIL
expect "default REF_FASTA, minimap2 BAM on the masked index (refused)" "${STOCK}.fasta" "$(mm2 "$MASKED")" FAIL
expect "default REF_FASTA, bwa-mem2 BAM on the masked FASTA (refused)" "${STOCK}.fasta" "$(bwa "$MASKED")" FAIL
expect "default REF_FASTA, minimap2 BAM on the default index: nothing to say" "${STOCK}.fasta" "$(mm2 "$STOCK")" none
expect "masked REF_FASTA, no aligner @PG line: a warning" "${MASKED}.fasta" "$(none)" WARN
expect "default REF_FASTA, no aligner @PG line: nothing to say" "${STOCK}.fasta" "$(none)" none
expect "masked REF_FASTA, minimap2 @PG renamed minimap2.1 by a merge" "${MASKED}.fasta" \
  "$(mm2 "$STOCK" | sed 's/ID:minimap2/ID:minimap2.1/; s/PN:minimap2/PN:mm/')" FAIL

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${FAILS} FAILED"; fi
[ "$FAILS" -eq 0 ]
