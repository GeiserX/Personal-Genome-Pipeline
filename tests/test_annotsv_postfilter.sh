#!/usr/bin/env bash
# The class 4 and 5 checks with the real tools: step 05 (AnnotSV, then
# bin/annotsv_postfilter.py in ANNOTSV_IMAGE) and step 22 (SURVIVOR_IMAGE,
# BCFTOOLS_IMAGE) on the made-up data of tests/fixtures/annotsv/class45.sh,
# then both reports from bin/render_report.py.
#
# Checks:
#   1. step 05 writes <sample>_sv_class4_5_checked.tsv from what AnnotSV
#      wrote: EXON (exon 1 end to the base before exon 3 of NM_900001) is
#      class 4 and on exon edges; NEAR, 10 bp inside the intron, is not; EXON
#      carries the benign-set AF 0.05 of the common deletion and is common;
#      every row names its caller (manta); DEL2 and DEL3 are not on exon edges;
#   2. step 22's consensus holds DEL2 (Manta and Delly) and nothing else;
#   3. the reports: DEL2 has 2 callers (manta, delly), DEL3 a single caller,
#      EXON its exon-edge flag and frequency;
#   4. red control: with the consensus emptied (its header only) DEL2 is a
#      single caller in the report, so check 3 can fail.
# Needs docker.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=../versions.env
. "${REPO}/versions.env"
# shellcheck source=fixtures/annotsv/class45.sh
. "${REPO}/tests/fixtures/annotsv/class45.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }
# expect DESC GOT WANT: GOT matches the extended regex WANT.
expect() {
  if [[ "$2" =~ $3 ]]; then pass "$1 ($2)"; else fail "$1: got '$2', want /$3/"; fi
}

command -v docker >/dev/null 2>&1 || { echo "FAIL: docker is needed to run this test"; exit 1; }

G="${WORK}/genome"
S="${G}/s1"
class45_bundle "${G}/annotsv_annotations"
class45_calls "$S" s1

# col TSV ID COLUMN: COLUMN of the checked row of the input record ID.
col() {
  awk -F'\t' -v id="$2" -v c="$3" 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    $h["ID"] == id { print (c in h) ? $h[c] : "no column " c; exit }' "$1"
}

# --- 1. step 05 ---------------------------------------------------------------------------
rc=0
GENOME_DIR="$G" "${REPO}/scripts/05-annotsv.sh" s1 > "${WORK}/step05.log" 2>&1 || rc=$?
expect "step 05 exits 0" "$rc" '^0$'
C="${S}/annotsv/s1_sv_class4_5_checked.tsv"
if [ -s "$C" ]; then
  cat "$C"
  expect "EXON is class 4" "$(col "$C" EXON ACMG_class)" '^4$'
  expect "EXON is on exon edges" "$(col "$C" EXON Exon_edge)" '^yes$'
  expect "EXON names the transcript and the exons" "$(col "$C" EXON Exon_edge_detail)" '^NM_900001 \(FIXG1\): exon 1 and exon 3 edges$'
  expect "EXON: the benign-set AF and its source" "$(col "$C" EXON Population_AF) $(col "$C" EXON Population_source)" '^0\.05 fixture_common_loss'
  expect "EXON is common" "$(col "$C" EXON Common)" '^yes$'
  expect "NEAR, 10 bp inside the intron, is not on exon edges" "$(col "$C" NEAR Exon_edge)" '^no$'
  expect "DEL2 and DEL3 are class 4 and not on exon edges" \
    "$(col "$C" DEL2 ACMG_class)$(col "$C" DEL2 Exon_edge) $(col "$C" DEL3 ACMG_class)$(col "$C" DEL3 Exon_edge)" '^4no 4no$'
  expect "every row names its caller" "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next } { print $h["Caller"] }' "$C" | sort -u | tr '\n' ' ')" '^manta $'
  expect "the step's log counts the checks" "$(grep -o 'Class 4-5 SVs checked: [0-9]*; deletions on exon edges: [0-9]*' "${WORK}/step05.log")" 'checked: 4; deletions on exon edges: 1$'
else
  fail "step 05 wrote no ${C}"
  tail -n 40 "${WORK}/step05.log"
fi

# --- 2. step 22 ---------------------------------------------------------------------------
rc=0
GENOME_DIR="$G" "${REPO}/scripts/22-survivor-merge.sh" s1 > "${WORK}/step22.log" 2>&1 || rc=$?
expect "step 22 exits 0" "$rc" '^0$'
CONS="${S}/sv_merged/s1_sv_consensus.vcf.gz"
expect "the consensus: DEL2 only, two callers" \
  "$(gzip -cd "$CONS" 2>/dev/null | grep -v '^#' | awk -F'\t' '{ s = $8; sub(/.*SUPP=/, "", s); sub(/;.*/, "", s); print $2 ":" s }' | tr '\n' ' ')" \
  '^200000[0-9]:2 $'
[ "$rc" -eq 0 ] || tail -n 30 "${WORK}/step22.log"

# --- 3. the reports -----------------------------------------------------------------------
report() {
  python3 "${REPO}/bin/render_report.py" --sample s1 --sample-dir "$S" --json "${WORK}/summary.json" \
    -o "${WORK}/report.txt" -o "${WORK}/report.html" > /dev/null
}
report
TXT=$(cat "${WORK}/report.txt")
HTML=$(cat "${WORK}/report.html")
sed -n '/Class 4 and 5 (AnnotSV)/,/^$/p' "${WORK}/report.txt"
expect "text: DEL2 has two callers" "$(grep -o 'DEL  20:200000[0-9]-2450000 .*' <<< "$TXT")" 'support: 2 callers \(manta, delly\)$'
expect "text: DEL3 is a single caller" "$(grep -o 'DEL  20:300000[0-9]-3450000 .*' <<< "$TXT")" 'support: single caller \(manta\)$'
expect "text: EXON's flags" "$(grep -c 'on exon edges: NM_900001 (FIXG1): exon 1 and exon 3 edges' <<< "$TXT")" '^1$'
expect "text: the counts" "$(grep -o 'Class 4-5 SVs: .*' <<< "$TXT")" '^Class 4-5 SVs: 4 \(single caller: 3; deletions on exon edges: 1; common: 2\)$'
expect "HTML: DEL2's callers" "$(grep -c '2 callers (manta, delly)' <<< "$HTML")" '^[1-9]'

# --- 4. red control -----------------------------------------------------------------------
gzip -cd "$CONS" | grep '^#' | gzip -c > "${CONS}.tmp" && mv "${CONS}.tmp" "$CONS"
report
expect "red control: with an empty consensus DEL2 is a single caller" \
  "$(grep -o 'DEL  20:200000[0-9]-2450000 .*' "${WORK}/report.txt")" 'support: single caller \(manta\)$'

echo
if [ "$FAILS" -eq 0 ]; then echo "RESULT: ALL PASS"; else
  echo "RESULT: ${FAILS} FAILED"
  for f in step05 step22; do echo "--- ${f}.log"; tail -n 25 "${WORK}/${f}.log" 2>/dev/null; done
  exit 1
fi
