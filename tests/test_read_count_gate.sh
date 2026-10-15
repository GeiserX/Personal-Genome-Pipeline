#!/usr/bin/env bash
# The alignment read-count gate in scripts/lib/common.sh: FASTQ_RECORDS_AWK
# (records of a FASTQ), fastp_total_reads (fastp's summary count) and
# read_count_check (the verdict steps 02, 02a and validate-setup.sh act on).
# It checks each verdict and exit status, the table the reports read, a
# shortfall, a surplus, READ_COUNT_CHECK=warn, unusable counts and counts
# above 2^31. It first proves it can fail: a read_count_check that ignores
# the counts must be caught.
#
# READ_COUNT_FROM=FILE takes the helpers from another file, for example the
# common.sh of a commit before the gate.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FROM=${READ_COUNT_FROM:-${REPO}/scripts/lib/common.sh}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# The helpers alone, so the test needs neither versions.env nor GENOME_DIR.
{
  sed -n '/^FASTQ_RECORDS_AWK=/p' "$FROM"
  for f in fastp_total_reads read_count_mode read_count_check; do
    sed -n "/^${f}() {/,/^}/p" "$FROM"
  done
} > "${WORK}/helpers.sh"
for f in FASTQ_RECORDS_AWK= 'fastp_total_reads() {' 'read_count_mode() {' 'read_count_check() {'; do
  grep -qF -- "$f" "${WORK}/helpers.sh" || { echo "FAIL: no ${f%% *} in ${FROM}"; exit 1; }
done
# shellcheck source=/dev/null
. "${WORK}/helpers.sh"
: "${FASTQ_RECORDS_AWK:?}"

BAD=0
ok()  { echo "ok   $*"; }
bad() { echo "FAIL $*"; BAD=$((BAD + 1)); }
expect_eq() { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1: got '$2', want '$3'"; fi; }

fastq() { local i; for ((i = 0; i < $1; i++)); do printf '@r%d\nACGT\n+\nIIII\n' "$i"; done; }

# gate_cases: every case of read_count_check, run against the function now
# defined. Returns the number of failed cases in BAD.
gate_cases() {
  local rc out t="${WORK}/table.tsv"
  BAD=0
  rc=0; out=$(read_count_check 10 "R1 5 + R2 5" "$t" 5 5 2>&1) || rc=$?
  expect_eq "equal counts return 0" "$rc" 0
  expect_eq "equal counts: the table says ok" "$(awk -F'\t' '$1 == "read_count_check" {print $2}' "$t" 2>/dev/null)" ok
  expect_eq "the table's input_reads is the sum of the counts" "$(awk -F'\t' '$1 == "input_reads" {print $2}' "$t" 2>/dev/null)" 10
  expect_eq "the table's header" "$(head -n 1 "$t" 2>/dev/null)" "$(printf 'key\tvalue')"
  rc=0; out=$(read_count_check 6 src "$t" 5 5 2>&1) || rc=$?
  expect_eq "a shortfall returns 1" "$rc" 1
  case "$out" in *"ERROR: the BAM holds 6 primary records, but 10 reads went in (src): 4 reads (40.0%) are missing."*) ok "a shortfall names what is missing" ;;
    *) bad "a shortfall names what is missing: ${out}" ;; esac
  expect_eq "a shortfall: the table says mismatch" "$(awk -F'\t' '$1 == "read_count_check" {print $2}' "$t" 2>/dev/null)" mismatch
  expect_eq "a shortfall: missing_reads" "$(awk -F'\t' '$1 == "missing_reads" {print $2}' "$t" 2>/dev/null)" 4
  rc=0; out=$(READ_COUNT_CHECK=warn read_count_check 6 src "" 5 5 2>&1) || rc=$?
  expect_eq "a shortfall with READ_COUNT_CHECK=warn returns 3" "$rc" 3
  case "$out" in WARNING:*) ok "READ_COUNT_CHECK=warn says WARNING" ;; *) bad "READ_COUNT_CHECK=warn says WARNING: ${out}" ;; esac
  rc=0; out=$(read_count_check 11 src "" 5 5 2>&1) || rc=$?
  expect_eq "more records than reads returns 1" "$rc" 1
  case "$out" in *"more than the 10 reads"*) ok "a surplus is named" ;; *) bad "a surplus is named: ${out}" ;; esac
  rc=0; read_count_check "" src "" 5 5 2>/dev/null || rc=$?
  expect_eq "an empty BAM count returns 2" "$rc" 2
  rc=0; read_count_check 10 src "" 5 x 2>/dev/null || rc=$?
  expect_eq "a count that is not a number returns 2" "$rc" 2
  rc=0; read_count_check 10 src "" 2>/dev/null || rc=$?
  expect_eq "no input count returns 2" "$rc" 2
  rc=0; READ_COUNT_CHECK=maybe read_count_check 10 src "" 10 2>/dev/null || rc=$?
  expect_eq "a READ_COUNT_CHECK other than fail or warn returns 2" "$rc" 2
  rc=0; read_count_check 0 src "" 0 0 >/dev/null 2>&1 || rc=$?
  expect_eq "zero reads in and zero records return 0" "$rc" 0
  rc=0; read_count_check 3000000000 src "" 1500000000 1500000000 >/dev/null 2>&1 || rc=$?
  expect_eq "counts above 2^31 compare" "$rc" 0
  rc=0; read_count_check 2999999999 src "" 1500000000 1500000000 >/dev/null 2>&1 || rc=$?
  expect_eq "one read short above 2^31 is a shortfall" "$rc" 1
}

# --- control: a check that ignores the counts must be caught ------------------
eval "$(declare -f read_count_check | sed '1s/read_count_check/real_read_count_check/')"
read_count_check() { echo "Read count: fine"; [ -z "${3:-}" ] || printf 'key\tvalue\nread_count_check\tok\ninput_reads\t10\n' > "$3"; return 0; }
gate_cases > "${WORK}/control.txt"
if [ "$BAD" -eq 0 ]; then
  cat "${WORK}/control.txt"
  echo "FAIL: a read_count_check that ignores the counts passed every case"
  exit 1
fi
echo "Control: a read_count_check that ignores the counts fails ${BAD} cases."
eval "$(declare -f real_read_count_check | sed '1s/real_read_count_check/read_count_check/')"

# --- the real helpers ---------------------------------------------------------
gate_cases
FAILED=$BAD

expect_eq "FASTQ_RECORDS_AWK counts records" "$(fastq 3 | awk "$FASTQ_RECORDS_AWK")" 3
expect_eq "FASTQ_RECORDS_AWK on no input" "$(: | awk "$FASTQ_RECORDS_AWK")" 0
rc=0; fastq 2 | head -n 7 | awk "$FASTQ_RECORDS_AWK" >/dev/null 2>&1 || rc=$?
expect_eq "FASTQ_RECORDS_AWK exits 3 on a record cut short" "$rc" 3
FAILED=$((FAILED + BAD))

# fastp's report, laid out as fastp writes it (one key per line, tabs): the
# count is summary.after_filtering's, not before_filtering's and not
# read1_after_filtering's.
cat > "${WORK}/fastp.json" <<'JSON'
{
	"summary": {
		"fastp_version": "1.4.0",
		"sequencing": "paired end (151 cycles + 151 cycles)",
		"before_filtering": {
			"total_reads":2000,
			"total_bases":302000
		},
		"after_filtering": {
			"total_reads":1960,
			"total_bases":290000
		}
	},
	"read1_after_filtering": {
		"total_reads":980
	}
}
JSON
BAD=0
expect_eq "fastp_total_reads reads summary.after_filtering.total_reads" "$(fastp_total_reads "${WORK}/fastp.json")" 1960
tr -d '\n\t ' < "${WORK}/fastp.json" > "${WORK}/oneline.json"
expect_eq "fastp_total_reads on a one-line report" "$(fastp_total_reads "${WORK}/oneline.json")" 1960
printf '{"summary": {"before_filtering": {"total_reads":5}}}\n' > "${WORK}/nofilter.json"
rc=0; fastp_total_reads "${WORK}/nofilter.json" >/dev/null || rc=$?
expect_eq "fastp_total_reads fails without an after_filtering count" "$rc" 1
: > "${WORK}/empty.json"
rc=0; fastp_total_reads "${WORK}/empty.json" >/dev/null || rc=$?
expect_eq "fastp_total_reads fails on an empty file" "$rc" 1
FAILED=$((FAILED + BAD))

if [ "$FAILED" -ne 0 ]; then
  echo "FAIL: ${FAILED} case(s)"
  exit 1
fi
echo "PASS: every read-count case"
