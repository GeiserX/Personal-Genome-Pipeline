#!/usr/bin/env bash
# Steps 02 and 02a count the BAM's primary records (samtools view -c -F 0x900)
# and compare them with the FASTQ records of R1 plus R2. Equal: the BAM and
# qc/<sample>_read_count.tsv (02a: _bwamem2_read_count.tsv) are kept. Fewer
# or more: the step exits 1, keeps no new BAM or table, and leaves an earlier
# BAM and table as they were. READ_COUNT_CHECK=warn keeps the BAM with the
# table saying mismatch; a READ_COUNT_CHECK other than fail or warn stops the
# step before it aligns; a FASTQ cut inside a record stops it too.
# validate-setup.sh runs the same comparison for an existing BAM with fastp's
# report beside it: [OK], [FAIL], or [WARN] under READ_COUNT_CHECK=warn, and
# nothing without the report.
# The fake samtools answers `view -c` with FAKE_PRIMARY (default 6).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
# fastq N: N four-line records
fastq() { local i; for ((i = 0; i < $1; i++)); do printf '@r%d\nACGT\n+\nIIII\n' "$i"; done; }
mkdir -p "${GENOME_DIR}/sample1/fastq"
for r in R1 R2; do fastq 3 | gzip -c > "${GENOME_DIR}/sample1/fastq/sample1_${r}.fastq.gz"; done

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=" $* "
last=${!#}
word_after() {
  local -a w
  read -r -a w <<<"$args"
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]}"; return 0; fi
  done
  return 1
}
case "$args" in
  *" minimap2 "*" -d "*) : > "$(host_path "$(word_after -d)")" ;;
  *" bwa-mem2 index "*)
    p=$(host_path "$(word_after -p)")
    for e in 0123 amb ann pac bwt.2bit.64; do : > "${p}.${e}"; done ;;
  *" minimap2 "*|*" bwa-mem2 mem "*) printf '@HD\tVN:1.6\tSO:unsorted\n' ;;
  *markdup*)
    cat > /dev/null
    out=$(host_path "$last")
    printf 'BAI\001\n' > "${out}.bai"
    printf 'BAM\001 run %s\n' "${FAKE_TAG:-1}" > "$out" ;;
  *" samtools view -c "*) echo "${FAKE_PRIMARY:-6}" ;;
  *" samtools quickcheck "*) : ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

S="${GENOME_DIR}/sample1"
leftovers() { find "${S}/aligned" "${S}/aligned_bwamem2" "${S}/qc" -mindepth 1 \( -name '*tmp*' -o -name '*.count' \) 2>/dev/null; }
val() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$2"; }

# --- step 02: equal counts ----------------------------------------------------------
BAM="${S}/aligned/sample1_sorted.bam"
T="${S}/qc/sample1_read_count.tsv"
: > "$FAKE_DOCKER_LOG"
run_expect 0 equal env THREADS=2 "${SCRIPTS}/02-alignment.sh" sample1
docker_log_has 'samtools view -c -F 0x900 -@ 2 /genome/sample1/aligned/sample1_sorted\.tmp\.bam' \
  "step 02 did not count the primary records of the temporary BAM on 2 threads"
[ -f "$BAM" ] && [ -f "${BAM}.bai" ] || fail "no BAM and index after equal counts"
[ -f "$T" ] || fail "no ${T} after equal counts"
[ "$(head -n 1 "$T")" = "$(printf 'key\tvalue')" ] || fail "${T} has no key/value header"
[ "$(val input_reads "$T")" = 6 ] || fail "input_reads is $(val input_reads "$T"), want 6 (3 records in R1 + 3 in R2)"
[ "$(val primary_records "$T")" = 6 ] || fail "primary_records is $(val primary_records "$T"), want 6"
[ "$(val missing_reads "$T")" = 0 ] || fail "missing_reads is $(val missing_reads "$T"), want 0"
[ "$(val read_count_check "$T")" = ok ] || fail "read_count_check is $(val read_count_check "$T"), want ok"
[ "$(val input_source "$T")" = "fastq/sample1_R1.fastq.gz 3 records + fastq/sample1_R2.fastq.gz 3 records" ] \
  || fail "input_source is '$(val input_source "$T")'"
output_has equal 'Read count: the BAM holds 6 primary records, the 6 reads that went in'
[ -z "$(leftovers)" ] || fail "temporary files left after equal counts: $(leftovers)"

# --- step 02: two reads short. The earlier BAM and table stay as they were -----------
before=$(cat "$BAM"); table_before=$(cat "$T")
run_expect 1 short env FAKE_PRIMARY=4 FAKE_TAG=second "${SCRIPTS}/02-alignment.sh" sample1
output_has short 'ERROR: the BAM holds 4 primary records, but 6 reads went in \(.*\): 2 reads \(33\.3%\) are missing\.'
output_has short 'READ_COUNT_CHECK=warn goes on anyway'
output_has short 'The earlier .*/sample1_sorted\.bam and its read count are left as they were'
[ "$(cat "$BAM")" = "$before" ] || fail "a BAM with reads missing replaced the earlier one"
[ "$(cat "$T")" = "$table_before" ] || fail "a BAM with reads missing replaced the earlier read-count table"
[ -z "$(leftovers)" ] || fail "temporary files left after the shortfall: $(leftovers)"

# --- no earlier BAM: nothing is left at all ---------------------------------------
rm -f "$BAM" "${BAM}.bai" "$T"
run_expect 1 short-fresh env FAKE_PRIMARY=4 "${SCRIPTS}/02-alignment.sh" sample1
[ ! -e "$BAM" ] && [ ! -e "${BAM}.bai" ] || fail "a BAM with reads missing was kept"
[ ! -e "$T" ] || fail "a read-count table was kept for a BAM that was not"
[ -z "$(leftovers)" ] || fail "temporary files left after the shortfall: $(leftovers)"

# --- more records than reads: the same stop ------------------------------------------
run_expect 1 surplus env FAKE_PRIMARY=7 "${SCRIPTS}/02-alignment.sh" sample1
output_has surplus 'more than the 6 reads that went in'
[ ! -e "$BAM" ] || fail "a BAM with more records than reads was kept"

# --- READ_COUNT_CHECK=warn: the BAM is kept, the table says mismatch ----------------
run_expect 0 warn env FAKE_PRIMARY=4 READ_COUNT_CHECK=warn "${SCRIPTS}/02-alignment.sh" sample1
output_has warn 'WARNING: the BAM holds 4 primary records, but 6 reads went in'
output_has warn 'READ_COUNT_CHECK=warn: going on'
[ -f "$BAM" ] && [ -f "${BAM}.bai" ] || fail "READ_COUNT_CHECK=warn did not keep the BAM"
[ "$(val read_count_check "$T")" = mismatch ] || fail "read_count_check is $(val read_count_check "$T"), want mismatch"
[ "$(val missing_reads "$T")" = 2 ] || fail "missing_reads is $(val missing_reads "$T"), want 2"

# --- a bad READ_COUNT_CHECK stops the step before it aligns ---------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 1 bad-mode env READ_COUNT_CHECK=maybe "${SCRIPTS}/02-alignment.sh" sample1
output_has bad-mode "READ_COUNT_CHECK must be 'fail' or 'warn', got 'maybe'"
if grep -q minimap2 "$FAKE_DOCKER_LOG"; then fail "step 02 aligned before it checked READ_COUNT_CHECK"; fi

# --- R2 cut inside a record: the count fails, no BAM ---------------------------------
rm -f "$BAM" "${BAM}.bai" "$T"
fastq 3 | head -n 10 | gzip -c > "${S}/fastq/sample1_R2.fastq.gz"
run_expect 1 cut-record "${SCRIPTS}/02-alignment.sh" sample1
output_has cut-record '10 lines, not a multiple of 4: a record is cut short'
output_has cut-record 'ERROR: could not count the records of'
[ ! -e "$BAM" ] || fail "a BAM was kept although its FASTQ could not be counted"
[ -z "$(leftovers)" ] || fail "temporary files left after the cut record: $(leftovers)"
fastq 3 | gzip -c > "${S}/fastq/sample1_R2.fastq.gz"

# --- step 02a: the same gate ---------------------------------------------------------
BAM="${S}/aligned_bwamem2/sample1_sorted.bam"
T="${S}/qc/sample1_bwamem2_read_count.tsv"
: > "$FAKE_DOCKER_LOG"
run_expect 0 bwa-equal env THREADS=2 "${SCRIPTS}/02a-alignment-bwamem2.sh" sample1
docker_log_has 'samtools view -c -F 0x900 -@ 2 /genome/sample1/aligned_bwamem2/sample1_sorted\.tmp\.bam' \
  "step 02a did not count the primary records of the temporary BAM"
[ -f "$BAM" ] || fail "no BAM from step 02a after equal counts"
[ "$(val read_count_check "$T")" = ok ] || fail "step 02a: read_count_check is '$(val read_count_check "$T")', want ok"
rm -f "$BAM" "${BAM}.bai" "$T"
run_expect 1 bwa-short env FAKE_PRIMARY=5 "${SCRIPTS}/02a-alignment-bwamem2.sh" sample1
output_has bwa-short '1 reads \(16\.7%\) are missing'
[ ! -e "$BAM" ] || fail "step 02a kept a BAM with a read missing"
[ ! -e "$T" ] || fail "step 02a kept a table for a BAM it did not keep"
run_expect 0 bwa-warn env FAKE_PRIMARY=5 READ_COUNT_CHECK=warn "${SCRIPTS}/02a-alignment-bwamem2.sh" sample1
[ "$(val read_count_check "$T")" = mismatch ] || fail "step 02a with warn: read_count_check is '$(val read_count_check "$T")'"
[ -z "$(leftovers)" ] || fail "temporary files left after step 02a: $(leftovers)"

# --- validate-setup.sh: an existing BAM with fastp's report beside it ----------------
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" v1
mkdir -p "${GENOME_DIR}/v1/fastq_trimmed"
cat > "${GENOME_DIR}/v1/fastq_trimmed/v1_fastp.json" <<'JSON'
{
	"summary": {
		"before_filtering": {
			"total_reads":8
		},
		"after_filtering": {
			"total_reads":6
		}
	}
}
JSON
: > "$FAKE_DOCKER_LOG"
run_rc v-equal "${SCRIPTS}/validate-setup.sh" v1
output_has v-equal '\[OK\].*BAM holds the 6 reads fastp wrote as primary records'
docker_log_has 'samtools view -c -F 0x900 -@ [0-9]+ /genome/v1/aligned/v1_sorted\.bam' \
  "validate-setup.sh did not count the BAM's primary records"
run_rc v-short env FAKE_PRIMARY=4 "${SCRIPTS}/validate-setup.sh" v1
expect_rc v-short 1
output_has v-short '\[FAIL\].*BAM read count: 4 primary records, but fastp wrote 6 reads'
output_has v-short '2 reads \(33\.3%\) are missing'
output_has v-short 'FASTQ_SUBDIR=fastq'
run_rc v-warn env FAKE_PRIMARY=4 READ_COUNT_CHECK=warn "${SCRIPTS}/validate-setup.sh" v1
output_has v-warn '\[WARN\].*BAM read count: 4 primary records, but fastp wrote 6 reads \(READ_COUNT_CHECK=warn\)'
output_lacks v-warn '\[FAIL\].*BAM read count'
# Without fastp's report there is nothing to compare with: no line, no count.
rm -f "${GENOME_DIR}/v1/fastq_trimmed/v1_fastp.json"
: > "$FAKE_DOCKER_LOG"
run_rc v-none env FAKE_PRIMARY=4 "${SCRIPTS}/validate-setup.sh" v1
output_lacks v-none 'BAM read count|reads fastp wrote'
if grep -q 'samtools view -c' "$FAKE_DOCKER_LOG"; then fail "validate-setup.sh counted a BAM without fastp's report"; fi
