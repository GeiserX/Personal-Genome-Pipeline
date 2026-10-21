#!/usr/bin/env bash
# Published index files are not older than their data files, so samtools and
# htslib do not warn "The index file is older than the data file" on every
# query. publishDir copies a task's files in parallel, and the small index
# used to land before its BAM or VCF; main.nf's completion handler now sets
# the index times after the last copy.
#   1. Every index the direct `nextflow run` published (case
#      nextflow-from-fastq-2-nextflow, publish mode copy) and the launcher's
#      run (case run-all-launcher, hard links) is not older than its data
#      file, and samtools and bcftools print no such warning on the BAM and
#      the VCF. Negative control: a copied pair with its index set an hour
#      older is counted, and samtools does print the warning for it.
#   2. The handler itself, on a short VCF-only run (roh) into a new outdir
#      with two pairs planted under the sample's folder: pair A's index is
#      newer than the run's start but older than its data (what a parallel
#      copy leaves), pair B's index is from two days before the run. Both
#      fail the check before the run (the red state). After it, A's index is
#      not older than its data, B's is left as it was, and the run's warning
#      names B and not A.
. "$(dirname "$0")/lib.sh"

command -v nextflow >/dev/null || { fail "nextflow is not on PATH"; finish; }
P=HG002P
DIRECT="${GENOME_DIR}/nf-fastq/${P}"
LAUNCHER="${E2E_WORK}/genome-runall/${P}"

# older_indexes DIR: each .bai, .crai, .csi and .tbi under DIR (not its nextflow/
# launch directory) whose data file beside it is newer, one per line
older_indexes() {
  find "$1" -path "$1/nextflow" -prune -o -type f \( -name '*.bai' -o -name '*.crai' -o -name '*.csi' -o -name '*.tbi' \) -print 2>/dev/null \
    | while IFS= read -r idx; do [ -e "${idx%.*}" ] && [ "$idx" -ot "${idx%.*}" ] && echo "$idx"; done
}
# index_pairs DIR: how many indexes under DIR have their data file beside them
index_pairs() {
  find "$1" -path "$1/nextflow" -prune -o -type f \( -name '*.bai' -o -name '*.crai' -o -name '*.csi' -o -name '*.tbi' \) -print 2>/dev/null \
    | while IFS= read -r idx; do [ -e "${idx%.*}" ] && echo "$idx"; done | wc -l | tr -d ' '
}
count() { grep -c . <<< "$1" || true; }

# --- 1. the published pairs of both earlier runs ----------------------------------------
for d in "$DIRECT" "$LAUNCHER"; do
  check "${d#"${E2E_WORK}"/} exists" test -d "$d"
  check_ge "index files with their data file under ${d##*/}/ (${d#"${E2E_WORK}"/})" "$(index_pairs "$d")" 3
  OLDER=$(older_indexes "$d")
  [ -z "$OLDER" ] || printf 'older than their data file:\n%s\n' "$OLDER"
  check_eq "published index files older than their data file (${d#"${E2E_WORK}"/})" "$(count "$OLDER")" 0
done
BAM="nf-fastq/${P}/aligned/${P}_sorted.bam" VCF="nf-fastq/${P}/vcf/${P}.vcf.gz"
WARN='index file is older than the data file'
REGION=$(awk 'NR == 1 {printf "%s:%d-%d", $1, $2 + 1, $2 + 1000}' "${FIXTURE_DIR}/regions.bed")
check "samtools reads the published BAM's index without the warning" lacks "$WARN" "$(sam view -c "$BAM" "$REGION" 2>&1)"
check "bcftools reads the published VCF's index without the warning" lacks "$WARN" "$(bcf view -H -r "$REGION" "$VCF" 2>&1 | head -n 5)"
RED="${GENOME_DIR}/nf-index-times-red"
rm -rf "$RED"; mkdir -p "$RED"
ln "${GENOME_DIR}/${BAM}" "$RED/"; cp -p "${GENOME_DIR}/${BAM}.bai" "$RED/"   # the index is a copy: only it changes
touch -d "@$(( $(stat -c %Y "${RED}/${P}_sorted.bam") - 3600 ))" "${RED}/${P}_sorted.bam.bai"
check_eq "negative control: a pair with its index an hour older is counted" "$(count "$(older_indexes "$RED")")" 1
check "negative control: samtools warns about that index" has "$WARN" \
  "$(sam view -c "nf-index-times-red/${P}_sorted.bam" "$REGION" 2>&1)"
rm -rf "$RED"

# --- 2. the completion handler on a planted pair ----------------------------------------
OUT="${GENOME_DIR}/nf-index-times"
PL="${OUT}/${SAMPLE}/planted"
rm -rf "$OUT"; mkdir -p "$PL" "${CASE_TMP}/run"
NOW=$(date +%s)
for x in a b; do printf 'data\n' > "${PL}/${x}.vcf.gz"; printf 'index\n' > "${PL}/${x}.vcf.gz.tbi"; done
touch -d "@$((NOW + 7200))" "${PL}/a.vcf.gz" "${PL}/b.vcf.gz"
touch -d "@$((NOW + 3600))" "${PL}/a.vcf.gz.tbi"          # A: after the run's start, before its data
touch -d "@$((NOW - 172800))" "${PL}/b.vcf.gz.tbi"        # B: two days before the run
check_eq "before the run: both planted indexes are older than their data" "$(count "$(older_indexes "$PL")")" 2
B_BEFORE=$(stat -c %Y "${PL}/b.vcf.gz.tbi")
printf 'sample,vcf,vcf_index\n%s,%s,%s\n' "$SAMPLE" "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
  "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi" > "${CASE_TMP}/sheet.csv"
echo "+ nextflow run main.nf (VCF only, --tools roh) --outdir ${OUT}"
(
  cd "${CASE_TMP}/run" && nextflow run "${REPO}/main.nf" -profile docker -ansi-log false \
    -work-dir "${E2E_WORK}/nf-work-index-times" --input "${CASE_TMP}/sheet.csv" \
    --reference "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta" --tools roh \
    --outdir "$OUT" --max_cpus 4 --max_memory 14.GB
) > "$STEP_LOG" 2>&1
NF_RC=$?
grep -vE 'Pulling|Waiting|Verifying|Download complete|Pull complete|Already exists' "$STEP_LOG" | sed 's/^/    | /'
cp "${CASE_TMP}/run/.nextflow.log" "${E2E_WORK}/logs/${CASE_NAME}.nextflow.log" 2>/dev/null
check_eq "nextflow run exits 0" "$NF_RC" 0
check "the run published ROH's output" test -n "$(find "${OUT}/${SAMPLE}/roh" -type f 2>/dev/null | head -n 1)"
LEFT=$(older_indexes "$PL")
check "pair A (written in the run, before its data) is no longer older" lacks 'a\.vcf\.gz\.tbi$' "$LEFT"
check "pair B (from before the run) is left older than its data" has 'b\.vcf\.gz\.tbi$' "$LEFT"
check_eq "pair B's index time is unchanged" "$(stat -c %Y "${PL}/b.vcf.gz.tbi")" "$B_BEFORE"
check "pair A's index still holds the index, not something else" has '^index$' "$(cat "${PL}/a.vcf.gz.tbi")"
LOG=$(cat "$STEP_LOG")
check "the run warns about the older index it left, naming B" has 'Index files older than their data file, left as they are: .*b\.vcf\.gz\.tbi' "$LOG"
check "the warning does not name A" lacks 'left as they are: .*a\.vcf\.gz\.tbi' "$LOG"
check "no handler error" lacks 'Could not check the times of the published index files|Failed to invoke .workflow.onComplete' "$LOG"
check_eq "index files older than their data under the run's outputs, the planted pairs aside" \
  "$(count "$(older_indexes "${OUT}/${SAMPLE}" | grep -v '/planted/')")" 0

finish
