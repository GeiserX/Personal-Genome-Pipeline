#!/usr/bin/env bash
# Step 02 marks duplicates, builds its minimap2 index with the sr preset under
# a name taken from the reference, and leaves no BAM behind when the aligner
# is killed half way. Reads the BAM case 20 aligned.
. "$(dirname "$0")/lib.sh"

BAM="${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "BAM passes samtools quickcheck" sam quickcheck -v "$BAM"
FS=$(sam flagstat "$BAM" 2>/dev/null)
echo "$FS"
DUPS=$(awk '/ duplicates$/ {print $1; exit}' <<< "$FS")
check_ge "reads flagged as duplicates (samtools flagstat)" "${DUPS:-0}" 1
check "the BAM header records samtools markdup" has '^@PG.*ID:samtools.*markdup' "$(sam view -H "$BAM" 2>/dev/null)"
check "the minimap2 index is named after the reference (.sr.mmi)" \
  nonempty reference/Homo_sapiens_assembly38.sr.mmi
check_eq "temporary index files left in reference/" \
  "$(find "${GENOME_DIR}/reference" -maxdepth 1 -name '*.mmi.tmp*' | wc -l | tr -d ' ')" 0

# Kill the aligner container while it maps, after it has written part of
# its output: no <sample>_sorted.bam may be left for run-all.sh to skip on.
# minimap2 maps 500 Mb of reads per batch and writes a batch when it is
# done, so the reads go in four times (760 Mb, two batches) and the kill
# comes after the first batch is written ("mapped N sequences"), with the
# second still mapping. Killed earlier, minimap2 would have written
# nothing, not even its header, and the old script left no BAM either.
K=HG002K
mkdir -p "${GENOME_DIR}/${K}/fastq"
for r in R1 R2; do
  src="${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz"
  cat "$src" "$src" "$src" "$src" > "${GENOME_DIR}/${K}/fastq/${K}_${r}.fastq.gz"
done
KLOG="${CASE_TMP}/kill.log"
"${REPO}/scripts/02-alignment.sh" "$K" > "$KLOG" 2>&1 &
PID=$!
KILLED=""
for _ in $(seq 1 900); do
  if grep -q 'worker_pipeline.*mapped' "$KLOG" 2>/dev/null; then
    sleep 3
    CID=$(docker ps -q --filter "ancestor=${MINIMAP2_IMAGE}" | head -n 1)
    [ -n "$CID" ] && KILLED=$(docker kill "$CID" 2>/dev/null)
    break
  fi
  kill -0 "$PID" 2>/dev/null || break
  sleep 1
done
KRC=0
wait "$PID" || KRC=$?
cat "$KLOG"
check "the aligner container was killed while mapping its second batch" test -n "$KILLED"
check_eq "batches minimap2 finished before the kill" "$(grep -c 'worker_pipeline.*mapped' "$KLOG")" 1
check "step 02 exits non-zero after the kill (exit ${KRC})" test "$KRC" -ne 0
check_eq "BAM left behind after the kill" "$(find "${GENOME_DIR}/${K}/aligned" -name "${K}_sorted.bam" 2>/dev/null | wc -l | tr -d ' ')" 0
check_eq "BAM index left behind after the kill" "$(find "${GENOME_DIR}/${K}/aligned" -name "${K}_sorted.bam.bai" 2>/dev/null | wc -l | tr -d ' ')" 0
check_eq "temporary files left behind after the kill" "$(find "${GENOME_DIR}/${K}/aligned" -mindepth 1 2>/dev/null | wc -l | tr -d ' ')" 0
rm -rf "${GENOME_DIR:?}/${K}"

finish
