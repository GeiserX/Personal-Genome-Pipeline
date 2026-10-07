#!/usr/bin/env bash
# Step 02 marks duplicates, builds its minimap2 index with the sr preset under
# a name taken from the reference, and leaves no BAM behind when the BAM
# writer or the aligner is killed half way. Reads the BAM case 20 aligned.
. "$(dirname "$0")/lib.sh"

BAM="${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
check "BAM passes samtools quickcheck" sam quickcheck -v "$BAM"
FS=$(sam flagstat "$BAM" 2>/dev/null)
echo "$FS"
DUPS=$(awk '/ duplicates$/ {print $1; exit}' <<< "$FS")
check_ge "reads flagged as duplicates (samtools flagstat)" "${DUPS:-0}" 1
check "the BAM header records samtools markdup" has '^@PG.*ID:samtools.*markdup' "$(sam view -H "$BAM" 2>/dev/null)"
check "the minimap2 index is named after the reference (.sr.mmi)" \
  nonempty reference/GRCh38_no_alt_analysis_set.sr.mmi
check_eq "temporary index files left in reference/" \
  "$(find "${GENOME_DIR}/reference" -maxdepth 1 -name '*.mmi.tmp*' | wc -l | tr -d ' ')" 0

# kill_run NAME COPIES WHEN IMAGE: run step 02 on a copy of the reads (the
# FASTQ COPIES times over) as sample NAME, and kill the container of IMAGE
# once WHEN is true (polled every 0.2 s). Sets KILLED (the container id) and
# KRC (the step's exit code), and checks that no BAM, index or temporary
# file is left behind.
kill_run() {
  local k=$1 copies=$2 when=$3 image=$4 r src i cid log="${CASE_TMP}/kill-$1.log"
  mkdir -p "${GENOME_DIR}/${k}/fastq"
  for r in R1 R2; do
    src="${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_${r}.fastq.gz"
    for ((i = 0; i < copies; i++)); do cat "$src"; done > "${GENOME_DIR}/${k}/fastq/${k}_${r}.fastq.gz"
  done
  "${REPO}/scripts/02-alignment.sh" "$k" > "$log" 2>&1 &
  local pid=$!
  KILLED=""
  for _ in $(seq 1 6000); do
    if eval "$when"; then
      cid=$(docker ps -q --filter "ancestor=${image}" | head -n 1)
      [ -n "$cid" ] && KILLED=$(docker kill "$cid" 2>/dev/null)
      break
    fi
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.2
  done
  KRC=0
  wait "$pid" || KRC=$?
  cat "$log"
  check "${k}: the container was killed mid-run" test -n "$KILLED"
  check "${k}: step 02 exits non-zero after the kill (exit ${KRC})" test "$KRC" -ne 0
  check_eq "${k}: BAM left behind after the kill" "$(find "${GENOME_DIR}/${k}/aligned" -name "${k}_sorted.bam" 2>/dev/null | wc -l | tr -d ' ')" 0
  check_eq "${k}: BAM index left behind after the kill" "$(find "${GENOME_DIR}/${k}/aligned" -name "${k}_sorted.bam.bai" 2>/dev/null | wc -l | tr -d ' ')" 0
  check_eq "${k}: temporary files left behind after the kill" "$(find "${GENOME_DIR}/${k}/aligned" -mindepth 1 2>/dev/null | wc -l | tr -d ' ')" 0
  KLOG=$log
  rm -rf "${GENOME_DIR:?}/${k}"
}

# 1. The BAM writer is killed while it writes the BAM: the old script wrote
#    straight to <sample>_sorted.bam and left the half-written file there,
#    under the name run-all.sh skips on.
W_FILE='[ -n "$(find "${GENOME_DIR}/HG002W/aligned" -maxdepth 1 -name "HG002W_sorted*.bam" -size +0 2>/dev/null)" ]'
kill_run HG002W 1 "$W_FILE" "$SAMTOOLS_IMAGE"

# 2. The aligner is killed while it maps its second batch (the reads go in
#    four times; minimap2 writes a batch when it is mapped). samtools sort
#    stops on the cut stream, so no BAM is left before or after this change;
#    the case keeps it so a change that loses that cannot pass.
A_MAPPED='grep -q "worker_pipeline.*mapped" "${CASE_TMP}/kill-HG002A.log" 2>/dev/null && sleep 3'
kill_run HG002A 4 "$A_MAPPED" "$MINIMAP2_IMAGE"
check_eq "HG002A: batches minimap2 finished before the kill" "$(grep -c 'worker_pipeline.*mapped' "$KLOG")" 1

finish
