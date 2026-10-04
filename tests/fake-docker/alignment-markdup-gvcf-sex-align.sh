#!/usr/bin/env bash
# Steps 02 and 02a: the minimap2 index uses the sr preset and is named after
# the reference; reads go through fixmate and markdup; THREADS reaches the
# aligner and samtools; the BAM, the minimap2 index and the BWA-MEM2 index
# reach their final names only when complete. An aligner that dies (exit
# 137), a BAM that fails quickcheck and an index build that dies each leave
# no file under the final name, and step 02a says how much memory the
# BWA-MEM2 index needs when the build is killed.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
mkdir -p "${GENOME_DIR}/sample1/fastq"
for r in R1 R2; do printf 'placeholder\n' | gzip -c > "${GENOME_DIR}/sample1/fastq/sample1_${r}.fastq.gz"; done

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
  *" minimap2 "*" -d "*)
    [ "${FAKE_INDEX:-ok}" = ok ] || exit 137
    : > "$(host_path "$(word_after -d)")" ;;
  *" bwa-mem2 index "*)
    [ "${FAKE_INDEX:-ok}" = ok ] || exit 137
    p=$(host_path "$(word_after -p)")
    for e in 0123 amb ann pac bwt.2bit.64; do : > "${p}.${e}"; done ;;
  *" minimap2 "*|*" bwa-mem2 mem "*)
    printf '@HD\tVN:1.6\tSO:unsorted\n'
    [ "${FAKE_ALIGN:-ok}" = ok ] || exit 137 ;;
  *markdup*)
    cat > /dev/null
    printf 'BAM\001 run %s\n' "${FAKE_TAG:-1}" > "$(host_path "$last")" ;;
  *" samtools index "*)
    : > "$(host_path "$last").bai" ;;
  *" samtools quickcheck "*)
    [ "${FAKE_QC:-ok}" = ok ] || exit 1 ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

A="${GENOME_DIR}/sample1/aligned"
REFD="${GENOME_DIR}/reference"
BAM="${A}/sample1_sorted.bam"
MMI="${REFD}/Homo_sapiens_assembly38.sr.mmi"
leftovers() { find "$A" "$REFD" -mindepth 1 \( -name '*tmp*' -o -name '*.part*' \) 2>/dev/null; }

# --- step 02: index, align, mark duplicates ------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 align env THREADS=2 "${SCRIPTS}/02-alignment.sh" sample1
docker_log_has 'minimap2 -x sr -t 2 -d /genome/reference/Homo_sapiens_assembly38\.sr\.mmi\.tmp\.[0-9]+ ' \
  "step 02 did not build the index with -x sr under a temporary name next to the reference"
docker_log_has 'minimap2 -t 2 -a -x sr .*/genome/reference/Homo_sapiens_assembly38\.sr\.mmi ' \
  "step 02 did not map with -t 2 against the .sr.mmi index"
docker_log_has 'samtools.*--cpus 2 .*fixmate.*-m.*sort.*markdup.* _ 2 /genome/sample1/aligned/sample1\.sort_tmp /genome/sample1/aligned/sample1_sorted\.tmp\.bam' \
  "step 02 did not run fixmate -m, sort and markdup on 2 threads into a temporary BAM"
[ -s "$MMI" ] || [ -f "$MMI" ] || fail "no ${MMI} after a successful build"
[ -f "$BAM" ] && [ -f "${BAM}.bai" ] || fail "no BAM and index after a successful run"
[ -z "$(leftovers)" ] || fail "temporary files left after a successful run: $(leftovers)"

# --- the aligner dies: no BAM -------------------------------------------------
rm -f "$BAM" "${BAM}.bai"
run_rc align-dies env FAKE_ALIGN=fail "${SCRIPTS}/02-alignment.sh" sample1
[ "$RC" -ne 0 ] || fail "step 02 exited 0 after the aligner died"
[ ! -e "$BAM" ] || fail "the aligner died and ${BAM} exists"
[ ! -e "${BAM}.bai" ] || fail "the aligner died and ${BAM}.bai exists"
[ -z "$(leftovers)" ] || fail "temporary files left after the aligner died: $(leftovers)"

# --- quickcheck fails: the earlier BAM stays as it was --------------------------
run_expect 0 align-again env FAKE_TAG=first "${SCRIPTS}/02-alignment.sh" sample1
before=$(cat "$BAM")
run_rc align-bad-bam env FAKE_TAG=second FAKE_QC=fail "${SCRIPTS}/02-alignment.sh" sample1
[ "$RC" -ne 0 ] || fail "step 02 exited 0 on a BAM that fails quickcheck"
[ "$(cat "$BAM")" = "$before" ] || fail "a BAM that failed quickcheck replaced the earlier one"
[ -z "$(leftovers)" ] || fail "temporary files left after quickcheck failed: $(leftovers)"

# --- the index build dies: no .sr.mmi -----------------------------------------
rm -f "$MMI"
run_rc index-dies env FAKE_INDEX=fail "${SCRIPTS}/02-alignment.sh" sample1
[ "$RC" -ne 0 ] || fail "step 02 exited 0 after the index build died"
[ ! -e "$MMI" ] || fail "the index build died and ${MMI} exists"
[ -z "$(leftovers)" ] || fail "temporary files left after the index build died: $(leftovers)"

# --- step 02a: BWA-MEM2 ---------------------------------------------------------
A="${GENOME_DIR}/sample1/aligned_bwamem2"
BAM="${A}/sample1_sorted.bam"
IDX="${REFD}/Homo_sapiens_assembly38.fasta"
: > "$FAKE_DOCKER_LOG"
run_expect 0 bwamem2 env THREADS=2 "${SCRIPTS}/02a-alignment-bwamem2.sh" sample1
docker_log_has 'bwa-mem2 index -p /genome/reference/Homo_sapiens_assembly38\.fasta\.tmp\.[0-9]+ ' \
  "step 02a did not build the index under a temporary prefix"
if awk '/bwa-mem2 index/ && /--memory/ { bad = 1 } END { exit !bad }' "$FAKE_DOCKER_LOG"; then
  fail "step 02a still caps the memory of the index build"
fi
docker_log_has 'bwa-mem2 mem -t 2 ' "step 02a did not align with -t 2"
docker_log_has 'samtools.*--cpus 2 .*fixmate.*markdup.* _ 2 /genome/sample1/aligned_bwamem2/sample1\.sort_tmp /genome/sample1/aligned_bwamem2/sample1_sorted\.tmp\.bam' \
  "step 02a did not pipe into fixmate, sort and markdup on 2 threads"
if awk '/bwa-mem2 mem/ && /\.sam( |$)/ { bad = 1 } END { exit !bad }' "$FAKE_DOCKER_LOG"; then
  fail "step 02a still writes a SAM file"
fi
for e in 0123 amb ann pac bwt.2bit.64; do [ -f "${IDX}.${e}" ] || fail "no ${IDX}.${e} after the index build"; done
[ -f "$BAM" ] && [ -f "${BAM}.bai" ] || fail "no BAM and index from step 02a"
[ -z "$(leftovers)" ] || fail "temporary files left after step 02a: $(leftovers)"

for e in 0123 amb ann pac bwt.2bit.64; do rm -f "${IDX}.${e}"; done
run_rc bwamem2-oom env FAKE_INDEX=fail "${SCRIPTS}/02a-alignment-bwamem2.sh" sample1
[ "$RC" -ne 0 ] || fail "step 02a exited 0 after the index build was killed"
output_has bwamem2-oom 'killed \(exit 137\)'
output_has bwamem2-oom 'about 90 GB'
[ ! -e "${IDX}.bwt.2bit.64" ] || fail "the index build was killed and ${IDX}.bwt.2bit.64 exists"
[ -z "$(leftovers)" ] || fail "temporary files left after the killed index build: $(leftovers)"
