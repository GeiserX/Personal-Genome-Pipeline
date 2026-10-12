#!/usr/bin/env bash
# The memory the PRS tasks need on the default score list, against what
# conf/base.config asks for.
#
# Two tasks of the PRS path read nothing but the scoring files of
# assets/pgs_scores.tsv, which are public PGS Catalog downloads, so their
# memory can be measured without a genome:
#   PRS_PREPARE        collect_summary.py prs-format, in PYTHON_IMAGE
#   FORMAT_SCOREFILES  pgsc_calc's first task, pgscatalog-format in
#                      PGSC_UTILS_IMAGE. pgsc_calc gives it 2 threads on its
#                      first attempt and 4 on a retry (its process_medium
#                      label: 2 CPUs times the attempt, capped by --max_cpus),
#                      in a container limited to its --max_memory, which is
#                      what the PRS process asks for.
#
#   bash tests/test_prs_memory.sh
#       No Docker, no network (guard.yml runs it on every pull request):
#       conf/base.config must ask for at least HEADROOM_PCT percent of the
#       peaks recorded below, for PRS and for PRS_PREPARE. A copy of the
#       config with a request under that, and one with no request, must both
#       be reported first.
#
#   bash tests/test_prs_memory.sh --measure DIR
#       The measurement (e2e.yml, job prs-memory, on dispatch; about 300 MB
#       of downloads into DIR). Runs both tasks on the default list in the
#       pinned images, each in a container limited to its request in
#       conf/base.config (capped below the runner's RAM), and reads its memory
#       from the container's cgroup each second: the highest of memory.peak
#       and anonymous memory plus swap. It fails when a task is killed, when
#       that figure times the headroom is over the request, or when it is
#       more than 5% over the one recorded here (memory.peak varies by about
#       0.1 GB between runs; past that, RECORDED_* and the docs need the new
#       figure). A container that allocates 2 GB under a 1 GB limit without
#       swap must be reported as killed first.
#       MEASURE_SWAP_MB lets each container swap that much on top of its
#       limit, so a request above the runner's RAM can still finish.
#       FORMAT_LIMIT_GB and PREPARE_LIMIT_GB replace the requests and
#       FORMAT_THREADS the thread count (default 1, what PRS gives it). With
#       EXPECT=killed every task must be killed instead: FORMAT_LIMIT_GB=8
#       PREPARE_LIMIT_GB=2 FORMAT_THREADS=2 is what the process labels and
#       pgsc_calc gave before the requests existed.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)

# The highest figures any run has given, in MB, on ubuntu-latest (4 CPUs,
# 16 GB; docs/25-prs.md, Memory). PRS_PREPARE: memory.peak
# of its container, 11,614 to 11,734 in ten runs. FORMAT_SCOREFILES on one
# thread: Nextflow's peak_rss in a whole pgsc_calc run (18.1 GB, the sum of
# its processes' peaks); this script's own figure was 14,467 to 14,544 with swap.
RECORDED_PREPARE_MB=11734
RECORDED_FORMAT_MB=18534
HEADROOM_PCT=125

# request PROCESS CONFIG: the GB of "memory = { N.GB * task.attempt }" in
# CONFIG's withName block for PROCESS; fails when there is none.
request() {
  local gb
  gb=$(awk -v p="withName: '$1' {" '
    index($0, p) { inside = 1; next }
    inside && /^[ \t]*}/ { exit }
    inside && match($0, /memory[ \t]*=[ \t]*\{[ \t]*[0-9]+\.GB \* task\.attempt/) {
      s = substr($0, RSTART, RLENGTH); sub(/\.GB.*/, "", s); sub(/.*[^0-9]/, "", s); print s; exit
    }' "$2")
  [ -n "$gb" ] || { echo "no 'memory = { N.GB * task.attempt }' in a withName: '$1' block of $2"; return 1; }
  echo "$gb"
}

# enough PROCESS PEAK_MB CONFIG: 0 when CONFIG's request for PROCESS is at
# least HEADROOM_PCT percent of PEAK_MB.
enough() {
  local gb
  gb=$(request "$1" "$3") || { echo "$gb"; return 1; }
  if [ $((gb * 1024 * 100)) -lt $(($2 * HEADROOM_PCT)) ]; then
    echo "$1 asks for ${gb} GB, under ${HEADROOM_PCT}% of the ${2} MB peak"
    return 1
  fi
  echo "$1 asks for ${gb} GB, at least ${HEADROOM_PCT}% of the ${2} MB peak"
}

if [ "${1:-}" != --measure ]; then
  WORK=$(mktemp -d)
  trap 'rm -rf "$WORK"' EXIT
  CONFIG="${REPO}/conf/base.config"
  FAILS=0
  # Controls: a request under the recorded peak and a missing request must fail.
  sed "/withName: 'PRS' {/,/}/ s/[0-9][0-9]*\.GB/1.GB/" "$CONFIG" > "${WORK}/low.config"
  if out=$(enough PRS "$RECORDED_FORMAT_MB" "${WORK}/low.config"); then
    echo "FAIL: a planted 1 GB request was not reported"; FAILS=$((FAILS + 1))
  else
    echo "ok:   planted low request reported: ${out}"
  fi
  sed "s/withName: 'PRS_PREPARE' {/withName: 'ELSE' {/" "$CONFIG" > "${WORK}/none.config"
  if out=$(enough PRS_PREPARE "$RECORDED_PREPARE_MB" "${WORK}/none.config"); then
    echo "FAIL: a missing request was not reported"; FAILS=$((FAILS + 1))
  else
    echo "ok:   planted missing request reported: ${out}"
  fi
  for pair in "PRS:${RECORDED_FORMAT_MB}" "PRS_PREPARE:${RECORDED_PREPARE_MB}"; do
    if out=$(enough "${pair%%:*}" "${pair#*:}" "$CONFIG"); then
      echo "ok:   ${out}"
    else
      echo "FAIL: ${out}"; FAILS=$((FAILS + 1))
    fi
  done
  [ "$FAILS" -eq 0 ]
  exit
fi

# --- the measurement -------------------------------------------------------------
DIR=${2:?Usage: $0 --measure DIR}
# shellcheck source=../scripts/lib/common.sh
. "${REPO}/scripts/lib/common.sh"
mkdir -p "${DIR}/scores"
DIR=$(cd "$DIR" && pwd)
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}
[ -r /sys/fs/cgroup/cgroup.controllers ] || { echo "ERROR: this needs cgroup v2 (memory.peak)" >&2; exit 1; }
HOST_MB=$(awk '/^MemTotal:/ {print int($2 / 1024)}' /proc/meminfo)
# A container limit above this would let the task take the whole machine down
# with it instead of being killed alone.
CAP_MB=$((HOST_MB - 1536))
FAILS=0

# measured LIMIT_MB DOCKER_RUN_ARGS...: runs the container with that memory
# limit (plus MEASURE_SWAP_MB of swap), and sets M_RC (exit code), M_KILLED
# (true when the kernel killed one of its processes for memory), M_PEAK_MB
# (memory.peak) and M_USED_MB (the highest memory.current minus page cache,
# plus swap, sampled each second).
measured() {
  local limit=$1 cid pid cg p u t0; shift
  M_PEAK_MB=0 M_USED_MB=0
  t0=$(date +%s)
  cid=$(docker run -d --memory "${limit}m" --memory-swap "$((limit + ${MEASURE_SWAP_MB:-0}))m" "$@")
  pid=$(docker inspect -f '{{.State.Pid}}' "$cid")
  cg="/sys/fs/cgroup$(awk -F: '$1 == 0 {print $3}' "/proc/${pid}/cgroup" 2>/dev/null || true)"
  while [ "$(docker inspect -f '{{.State.Running}}' "$cid")" = true ]; do
    p=$(cat "${cg}/memory.peak" 2>/dev/null || echo 0)
    u=$(awk 'NR == FNR {cur = $1; next} $1 == "file" {f = $2} $1 == "shmem" {s = $2} END {print cur - f + s}' \
          "${cg}/memory.current" "${cg}/memory.stat" 2>/dev/null || echo 0)
    u=$((u + $(cat "${cg}/memory.swap.current" 2>/dev/null || echo 0)))
    [ "$((p / 1048576))" -le "$M_PEAK_MB" ] || M_PEAK_MB=$((p / 1048576))
    [ "$((u / 1048576))" -le "$M_USED_MB" ] || M_USED_MB=$((u / 1048576))
    sleep 1
  done
  M_RC=$(docker wait "$cid")
  M_SECONDS=$(($(date +%s) - t0))
  M_KILLED=$(docker inspect -f '{{.State.OOMKilled}}' "$cid")
  docker logs "$cid" 2>&1 | tail -n 15
  docker rm "$cid" >/dev/null
}

# row TASK LIMIT_MB REQUEST_GB RECORDED_MB EXPECT: one line of the table, and
# the checks on the run `measured` just made. EXPECT is pass, killed, or
# finish (exit 0 is enough: a rerun that only makes the next task's input).
row() {
  local verdict=ok need=$M_PEAK_MB killed=false
  [ "$M_USED_MB" -le "$need" ] || need=$M_USED_MB
  if [ "$M_RC" -eq 137 ] || [ "$M_KILLED" = true ]; then killed=true; fi
  if [ "$5" = killed ]; then
    $killed || verdict="FAIL: not killed for memory under ${2} MB (exit ${M_RC}, peak ${need} MB)"
  elif [ "$M_RC" -ne 0 ] || $killed; then
    verdict="FAIL: exit ${M_RC}, killed for memory: ${M_KILLED}"
  elif [ "$5" = finish ]; then
    :
  elif [ "$need" -eq 0 ]; then
    verdict="FAIL: no memory sample"
  elif [ $(($3 * 1024 * 100)) -lt $((need * HEADROOM_PCT)) ]; then
    verdict="FAIL: ${HEADROOM_PCT}% of ${need} MB is over the ${3} GB request"
  elif [ $((need * 100)) -gt $(($4 * 105)) ]; then
    verdict="FAIL: over the recorded ${4} MB by more than 5%; record the new figure in tests/test_prs_memory.sh and the docs"
  fi
  case "$verdict" in FAIL*) FAILS=$((FAILS + 1)) ;; esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "${MEASURE_SWAP_MB:-0}" "$M_RC" "$M_KILLED" "$M_PEAK_MB" "$M_USED_MB" "$M_SECONDS" "$verdict" | tee -a "${DIR}/table.tsv"
}

# limit_for PROCESS OVERRIDE_GB: the container limit in MB (the request, or the
# override, capped at CAP_MB) and the request in GB.
limit_for() {
  local gb=$2
  [ -n "$gb" ] || gb=$(request "$1" "${REPO}/conf/base.config") || { echo "ERROR: ${gb}" >&2; exit 1; }
  local mb=$((gb * 1024))
  [ "$mb" -le "$CAP_MB" ] || mb=$CAP_MB
  echo "$mb $gb"
}

echo "Host: $(nproc) CPUs, ${HOST_MB} MB; container limits are capped at ${CAP_MB} MB"
EXPECT=${EXPECT:-pass}
case "$EXPECT" in pass|killed) ;; *) echo "ERROR: EXPECT must be pass or killed, got '${EXPECT}'" >&2; exit 1 ;; esac
printf 'task\tlimit_mb\tswap_mb\texit\tkilled\tpeak_mb\tanon_plus_swap_mb\tseconds\tverdict\n' > "${DIR}/table.tsv"

# Control: the kill must be seen. 2 GB touched under a 1 GB limit.
MEASURE_SWAP_MB=0 measured 1024 "$PYTHON_IMAGE" python3 -c 'b = bytearray(2 << 30); print(len(b))'
# Docker's OOMKilled flag is not always set on cgroup v2 (one run had exit
# 137 and false), so exit 137 counts as a kill here and in row.
if [ "$M_RC" -eq 137 ]; then
  echo "Control: a 2 GB allocation under a 1 GB limit was killed (exit ${M_RC}, OOMKilled ${M_KILLED}) and reported."
else
  echo "ERROR: a 2 GB allocation under a 1 GB limit was not reported as killed (exit ${M_RC}, killed ${M_KILLED})" >&2
  exit 1
fi

# The scoring files, as setup.sh and step 25 fetch them.
PGS_BASE_URL=${PGS_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores}
N_SCORES=0
while IFS= read -r id; do
  url="${PGS_BASE_URL}/${id}/ScoringFiles/Harmonized/${id}_hmPOS_GRCh38.txt.gz"
  [ -f "${DIR}/scores/${id}.txt.gz" ] || fetch "$url" "${DIR}/scores/${id}.txt.gz" md5 "${url}.md5"
  N_SCORES=$((N_SCORES + 1))
done < <(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO}/assets/pgs_scores.tsv")
echo "${N_SCORES} scoring files, $(du -sm "${DIR}/scores" | cut -f1) MB"

# PRS_PREPARE, with the chrX rows (the larger of its two sets).
prepare() {
  rm -rf "${DIR}/prepare"
  mkdir -p "${DIR}/prepare"
  measured "$1" --network none -u "$(id -u):$(id -g)" \
    -v "${REPO}/bin:/pgp-bin:ro" -v "${REPO}/assets/pgs_scores.tsv:/pgs_scores.tsv:ro" \
    -v "${DIR}/scores:/scores:ro" -v "${DIR}/prepare:/out" "$PYTHON_IMAGE" \
    python3 /pgp-bin/collect_summary.py prs-format --scores /scores --labels /pgs_scores.tsv \
      --out /out/pgs --alleles /out/score_alleles.tsv --keep-x
}
read -r LIMIT GB < <(limit_for PRS_PREPARE "${PREPARE_LIMIT_GB:-}")
prepare "$LIMIT"
row PRS_PREPARE "$LIMIT" "$GB" "$RECORDED_PREPARE_MB" "$EXPECT"
if [ "$M_RC" -ne 0 ]; then
  echo "PRS_PREPARE failed under its limit; running it again with ${CAP_MB} MB for FORMAT_SCOREFILES' input."
  prepare "$CAP_MB"
  row "PRS_PREPARE (again, for the input)" "$CAP_MB" "$((CAP_MB / 1024))" "$CAP_MB" finish
  [ "$M_RC" -eq 0 ] || { echo "ERROR: PRS_PREPARE does not fit this machine" >&2; exit 1; }
fi

# FORMAT_SCOREFILES, the command of pgsc_calc's module at PGSC_CALC_VERSION.
read -r LIMIT GB < <(limit_for PRS "${FORMAT_LIMIT_GB:-}")
for threads in ${FORMAT_THREADS:-1}; do
  rm -rf "${DIR}/format"
  mkdir -p "${DIR}/format"
  measured "$LIMIT" --network none -u "$(id -u):$(id -g)" \
    -v "${DIR}/prepare/pgs:/pgs:ro" -v "${DIR}/format:/work" -w /work "$PGSC_UTILS_IMAGE" \
    /bin/bash -euo pipefail -c "mkdir formatted; pgscatalog-format -s /pgs/*.txt.gz -t GRCh38 --threads ${threads} -o formatted/ -l log_scorefiles.json -v"
  row "FORMAT_SCOREFILES, ${threads} thread(s)" "$LIMIT" "$GB" "$RECORDED_FORMAT_MB" "$EXPECT"
  if [ "$M_RC" -eq 0 ] && [ "$EXPECT" = pass ]; then
    n=$(find "${DIR}/format/formatted" -name 'normalised_*' | wc -l | tr -d ' ')
    [ "$n" -eq "$N_SCORES" ] || { echo "FAIL: ${n} normalised files, want ${N_SCORES}"; FAILS=$((FAILS + 1)); }
  fi
done

{
  echo "#### PRS memory on the default score list (${N_SCORES} scores), $(nproc) CPUs, ${HOST_MB} MB, expect ${EXPECT}"
  echo
  echo '```'
  cat "${DIR}/table.tsv"
  echo '```'
} | tee -a "$SUMMARY"
echo "${FAILS} failure(s)"
[ "$FAILS" -eq 0 ]
