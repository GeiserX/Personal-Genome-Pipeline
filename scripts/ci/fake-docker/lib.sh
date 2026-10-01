#!/usr/bin/env bash
# lib.sh: helpers for the case files in tests/fake-docker/*.sh.
#
# scripts/ci/run-fake-docker-suite.sh runs each case in a clean environment
# (env -i) with:
#   REPO_ROOT        the repository checkout
#   CASE_WORK        an empty temp directory owned by this case
#   FAKE_DOCKER_LOG  the log the fake docker, wget and curl append to
#   PATH             scripts/ci/fake-docker first, then the system paths
# A case passes when it exits 0. The runner also fails it when the docker log
# shows an empty image name, a missing image or an unknown docker option.
set -euo pipefail

: "${REPO_ROOT:?run cases through scripts/ci/run-fake-docker-suite.sh}"
: "${CASE_WORK:?run cases through scripts/ci/run-fake-docker-suite.sh}"
: "${FAKE_DOCKER_LOG:?run cases through scripts/ci/run-fake-docker-suite.sh}"

# shellcheck disable=SC2034  # read by the case files
SCRIPTS="${REPO_ROOT}/scripts"

fail() {
  echo "ASSERT FAIL: $*" >&2
  exit 1
}

# seed_reference GENOME_DIR: placeholder GRCh38 FASTA and .fai. The FASTA is a
# sparse file with an apparent size of 3.2 GB (no disk used) because
# validate-setup.sh rejects a reference smaller than 3 GB.
seed_reference() {
  mkdir -p "$1/reference"
  truncate -s 3200000000 "$1/reference/Homo_sapiens_assembly38.fasta"
  printf 'chr1\t248956422\t6\t60\t61\n' > "$1/reference/Homo_sapiens_assembly38.fasta.fai"
}

# seed_clinvar GENOME_DIR: the raw ClinVar download and the two files setup.sh
# derives from it, each with its index.
seed_clinvar() {
  mkdir -p "$1/clinvar"
  local f
  for f in clinvar clinvar_chr clinvar_pathogenic_chr; do
    printf 'placeholder\n' > "$1/clinvar/${f}.vcf.gz"
    printf 'placeholder\n' > "$1/clinvar/${f}.vcf.gz.tbi"
  done
}

# seed_sample GENOME_DIR SAMPLE: placeholder sorted BAM and VCF with indexes,
# so run-all.sh skips alignment and variant calling.
seed_sample() {
  local g=$1 s=$2
  mkdir -p "${g}/${s}/aligned" "${g}/${s}/vcf"
  printf 'placeholder\n' > "${g}/${s}/aligned/${s}_sorted.bam"
  printf 'placeholder\n' > "${g}/${s}/aligned/${s}_sorted.bam.bai"
  printf 'placeholder\n' > "${g}/${s}/vcf/${s}.vcf.gz"
  printf 'placeholder\n' > "${g}/${s}/vcf/${s}.vcf.gz.tbi"
}

# run_expect CODE NAME COMMAND [ARGS...]: run COMMAND, keep its output in
# ${CASE_WORK}/NAME.out, print it, and fail unless it exited with CODE.
run_expect() {
  local want=$1 name=$2 rc=0
  shift 2
  echo "--- ${name}: $* (expect exit ${want})"
  "$@" > "${CASE_WORK}/${name}.out" 2>&1 || rc=$?
  cat "${CASE_WORK}/${name}.out"
  echo "--- ${name}: exit ${rc}"
  [ "$rc" -eq "$want" ] || fail "${name} exited ${rc}, expected ${want}"
}

# run_rc NAME COMMAND [ARGS...]: like run_expect, but only records the exit
# code in RC, so a case can check what the command did before its exit code.
# expect_rc NAME CODE: fail unless the last run_rc exited with CODE.
RC=0
run_rc() {
  local name=$1
  shift
  RC=0
  echo "--- ${name}: $*"
  "$@" > "${CASE_WORK}/${name}.out" 2>&1 || RC=$?
  cat "${CASE_WORK}/${name}.out"
  echo "--- ${name}: exit ${RC}"
}
expect_rc() {
  [ "$RC" -eq "$2" ] || fail "$1 exited ${RC}, expected $2"
}

# docker_log_has AWK_REGEX MESSAGE: fail with MESSAGE unless a line of the
# docker/wget/curl log matches. One awk process, so no pipe can lose a match.
docker_log_has() {
  awk -v re="$1" '$0 ~ re { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG" || fail "$2"
}

# hide_commands NAME...: point PATH at one directory of links to every command
# on the current PATH except NAME..., so `command -v NAME` fails.
hide_commands() {
  local bin="${CASE_WORK}/path-without" d f n h skip
  mkdir -p "$bin"
  local -a dirs
  IFS=: read -r -a dirs <<<"$PATH"
  for d in "${dirs[@]}"; do
    for f in "$d"/*; do
      n=${f##*/}
      if [ -d "$f" ] || [ ! -x "$f" ] || [ -e "${bin}/${n}" ] || [ -L "${bin}/${n}" ]; then continue; fi
      skip=false
      for h in "$@"; do [ "$n" != "$h" ] || skip=true; done
      $skip || ln -s "$f" "${bin}/${n}"
    done
  done
  export PATH="$bin"
  hash -r
  for h in "$@"; do
    if command -v "$h" >/dev/null 2>&1; then fail "hide_commands: ${h} is still on PATH"; fi
  done
}

# use_output_hook: make every `docker run` create the files its command
# writes under /genome (after -o, --output or >), plus the .tbi of every
# `index` target, in GENOME_DIR. For scripts that check a step's output.
use_output_hook() {
  cat > "${CASE_WORK}/hook-outputs" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
shift   # the image
args="$*"
mk() {
  local h="${GENOME_DIR:?}${1#/genome}"
  mkdir -p "$(dirname "$h")"
  case "$h" in
    *.gz) printf '##fileformat=VCFv4.2\n' | gzip -c > "$h" ;;
    *) : > "$h" ;;
  esac
}
path='/genome/[^[:space:];&|"'\'']+'
while read -r p; do
  [ -n "$p" ] && mk "$p"
done < <(grep -oE -- "(-o|--output|>)[[:space:]]*${path}" <<<"$args" | grep -oE -- "${path}" || true)
while read -r p; do
  [ -n "$p" ] && mk "${p}.tbi"
done < <(grep -oE -- "index[[:space:]]+(-t[[:space:]]+)?${path}" <<<"$args" | grep -oE -- "${path}" || true)
exit 0
HOOK
  chmod +x "${CASE_WORK}/hook-outputs"
  export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-outputs"
}

# output_has NAME REGEX / output_lacks NAME REGEX: check a run_expect output.
output_has() {
  grep -qE -- "$2" "${CASE_WORK}/$1.out" || fail "$1 output has no line matching: $2"
}
output_lacks() {
  local hit
  if hit=$(grep -nE -- "$2" "${CASE_WORK}/$1.out"); then
    fail "$1 output matches '$2': ${hit}"
  fi
}
