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
