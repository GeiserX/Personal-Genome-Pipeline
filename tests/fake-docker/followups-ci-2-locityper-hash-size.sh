#!/usr/bin/env bash
# setup.sh --locityper-data sizes jellyfish's hash to the reference: the sum
# of the .fai lengths, at least 10 million entries. A reference over 2^31
# bases (GRCh38 is about 3.1G) keeps its full size: awk must not clamp it to
# a 32-bit integer. The run stops later, at the panel's index (the fake
# download does not match its sha256), after the count this case checks.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" FETCH_TRIES=1 FETCH_WAIT=0
seed_reference "$GENOME_DIR"
FAI="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
export LOCITYPER_LOCI="${CASE_WORK}/loci.bed"
printf 'chr1\t100\t200\tTEST\n' > "$LOCITYPER_LOCI"
use_output_hook

# size_for NAME FAI_ROWS WANT: run the install on a .fai of FAI_ROWS and check
# the --size jellyfish got.
size_for() {
  printf '%b' "$2" > "$FAI"
  : > "$FAKE_DOCKER_LOG"
  run_rc "$1" "${SCRIPTS}/setup.sh" --locityper-data "$GENOME_DIR"
  expect_rc "$1" 1
  output_has "$1" "Could not download the panel's index"
  docker_log_has "^run image=[^ ]*locityper.* jellyfish count .*--size $3 " \
    "$1: jellyfish did not get --size $3 ($(grep -o -- '--size [^ ]*' "$FAKE_DOCKER_LOG" || echo 'no jellyfish call'))"
}
size_for two-contigs 'chr1\t248956422\t6\t60\t61\nchr2\t242193529\t253105758\t60\t61\n' 491149951
size_for over-2g 'chr1\t2000000000\t6\t60\t61\nchr2\t1100000000\t2033333340\t60\t61\n' 3100000000
size_for floor 'chr1\t1000000\t6\t60\t61\n' 10000000
echo "PASS jellyfish's hash follows the reference length"
