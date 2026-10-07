#!/usr/bin/env bash
# Step 26 without the ancestry panel prints one line, exits 0, starts no
# container and downloads nothing. With a panel but no site list beside it,
# step 25 (which step 26 runs) stops before it starts pgsc_calc, naming
# setup.sh --ancestry-panel.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

G="${CASE_WORK}/genome"
mkdir -p "$G"
export GENOME_DIR="$G"
seed_sample "$G" sample1

# --- no panel ------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 no-panel "${SCRIPTS}/26-ancestry.sh" sample1
[ "$(grep -c . "${CASE_WORK}/no-panel.out")" -eq 1 ] || fail "step 26 without a panel printed more than one line: $(cat "${CASE_WORK}/no-panel.out")"
output_has no-panel '^Step 26 skipped: no ancestry reference panel at .*/reference/pgsc_calc/'"${PGSC_PANEL}"'\.tar\.zst; install it with scripts/setup\.sh --ancestry-panel '
if grep -qE '^(curl|wget|run) ' "$FAKE_DOCKER_LOG"; then
  fail "step 26 without a panel downloaded or started something: $(grep -E '^(curl|wget|run) ' "$FAKE_DOCKER_LOG" | head -3)"
fi
[ ! -e "${G}/sample1/ancestry" ] || fail "step 26 without a panel made ${G}/sample1/ancestry"

# ANCESTRY_PANEL=none skips too, even with a panel installed.
mkdir -p "${G}/reference/pgsc_calc"
printf 'placeholder\n' > "${G}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst"
run_expect 0 panel-none env ANCESTRY_PANEL=none "${SCRIPTS}/26-ancestry.sh" sample1
output_has panel-none '^Step 26 skipped: '

# --- a panel without its site list ----------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 1 no-sites "${SCRIPTS}/26-ancestry.sh" sample1
output_has no-sites 'has no site list beside it'
output_has no-sites 'setup\.sh --ancestry-panel'
if grep -q '^nextflow ' "$FAKE_DOCKER_LOG"; then
  fail "step 25 started pgsc_calc although the panel has no site list"
fi
