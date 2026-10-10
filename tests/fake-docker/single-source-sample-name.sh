#!/usr/bin/env bash
# Every script that takes a sample name refuses one that is not a plain name,
# with exit 2, before it starts a container or creates a directory. The name
# goes into container paths and `bash -c` bodies, so '../x' or 'x;id' must
# never get that far. Nor may '-x': CPSR, bin/pgx_parse.py and
# bin/collect_summary.py take the name as an argument value, and argparse reads
# '-x' as an option. A plain name with '.', '_' and '-' is accepted, and so is
# a two-character one (CPSR's 3-character minimum is handled at step 17).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
mkdir -p "$GENOME_DIR"
export PLATFORM=ont   # the long-read scripts ask for it before anything else

checked=0
for s in "${SCRIPTS}"/*.sh; do
  name=$(basename "$s")
  [ "$name" != setup.sh ] || continue   # takes a directory, not a sample
  for bad in '../x' '..' '.' 'a b' 'x;id' 'x$(id)' '-x' ''; do
    [ -n "$bad" ] || [ "$name" != validate-setup.sh ] || continue   # no sample is valid there
    rc=0
    "$s" "$bad" male extra > "${CASE_WORK}/name.out" 2>&1 || rc=$?
    if [ -z "$bad" ]; then
      # An empty name is refused by the usage check or by validate_sample.
      [ "$rc" -ne 0 ] || fail "${name} accepted an empty sample name"
      continue
    fi
    if [ "$rc" -ne 2 ] || ! grep -q 'invalid sample name' "${CASE_WORK}/name.out"; then
      cat "${CASE_WORK}/name.out"
      fail "${name} '${bad}' exited ${rc} without 'invalid sample name', expected exit 2"
    fi
  done
  checked=$((checked + 1))
done
echo "checked ${checked} scripts"
[ "$checked" -ge 50 ] || fail "only ${checked} scripts were checked; the loop lost most of scripts/*.sh"

if awk '/^(run|pull) / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "a script started a container although the sample name was invalid"
fi
[ -z "$(ls -A "$GENOME_DIR")" ] || fail "a script created $(ls -A "$GENOME_DIR") in GENOME_DIR for an invalid sample name"

# Control: a valid name is not refused (step 11 then stops on the missing VCF).
run_rc good-name "${SCRIPTS}/11-roh-analysis.sh" 'sample.1_A-b'
output_lacks good-name 'invalid sample name'
[ "$RC" -ne 2 ] || fail "a valid sample name was refused"
run_rc short-name "${SCRIPTS}/11-roh-analysis.sh" S1
output_lacks short-name 'invalid sample name'
[ "$RC" -ne 2 ] || fail "a two-character sample name was refused"
