#!/usr/bin/env bash
# finish in tests/e2e/lib.sh: an e2e case that ran no check fails, unless it
# set NO_CHECKS_REASON, which finish prints and adds to the notes. A case
# with a passing check passes, one with a failed check fails. Planted cases,
# no docker, no network.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
failed=0
mkdir -p "${T}/work" "${T}/genome" "${T}/fixture" "${T}/cases"
# expect DESC EXIT BODY: a case of BODY between lib.sh and finish exits EXIT.
expect() {
  local f="${T}/cases/planted-$RANDOM.sh" rc=0
  printf '. "%s/tests/e2e/lib.sh"\n%s\nfinish\n' "$REPO" "$3" > "$f"
  : > "${T}/notes.md"
  REPO=$REPO GENOME_DIR="${T}/genome" SAMPLE=HG002 FIXTURE_DIR="${T}/fixture" E2E_WORK="${T}/work" \
    E2E_NOTES="${T}/notes.md" bash "$f" > "${T}/out" 2>&1 || rc=$?
  if [ "$rc" = "$2" ]; then echo "PASS: $1 (exit ${rc}: $(tail -n 1 "${T}/out"))"
  else echo "FAIL: $1 (exit ${rc}, want $2): $(cat "${T}/out")"; failed=1; fi
}
expect "a case that ran no check fails" 1 ':'
if grep -q "no check ran, so the case proves nothing" "${T}/out"; then echo "PASS: and says it ran no check"
else echo "FAIL: no message for a case with no check"; failed=1; fi
expect "a case with no check and a reason passes" 0 'NO_CHECKS_REASON="nothing to measure on this run"'
if grep -q "no checks run: nothing to measure on this run" "${T}/notes.md"; then echo "PASS: and its reason is in the notes"
else echo "FAIL: the reason is not in the notes"; failed=1; fi
expect "a case with one passing check passes" 0 'check_eq "one" 1 1'
expect "a case with a failed check fails" 1 'check_eq "one" 1 2'
expect "a failed check fails the case even with a reason set" 1 'NO_CHECKS_REASON=x; check_eq "one" 1 2'
exit "$failed"
