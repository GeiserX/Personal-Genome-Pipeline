#!/usr/bin/env bash
# The exit codes that mean "killed for memory or time" are retried by two
# errorStrategy closures: conf/base.config's process-wide one and
# nextflow.config's withName block for the report-only tools. Nextflow has no
# shared constant across the two files, so this test is the single check that
# both lists are the same. It first proves it can fail: a copy with one list
# changed must be reported as a mismatch, and a copy with no list at all must
# be reported as missing.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# lists FILE: every "task.exitStatus in [..]" list in FILE, one per line, spaces removed
lists() { grep -o 'task\.exitStatus in \[[^]]*\]' "$1" | sed -e 's/.*\[//' -e 's/\]//' -e 's/ //g' || true; }

# compare BASE NEXTFLOW: 0 when each file has at least one list and all lists are equal
compare() {
  local a b
  a=$(lists "$1"); b=$(lists "$2")
  [ -n "$a" ] || { echo "no 'task.exitStatus in [...]' list in $1"; return 1; }
  [ -n "$b" ] || { echo "no 'task.exitStatus in [...]' list in $2"; return 1; }
  if [ "$(printf '%s\n%s\n' "$a" "$b" | sort -u | wc -l)" -ne 1 ]; then
    echo "retry exit codes differ:"; echo "  $1: $(echo "$a" | paste -sd' ' -)"; echo "  $2: $(echo "$b" | paste -sd' ' -)"
    return 1
  fi
  echo "both files retry exit codes [$a]"
}

FAILS=0
# Controls: a planted mismatch and a planted missing list must both fail.
sed 's/task\.exitStatus in \[[^]]*\]/task.exitStatus in [1]/' "${REPO}/nextflow.config" > "${WORK}/mismatch.config"
if out=$(compare "${REPO}/conf/base.config" "${WORK}/mismatch.config"); then
  echo "FAIL: a planted mismatch was not reported"; FAILS=$((FAILS + 1))
else
  echo "ok:   planted mismatch reported: ${out//$'\n'/ }"
fi
sed 's/task\.exitStatus in/task.exitStatus inside/' "${REPO}/conf/base.config" > "${WORK}/none.config"
if out=$(compare "${WORK}/none.config" "${REPO}/nextflow.config"); then
  echo "FAIL: a file with no list was not reported"; FAILS=$((FAILS + 1))
else
  echo "ok:   planted missing list reported: ${out}"
fi

if out=$(compare "${REPO}/conf/base.config" "${REPO}/nextflow.config"); then
  echo "ok:   ${out}"
else
  echo "FAIL: ${out//$'\n'/ }"; FAILS=$((FAILS + 1))
fi
[ "$FAILS" -eq 0 ]
