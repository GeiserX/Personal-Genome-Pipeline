#!/usr/bin/env bash
# check-exec-bits.sh [ROOT]: every bin/ helper the pipeline calls by name is
# committed with its executable bit, and starts with a #! line.
#
# Nextflow puts bin/ on the PATH of every task, so a module calls a helper by
# its bare name (`cdna_scan.py --out ...`). A helper committed as 100644 then
# fails only inside a real task ("Permission denied", exit 126): the stub runs
# and the unit tests, which call it through python3, never see it. That is how
# bin/cdna_scan.py reached main without the bit.
#
# A helper counts as called when its file name appears in main.nf,
# workflows/*.nf or modules/**/*.nf as a word that no "/" or word character
# comes before (`bin/foo.py` in a comment is a path, not a call). The mode is
# the one git records (git ls-files -s), which is what a checkout gets.
# ROOT is a git work tree, by default the checkout this script is in.
set -euo pipefail

ROOT=${1:-"$(cd "$(dirname "$0")/../.." && pwd)"}
cd "$ROOT"

mapfile -t nf < <(git ls-files -- main.nf 'workflows/*.nf' 'modules/*.nf' 'modules/**/*.nf' | LC_ALL=C sort -u)
[ "${#nf[@]}" -gt 0 ] || { echo "ERROR: no Nextflow file under ${ROOT}: nothing could call a helper" >&2; exit 1; }

bad=0 called=0
while read -r mode _ _ path; do
  name=${path#bin/}
  case "$name" in */*) continue ;; esac
  # The name as a regex: dots are literal.
  re="(^|[^/A-Za-z0-9_.-])${name//./\\.}([^A-Za-z0-9_.-]|\$)"
  hits=$(grep -nE -- "$re" "${nf[@]}" 2> /dev/null | head -n 1 || true)
  [ -n "$hits" ] || continue
  called=$((called + 1))
  if [ "$mode" != 100755 ]; then
    echo "FAIL ${path}: mode ${mode}, but ${hits%%:*} calls it by name (${hits#*:}); run: git update-index --chmod=+x ${path}"
    bad=$((bad + 1))
  elif [ "$(head -c 2 "$path")" != '#!' ]; then
    echo "FAIL ${path}: executable, but its first line is not a #! line, so the task cannot run it"
    bad=$((bad + 1))
  fi
done < <(git ls-files -s -- bin)

if [ "$called" -eq 0 ]; then
  echo "ERROR: no bin/ helper is called by name from ${#nf[@]} Nextflow files: the parser no longer sees the calls" >&2
  exit 1
fi
if [ "$bad" -gt 0 ]; then
  echo "${bad} of ${called} bin/ helpers the pipeline calls cannot be run by a task"
  exit 1
fi
echo "OK: the ${called} bin/ helpers the pipeline calls by name are executable"
