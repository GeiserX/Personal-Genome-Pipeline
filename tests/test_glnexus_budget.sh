#!/usr/bin/env bash
# GLnexus's --mem-gbytes sizes only its RocksDB (block cache and bulk-load
# write buffers); the alleles it discovers and the unifier's sites are held
# outside it, so a budget near the limit leaves no room for them. Both entry
# points must give GLnexus at most half its container's memory, by the same
# rule:
#   scripts/38-family.sh           --mem-gbytes "$((<expr of GLNEXUS_MEM>))"
#   modules/local/glnexus/main.nf  def mem = <expr of task.memory.toGiga()>
# This test reads each expression from the source, evaluates it for container
# limits of 1 to 64 GB, and checks both. It first proves it can fail: the
# limit-minus-1 rule, and a module under half but out of step with the
# script, must both be reported.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# bash_expr SCRIPT: the arithmetic passed to --mem-gbytes, without $(( ))
bash_expr() {
  # shellcheck disable=SC2016 # the $(( )) is matched, not expanded
  grep -o -- '--mem-gbytes "\$((.*))"' "$1" | sed -e 's/^--mem-gbytes "\$((//' -e 's/))"$//' || true
}
# nf_expr MODULE: the right-hand side of "def mem = ...", as Python over L
nf_expr() {
  grep -E '^[[:space:]]*def mem = ' "$1" | sed -e 's/^[[:space:]]*def mem = //' \
    -e 's/Math\.max(/max(/g' -e 's/(task\.memory\.toGiga() as int)/L/g' \
    -e 's/\.intdiv(\([0-9]*\))/ \/\/ \1/g' || true
}

# check SCRIPT MODULE: 0 when both budgets are at most half the limit (at
# least 1) and equal for every limit from 1 to 64 GB
check() {
  local be ne l b n bad=0
  be=$(bash_expr "$1"); ne=$(nf_expr "$2")
  [ -n "$be" ] || { echo "no --mem-gbytes \"\$((...))\" in $1"; return 1; }
  [ -n "$ne" ] || { echo "no 'def mem = ...' in $2"; return 1; }
  # Only max, L, integers and arithmetic may be left of the Groovy expression.
  if [ -n "$(printf '%s' "$ne" | sed 's/max//g' | tr -cd 'A-Za-z_' | tr -d 'L')" ]; then
    echo "cannot read the Groovy budget '$ne' in $2"; return 1
  fi
  for l in $(seq 1 64); do
    b=$(GLNEXUS_MEM=$l bash -c "echo \$(( $be ))")
    n=$(python3 -c "L=$l; print($ne)")
    if [ "$b" -lt 1 ] || [ "$n" -lt 1 ]; then
      echo "limit ${l} GB: budget below 1 GB (script ${b}, module ${n})"; bad=1
    elif [ "$l" -ge 2 ] && { [ $((b * 2)) -gt "$l" ] || [ $((n * 2)) -gt "$l" ]; }; then
      echo "limit ${l} GB: budget over half the limit (script ${b}, module ${n})"; bad=1
    elif [ "$b" -ne "$n" ]; then
      echo "limit ${l} GB: the two entry points differ (script ${b}, module ${n})"; bad=1
    fi
  done
  [ "$bad" -eq 0 ] || return 1
  echo "budget at 8 / 16 / 32 GB: $(for l in 8 16 32; do GLNEXUS_MEM=$l bash -c "echo -n \$(( $be ))"; echo -n ' '; done)GB in both"
}

FAILS=0
# Control: the rule this test replaced (limit minus 1) must be reported.
# shellcheck disable=SC2016 # the $(( )) is written into the copy, not expanded
sed 's/--mem-gbytes "\$((.*))"/--mem-gbytes "$((GLNEXUS_MEM > 1 ? GLNEXUS_MEM - 1 : 1))"/' \
  "${REPO}/scripts/38-family.sh" > "${WORK}/old.sh"
sed 's/^\([[:space:]]*def mem = \).*/\1Math.max(1, (task.memory.toGiga() as int) - 1)/' \
  "${REPO}/modules/local/glnexus/main.nf" > "${WORK}/old.nf"
if out=$(check "${WORK}/old.sh" "${WORK}/old.nf"); then
  echo "FAIL: the limit-minus-1 budget was not reported"; FAILS=$((FAILS + 1))
else
  echo "ok:   limit-minus-1 budget reported: $(echo "$out" | head -1)"
fi
# Control: a module under half but out of step with the script must be
# reported by the equality check, not the half rule.
sed 's/^\([[:space:]]*def mem = \).*/\1Math.max(1, (task.memory.toGiga() as int).intdiv(3))/' \
  "${REPO}/modules/local/glnexus/main.nf" > "${WORK}/third.nf"
if out=$(check "${REPO}/scripts/38-family.sh" "${WORK}/third.nf"); then
  echo "FAIL: a module out of step with the script was not reported"; FAILS=$((FAILS + 1))
elif ! printf '%s\n' "$out" | grep -q 'differ'; then
  echo "FAIL: a module out of step was not reported as differing: $(echo "$out" | head -1)"; FAILS=$((FAILS + 1))
else
  echo "ok:   module out of step reported: $(echo "$out" | head -1)"
fi

if out=$(check "${REPO}/scripts/38-family.sh" "${REPO}/modules/local/glnexus/main.nf"); then
  echo "ok:   ${out}"
else
  echo "FAIL:"; echo "$out" | head -5; FAILS=$((FAILS + 1))
fi
[ "$FAILS" -eq 0 ]
