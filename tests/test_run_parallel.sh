#!/usr/bin/env bash
# run_parallel (scripts/lib/common.sh) runs the scattered units of steps 03a,
# 03b and 29. This test runs it on units that exit 0, exit 1 and get killed,
# and checks what it returns, that every unit ran, that no more than MAX ran
# at once, and that it returned only after all had ended. The last check runs
# many short units that end together, the case where `wait -n` once returned
# 127 for a unit that had exited 0. It first proves it can fail: a copy of
# run_parallel that ignores each unit's status must be caught.
#
# RUN_PARALLEL_FROM=FILE takes the function from another file, for example an
# older common.sh.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FROM=${RUN_PARALLEL_FROM:-${REPO}/scripts/lib/common.sh}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
shopt -s nullglob

SRC=$(sed -n '/^run_parallel() {/,/^}/p' "$FROM")
[ -n "$SRC" ] || { echo "FAIL: no run_parallel() in ${FROM}"; exit 1; }

# unit SPEC: SPEC is NAME:ACTION. Marks itself live, records how many units are
# live, then exits 0 ("ok"), exits 1 ("fail") or kills itself with SIGKILL
# ("kill"); "fast" exits 0 at once. The live mark goes before the exit, so a
# count above MAX means run_parallel really ran more than MAX at once.
unit() {
  local name=${1%%:*} action=${1#*:} live
  : > "${D}/ran.${name}"
  [ "$action" = fast ] && exit 0
  : > "${D}/live.${name}"
  live=("$D"/live.*)
  echo "${#live[@]}" >> "${D}/counts"
  sleep 0.2
  rm -f "${D}/live.${name}"
  case "$action" in
    ok) exit 0 ;;
    fail) exit 1 ;;
    kill) kill -KILL "$BASHPID" ;;
  esac
}

# case_run LABEL WANT MAX SPEC...: run_parallel MAX unit SPEC..., then check it
# returned WANT (0 or 1), ran every unit, never more than MAX at once, and left
# none running. Prints one line; returns 1 on a mismatch.
case_run() {
  local label=$1 want=$2 max=$3 got=0 s name peak bad="" left
  shift 3
  D="${WORK}/${label}"
  rm -rf "$D"; mkdir -p "$D"; : > "${D}/counts"
  run_parallel "$max" unit "$@" 2>/dev/null || got=$?
  [ "$got" -eq "$want" ] || bad="returned ${got}, want ${want}"
  for s in "$@"; do
    name=${s%%:*}
    [ -e "${D}/ran.${name}" ] || bad="${bad:+${bad}; }unit ${name} never ran"
  done
  peak=$(sort -n "${D}/counts" | tail -1)
  [ "${peak:-0}" -le "$max" ] || bad="${bad:+${bad}; }${peak} units live at once, MAX ${max}"
  left=("$D"/live.*)
  if [ "${#left[@]}" -gt 0 ]; then bad="${bad:+${bad}; }returned with units still live"; fi
  if [ -n "$bad" ]; then echo "  not ok: ${label}: ${bad}"; return 1; fi
  echo "  ok:     ${label}: returned ${got}"
}

# check: every case against the run_parallel defined now; prints the number of
# failed cases and returns 1 when there is any.
check() {
  local f=0 i specs
  case_run all-exit-0 0 2 a:ok b:ok c:ok d:ok e:ok || f=$((f + 1))
  case_run one-exit-1-while-throttled 1 2 a:fail b:ok c:ok d:ok e:ok || f=$((f + 1))
  case_run one-exit-1-last 1 2 a:ok b:ok c:ok d:ok e:fail || f=$((f + 1))
  case_run one-killed 1 2 a:ok b:kill c:ok d:ok || f=$((f + 1))
  case_run exit-1-and-killed 1 3 a:ok b:fail c:kill d:ok e:ok f:ok || f=$((f + 1))
  case_run one-at-a-time 1 1 a:ok b:fail c:ok || f=$((f + 1))
  case_run more-slots-than-units 0 8 a:ok b:ok || f=$((f + 1))
  # 25 rounds of 24 units that exit at once, 4 at a time: all must count as 0.
  specs=()
  for i in $(seq 24); do specs+=("u${i}:fast"); done
  for i in $(seq 25); do
    case_run "fast-units-round-${i}" 0 4 "${specs[@]}" > "${WORK}/fast.out" || { cat "${WORK}/fast.out"; f=$((f + 1)); }
  done
  echo "  (25 rounds of 24 fast units checked)"
  echo "$f failed case(s)"
  [ "$f" -eq 0 ]
}

# Control: a run_parallel that drops each unit's status must be caught.
BROKEN=${SRC//|| rc=1/|| true}
if [ "$BROKEN" = "$SRC" ]; then
  BROKEN=${SRC//return \"\$rc\"/return 0}
fi
[ "$BROKEN" != "$SRC" ] || { echo "FAIL: could not plant the control (no '|| rc=1' or 'return \"\$rc\"' in run_parallel)"; exit 1; }
echo "Control: run_parallel with each unit's status ignored"
eval "$BROKEN"
if check; then
  echo "FAIL: a run_parallel that ignores failing units passed every case"
  exit 1
fi
echo "Control caught."

echo "run_parallel from ${FROM}"
eval "$SRC"
check
echo "PASS: run_parallel counts every unit's status"
