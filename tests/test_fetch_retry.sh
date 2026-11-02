#!/usr/bin/env bash
# fetch (scripts/lib/common.sh) downloads a file again when a finished
# download fails its check, and gives up only after FETCH_TRIES tries. A fake
# _get serves each try from a script of answers: "good" (the right file),
# "bad" (a whole gzip file with the wrong content), "trunc" (the first bytes
# of the file, with no transfer error), "broken" (the first bytes, then a
# transfer error). The same cases run against a planted mutant that returns
# at the first failed check, the behaviour before this test: it must fail the
# "bad once, then good" case, or the test could not catch the old bug.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
export FETCH_WAIT=0

# The served file, a real gzip so the gzip -t check has something to pass.
printf 'line %s\n' $(seq 1 2000) | gzip -c > "${WORK}/good.gz"
printf 'not the file\n' | gzip -c > "${WORK}/bad.gz"
GOOD_MD5=$(bash -c '. "$1" && _digest md5 "$2"' _ "${REPO}/scripts/lib/common.sh" "${WORK}/good.gz")

# A mutant copy of the library: a failed check returns at once, no new try.
mkdir -p "${WORK}/mutant/scripts/lib"
cp "${REPO}/versions.env" "${WORK}/mutant/"
# shellcheck disable=SC2016  # the $ is sed's, not the shell's
sed 's/^\(      rm -f "\$part"\)$/\1; return 1/' "${REPO}/scripts/lib/common.sh" > "${WORK}/mutant/scripts/lib/common.sh"
if cmp -s "${REPO}/scripts/lib/common.sh" "${WORK}/mutant/scripts/lib/common.sh"; then
  echo "FAIL: the mutant is the same as common.sh; the sed no longer matches fetch"
  exit 1
fi

# run_case LIB NAME ANSWERS KIND: fetch with the fake _get answering ANSWERS
# (space-separated, one per try) and check KIND (md5 or gz). Prints
# "rc=<fetch exit> calls=<tries used> dest=<good|bad|none> part=<yes|no>".
run_case() {
  local lib=$1 case_id=$2 answers=$3 kind=$4
  (
    # shellcheck source=../scripts/lib/common.sh
    . "$lib"
    n=0
    # shellcheck disable=SC2329,SC2317  # called by fetch
    _get() {
      local a
      n=$((n + 1))
      a=$(echo "$answers" | awk -v i="$n" '{print $i}')
      echo "$n" > "${WORK}/${case_id}.calls"
      case "$a" in
        good) cat "${WORK}/good.gz" > "$2" ;;
        bad) cat "${WORK}/bad.gz" > "$2" ;;
        trunc) head -c 200 "${WORK}/good.gz" > "$2" ;;
        broken) head -c 100 "${WORK}/good.gz" > "$2"; return 18 ;;
        *) return 22 ;;
      esac
    }
    dest="${WORK}/${case_id}/file.gz"
    rc=0
    if [ "$kind" = md5 ]; then
      fetch https://example.org/file.gz "$dest" md5 "$GOOD_MD5" 2>"${WORK}/${case_id}.err" || rc=$?
    else
      fetch https://example.org/file.gz "$dest" 2>"${WORK}/${case_id}.err" || rc=$?
    fi
    d=none
    if [ -e "$dest" ]; then
      if cmp -s "$dest" "${WORK}/good.gz"; then d=good; else d=bad; fi
    fi
    p=no; [ -e "${dest}.part" ] && p=yes
    echo "rc=${rc} calls=$(cat "${WORK}/${case_id}.calls" 2>/dev/null || echo 0) dest=${d} part=${p}"
  )
}

FAILS=0
expect() {
  local label=$1 want=$2 got=$3
  if [ "$got" = "$want" ]; then
    echo "ok:   ${label}: ${got}"
  else
    echo "FAIL: ${label}: got '${got}', expected '${want}'"; FAILS=$((FAILS + 1))
  fi
}

LIB="${REPO}/scripts/lib/common.sh"
export FETCH_TRIES=3
expect "md5 bad once, then good" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" md5-once 'bad good' md5)"
expect "md5 bad every time" "rc=1 calls=3 dest=none part=no" "$(run_case "$LIB" md5-always 'bad bad bad' md5)"
expect "gzip bad once, then good" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" gz-once 'trunc good' gz)"
expect "gzip bad every time" "rc=1 calls=3 dest=none part=no" "$(run_case "$LIB" gz-always 'trunc trunc trunc' gz)"
expect "broken transfer, then good" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" broken-once 'broken good' md5)"
expect "broken every time keeps the .part" "rc=1 calls=3 dest=none part=yes" "$(run_case "$LIB" broken-always 'broken broken broken' md5)"
expect "one try, bad: no second download" "rc=1 calls=1 dest=none part=no" "$(FETCH_TRIES=1 run_case "$LIB" one-try 'bad good' md5)"

if grep -q 'ERROR: md5 checksum of file\.gz is [0-9a-f]*, expected '"${GOOD_MD5}" "${WORK}/md5-always.err"; then
  echo "ok:   the last failed check names both checksums"
else
  echo "FAIL: no checksum error after the last try:"; cat "${WORK}/md5-always.err"; FAILS=$((FAILS + 1))
fi
if [ "$(grep -c 'failed its check' "${WORK}/md5-always.err")" -eq 3 ]; then
  echo "ok:   each failed check is reported"
else
  echo "FAIL: expected 3 'failed its check' lines:"; cat "${WORK}/md5-always.err"; FAILS=$((FAILS + 1))
fi

# Control: the mutant must fail the retry case.
got=$(run_case "${WORK}/mutant/scripts/lib/common.sh" mutant 'bad good' md5)
if [ "$got" = "rc=1 calls=1 dest=none part=no" ]; then
  echo "ok:   control: the no-retry mutant fails 'bad once, then good' (${got})"
else
  echo "FAIL: control: the no-retry mutant gave '${got}'; the test cannot tell it from the fix"; FAILS=$((FAILS + 1))
fi

[ "$FAILS" -eq 0 ]
