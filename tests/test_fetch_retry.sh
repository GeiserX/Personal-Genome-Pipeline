#!/usr/bin/env bash
# fetch (scripts/lib/common.sh) downloads a file again when a finished
# download fails its check, and gives up only after FETCH_TRIES tries. A fake
# _get serves each try from a script of answers: "good" (the right file),
# "bad" (a whole gzip file with the wrong content), "trunc" (the first bytes
# of the file, with no transfer error), "broken" (the first bytes, then a
# transfer error), "fresh" (the right file, only if no .part is there) and
# "resume" (the right file, only if the 100-byte .part of "broken" is there).
# Two planted mutants show the test can fail: one returns at the first failed
# check, the behaviour before this test, and one drops the .part of a broken
# transfer instead of resuming it.
#
# Stalls: _get, which fetch calls for each try, gives up on a transfer that
# stops sending data, so fetch can try again and resume the .part. A local
# HTTP server sends the start of a file and then nothing, holding the
# connection open. With FETCH_STALL=2, fetch with curl and fetch with wget
# (curl hidden from PATH) each fail within seconds, keep the .part, and try
# FETCH_TRIES times; the same server's whole file downloads under the same
# bounds (at most 20 s). Two planted mutants, curl and wget without the
# bounds, still hang at 25 s, when this test stops them: the behaviour before
# the bounds.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
SRV=""
trap '[ -z "$SRV" ] || kill "$SRV" 2>/dev/null; rm -rf "$WORK"' EXIT
export FETCH_WAIT=0

# The served file, a real gzip so the gzip -t check has something to pass.
printf 'line %s\n' $(seq 1 2000) | gzip -c > "${WORK}/good.gz"
printf 'not the file\n' | gzip -c > "${WORK}/bad.gz"
GOOD_MD5=$(bash -c '. "$1" && _digest md5 "$2"' _ "${REPO}/scripts/lib/common.sh" "${WORK}/good.gz")

# mutant NAME SED: a copy of the library with SED applied to it, at
# ${WORK}/NAME/scripts/lib/common.sh; exits when SED changes nothing.
mutant() {
  mkdir -p "${WORK}/$1/scripts/lib"
  cp "${REPO}/versions.env" "${WORK}/$1/"
  sed "$2" "${REPO}/scripts/lib/common.sh" > "${WORK}/$1/scripts/lib/common.sh"
  if cmp -s "${REPO}/scripts/lib/common.sh" "${WORK}/$1/scripts/lib/common.sh"; then
    echo "FAIL: mutant $1 is the same as common.sh; its sed no longer matches fetch"
    exit 1
  fi
}
# shellcheck disable=SC2016  # the $ is sed's, not the shell's
mutant no-retry 's/^\(      rm -f "\$part"\)$/\1; return 1/'
# shellcheck disable=SC2016
mutant no-resume 's/^\(      bad=""\)$/\1; rm -f "$part"/'

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
        fresh) [ ! -e "$2" ] || return 22; cat "${WORK}/good.gz" > "$2" ;;
        resume) [ "$(wc -c < "$2" 2>/dev/null)" -eq 100 ] || return 22; cat "${WORK}/good.gz" > "$2" ;;
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
expect "md5 bad once, then a fresh download" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" md5-once 'bad fresh' md5)"
expect "md5 bad every time" "rc=1 calls=3 dest=none part=no" "$(run_case "$LIB" md5-always 'bad bad bad' md5)"
expect "gzip bad once, then a fresh download" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" gz-once 'trunc fresh' gz)"
expect "gzip bad every time" "rc=1 calls=3 dest=none part=no" "$(run_case "$LIB" gz-always 'trunc trunc trunc' gz)"
expect "broken transfer, then resumed" "rc=0 calls=2 dest=good part=no" "$(run_case "$LIB" broken-once 'broken resume' md5)"
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

# Controls: each mutant must fail the case it breaks.
got=$(run_case "${WORK}/no-retry/scripts/lib/common.sh" no-retry 'bad fresh' md5)
if [ "$got" = "rc=1 calls=1 dest=none part=no" ]; then
  echo "ok:   control: the no-retry mutant fails 'bad once, then a fresh download' (${got})"
else
  echo "FAIL: control: the no-retry mutant gave '${got}'; the test cannot tell it from the fix"; FAILS=$((FAILS + 1))
fi
got=$(run_case "${WORK}/no-resume/scripts/lib/common.sh" no-resume 'broken resume' md5)
if [ "$got" != "rc=0 calls=2 dest=good part=no" ]; then
  echo "ok:   control: the no-resume mutant fails 'broken transfer, then resumed' (${got})"
else
  echo "FAIL: control: the no-resume mutant passed; the test cannot see a dropped .part"; FAILS=$((FAILS + 1))
fi

# --- Stalls -------------------------------------------------------------------
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }
command -v curl >/dev/null || { echo "FAIL: curl is needed to run this test"; exit 1; }
command -v wget >/dev/null || { echo "FAIL: wget is needed to run this test"; exit 1; }

# The server: /stall.bin sends a 100000-byte Content-Length, 1000 bytes, then
# holds the connection for 60 s; /whole.bin sends the whole file.
cat > "${WORK}/server.py" <<'PY'
import http.server, sys, time
class H(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-Length", "100000")
        self.end_headers()
        if self.path == "/stall.bin":
            self.wfile.write(b"x" * 1000); self.wfile.flush(); time.sleep(60)
        else:
            self.wfile.write(b"x" * 100000)
s = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
open(sys.argv[1], "w").write(str(s.server_address[1]))
s.serve_forever()
PY
python3 "${WORK}/server.py" "${WORK}/port" & SRV=$!
for _ in $(seq 1 50); do [ -s "${WORK}/port" ] && break; sleep 0.1; done
[ -s "${WORK}/port" ] || { echo "FAIL: the test server did not start"; exit 1; }
URL="http://127.0.0.1:$(cat "${WORK}/port")"

# A PATH without curl, so _get takes its wget branch.
mkdir -p "${WORK}/nocurl"
for c in bash wget awk sed mv rm mkdir cat sleep basename dirname gzip md5sum sha256sum wc head tail tr cut date env \
         grep sort stat readlink uname id timeout nproc getconf printf; do
  p=$(command -v "$c") && ln -sf "$p" "${WORK}/nocurl/$c"
done

# stall_mutant NAME SED: a copy of the library with SED applied, at
# ${WORK}/NAME/scripts/lib/common.sh; stops when SED changes nothing.
stall_mutant() {
  mkdir -p "${WORK}/$1/scripts/lib"
  cp "${REPO}/versions.env" "${WORK}/$1/"
  sed "$2" "${REPO}/scripts/lib/common.sh" > "${WORK}/$1/scripts/lib/common.sh"
  if cmp -s "${REPO}/scripts/lib/common.sh" "${WORK}/$1/scripts/lib/common.sh"; then
    echo "FAIL: mutant $1 is the same as common.sh; its sed no longer matches _get"
    exit 1
  fi
}
# shellcheck disable=SC2016  # the $ and braces are sed's, not the shell's
stall_mutant curl-unbounded 's/ "\${cbound\[@\]}"//g'
# shellcheck disable=SC2016
stall_mutant wget-unbounded 's/ "\${wbound\[@\]}"//g'

# fetch_in LIB PATH NAME FILE [LIMIT]: fetch URL/FILE with LIB, PATH and
# FETCH_STALL=2, two tries, under a LIMIT-second timeout (default 30).
# Prints "rc=<exit> secs=<seconds> part=<bytes>".
fetch_in() {
  local start rc=0 dest="${WORK}/$3/$4"
  start=$(date +%s)
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  PATH="$2" FETCH_STALL=2 FETCH_TRIES=2 FETCH_WAIT=0 timeout -k 2 "${5:-30}" bash -c '. "$1" && fetch "$2" "$3"' \
    _ "$1" "${URL}/$4" "$dest" > "${WORK}/$3.log" 2>&1 || rc=$?
  echo "rc=${rc} secs=$(( $(date +%s) - start )) part=$( { wc -c < "${dest}.part" || echo 0; } 2>/dev/null | tr -d ' ')"
}
# expect_re DESC GOT WANT [LOG]: GOT matches the extended regex WANT.
expect_re() {
  if [[ "$2" =~ $3 ]]; then pass "$1 ($2)"; else fail "$1: got '$2', want /$3/"; sed 's/^/  /' "${WORK}/${4:-none}.log" 2>/dev/null | tail -n 5; fi
}

expect_re "curl: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "$PATH" curl-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' curl-stall
expect_re "curl: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "$PATH" curl-whole whole.bin)" '^rc=0 ' curl-whole
expect_re "wget: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' wget-stall
expect_re "wget: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-whole whole.bin)" '^rc=0 ' wget-whole
expect_re "control: curl without the bounds hangs until the timeout" \
  "$(fetch_in "${WORK}/curl-unbounded/scripts/lib/common.sh" "$PATH" curl-mutant stall.bin 25)" '^rc=124 ' curl-mutant
expect_re "control: wget without the bounds hangs until the timeout" \
  "$(fetch_in "${WORK}/wget-unbounded/scripts/lib/common.sh" "${WORK}/nocurl" wget-mutant stall.bin 25)" '^rc=124 ' wget-mutant
for z in 0 00 x; do
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  if FETCH_STALL=$z bash -c '. "$1" && _get "$2" "$3"' _ "$LIB" "${URL}/whole.bin" "${WORK}/zero" > "${WORK}/zero.log" 2>&1; then
    fail "FETCH_STALL=${z} was accepted"
  else
    expect_re "FETCH_STALL=${z} is refused" "$(cat "${WORK}/zero.log")" 'FETCH_STALL must be'
  fi
done

if [ "$FAILS" -gt 0 ]; then echo "${FAILS} check(s) failed"; exit 1; fi
echo "All fetch checks passed."
