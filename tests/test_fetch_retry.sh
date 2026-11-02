#!/usr/bin/env bash
# fetch (scripts/lib/common.sh) and the downloads it retries.
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
FAILS=0
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
# expect DESC GOT WANT: GOT matches the extended regex WANT.
expect() {
  if [[ "$2" =~ $3 ]]; then pass "$1 ($2)"; else fail "$1: got '$2', want /$3/"; sed 's/^/  /' "${WORK}/${4:-none}.log" 2>/dev/null | tail -n 5; fi
}

LIB="${REPO}/scripts/lib/common.sh"
expect "curl: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "$PATH" curl-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' curl-stall
expect "curl: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "$PATH" curl-whole whole.bin)" '^rc=0 ' curl-whole
expect "wget: a stalled transfer fails both tries in seconds and keeps the .part" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-stall stall.bin)" '^rc=1 secs=([0-9]|1[0-9]|20) part=1000$' wget-stall
expect "wget: a whole file downloads under the same bounds" \
  "$(fetch_in "$LIB" "${WORK}/nocurl" wget-whole whole.bin)" '^rc=0 ' wget-whole
expect "control: curl without the bounds hangs until the timeout" \
  "$(fetch_in "${WORK}/curl-unbounded/scripts/lib/common.sh" "$PATH" curl-mutant stall.bin 25)" '^rc=124 ' curl-mutant
expect "control: wget without the bounds hangs until the timeout" \
  "$(fetch_in "${WORK}/wget-unbounded/scripts/lib/common.sh" "${WORK}/nocurl" wget-mutant stall.bin 25)" '^rc=124 ' wget-mutant
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
if FETCH_STALL=0 bash -c '. "$1" && _get "$2" "$3"' _ "$LIB" "${URL}/whole.bin" "${WORK}/zero" > "${WORK}/zero.log" 2>&1; then
  fail "FETCH_STALL=0 was accepted"
else
  expect "FETCH_STALL=0 is refused" "$(cat "${WORK}/zero.log")" 'FETCH_STALL must be'
fi

if [ "$FAILS" -gt 0 ]; then echo "${FAILS} check(s) failed"; exit 1; fi
echo "All fetch stall checks passed."
