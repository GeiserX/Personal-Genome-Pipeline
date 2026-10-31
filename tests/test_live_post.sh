#!/usr/bin/env bash
# live_post (scripts/lib/common.sh) sends run-all.sh's step messages to the
# live site (bin/pgp_live.py). This test starts the site on loopback and checks
# that a message reaches events.jsonl, that without PGP_LIVE_URL nothing is
# sent, and that with the site down it warns once and returns 0 every time, so
# a run never stops for it. It first proves the store check can fail: a
# message to another run must not show up in this run's file.
#
# LIVE_POST_FROM=FILE takes the function from another file.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
FROM=${LIVE_POST_FROM:-${REPO}/scripts/lib/common.sh}
WORK=$(mktemp -d)
PID=""
cleanup() { [ -z "$PID" ] || kill "$PID" 2>/dev/null || true; rm -rf "$WORK"; }
trap cleanup EXIT
fail() { echo "FAIL: $*"; exit 1; }

SRC=$(sed -n '/^live_post() {/,/^}/p' "$FROM")
[ -n "$SRC" ] || fail "no live_post() in ${FROM}"
eval "$SRC"
command -v curl >/dev/null || fail "curl is needed"

free_port() { python3 -c 'import socket; s=socket.socket(); s.bind(("127.0.0.1", 0)); print(s.getsockname()[1])'; }
WEB=$(free_port) ING=$(free_port)
mkdir -p "${WORK}/genome"
python3 "${REPO}/bin/pgp_live.py" --root "${WORK}/genome" --state-dir "${WORK}/state" --port "$WEB" \
  --ingest-port "$ING" --quiet 2> "${WORK}/server.err" &
PID=$!
for _ in $(seq 50); do curl -s -o /dev/null "http://127.0.0.1:${ING}/" && break; sleep 0.1; done
EV="${WORK}/state/runs/demo1/events.jsonl"
MSG='{"event":"pgp_step","step":"14","state":"started"}'

# Without PGP_LIVE_URL nothing is sent
unset PGP_LIVE_URL PGP_LIVE_WARNED
live_post demo1 "$MSG" || fail "live_post without PGP_LIVE_URL returned non-zero"
[ ! -e "$EV" ] || fail "a message was stored without PGP_LIVE_URL"

# Control: a message to another run is not in this run's file
export PGP_LIVE_URL="http://127.0.0.1:${ING}/"
live_post other1 "$MSG"
if [ -e "$EV" ] && grep -q '"step":"14"' "$EV"; then fail "control: another run's message is in demo1's file"; fi
[ -s "${WORK}/state/runs/other1/events.jsonl" ] || fail "control: the other run's message was not stored"
echo "Control: the store check tells one run's file from another's."

# The site up: the message is stored (a URL ending in '/' works too)
live_post demo1 "$MSG"
grep -q '"event":"pgp_step","step":"14","state":"started"' "$EV" || fail "the message is not in ${EV}: $(cat "$EV" 2>/dev/null)"
echo "PASS: live_post stores the message"

# The site down: one warning, exit 0 each time
kill "$PID"; wait "$PID" 2>/dev/null || true; PID=""
unset PGP_LIVE_WARNED
rc=0; { live_post demo1 "$MSG"; live_post demo1 "$MSG"; } 2> "${WORK}/warn.txt" || rc=$?
[ "$rc" -eq 0 ] || fail "live_post returned ${rc} with the site down"
[ "$(grep -c 'WARNING: the live site' "${WORK}/warn.txt")" -eq 1 ] || fail "expected one warning: $(cat "${WORK}/warn.txt")"
echo "PASS: with the site down live_post warns once and returns 0"

# Control: a live_post that passes the failure on must be caught by the check above
eval "$(printf '%s\n' "$SRC" | sed -e 's/^live_post()/live_post_bad()/' -e 's/^  return 0$/  return 1/')"
rc=0; live_post_bad demo1 "$MSG" 2>/dev/null || rc=$?
[ "$rc" -ne 0 ] || fail "control: the planted live_post_bad returned 0, so the exit check cannot fail"
echo "Control: a live_post that returns the failure is caught."
