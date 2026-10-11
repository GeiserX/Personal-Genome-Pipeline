#!/usr/bin/env bash
# scripts/ci/check-links.sh --self-test asks github.com for a page that
# exists (200) and one that does not (404). When github.com answers 5xx for a
# while, the self-test must try again, and if the host still answers 5xx it
# must fail with a message that names the outage. A 5xx is never a pass.
#
# A fake `curl` on PATH answers instead of the network: 503 for the first
# FAKE_CURL_5XX_CALLS requests of each URL, then the real answer (404 for the
# missing page, 200 otherwise, no answer for 192.0.2.1).
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/curl" <<'FAKE'
#!/usr/bin/env bash
url=${*: -1}
key=$(printf '%s' "$url" | tr -c 'A-Za-z0-9' _)
n=$(( $(cat "${FAKE_CURL_STATE}/${key}" 2>/dev/null || echo 0) + 1 ))
echo "$n" > "${FAKE_CURL_STATE}/${key}"
printf '%s\n' "$url" >> "${FAKE_CURL_STATE}/calls"
case "$url" in
  https://192.0.2.1/*) printf '000'; exit 28 ;;
esac
if [ "$n" -le "${FAKE_CURL_5XX_CALLS:-0}" ]; then printf '503'; exit 0; fi
case "$url" in
  *no-such-file*) printf '404' ;;
  *) printf '200' ;;
esac
FAKE
chmod +x "${WORK}/bin/curl"

GOOD=https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/LICENSE

# run_self_test NAME 5XX_CALLS: the self-test with the fake curl; sets RC and OUT
run_self_test() {
  local state="${WORK}/$1"
  mkdir -p "$state"
  RC=0
  OUT=$(PATH="${WORK}/bin:${PATH}" FAKE_CURL_STATE="$state" FAKE_CURL_5XX_CALLS="$2" \
    LINKS_5XX_TRIES=3 LINKS_5XX_DELAY=0 bash "${REPO}/scripts/ci/check-links.sh" --self-test 2>&1) || RC=$?
  CALLS=$(grep -cxF "$GOOD" "${state}/calls" || true)
}

# 503 on the first try (HEAD and GET), then the real answers: a pass, after a retry.
run_self_test flaky 2
if [ "$RC" -eq 0 ] && grep -qF 'SELF-TEST OK' <<<"$OUT"; then
  pass "503 then 200: the self-test passes"
else
  echo "$OUT"; fail "503 then 200: the self-test exited ${RC}"
fi
if [ "$CALLS" -ge 3 ]; then pass "the 200 URL was asked again after its 503 (${CALLS} requests)"
else fail "the 200 URL was not asked again after its 503 (${CALLS} requests)"; fi

# 503 on every try: a failure that names the outage, after a bounded number of tries.
run_self_test down 1000
if [ "$RC" -ne 0 ]; then pass "a persistent 503 fails the self-test"; else fail "a persistent 503 passed the self-test"; fi
if grep -qE 'SELF-TEST FAIL: github\.com answered 503 .*on all 3 tries.*outage at github\.com' <<<"$OUT"; then
  pass "the failure names the github.com outage"
else
  echo "$OUT"; fail "the failure does not name the github.com outage"
fi
if [ "$CALLS" -eq 6 ]; then pass "3 tries, HEAD and GET each (${CALLS} requests)"
else fail "expected 6 requests of the 200 URL (3 tries, HEAD and GET), got ${CALLS}"; fi
if grep -qF 'SELF-TEST OK' <<<"$OUT"; then fail "a persistent 503 printed SELF-TEST OK"; fi

[ "$FAILS" -eq 0 ]
