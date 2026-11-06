#!/usr/bin/env bash
# docs/not-reported.md says what the pipeline leaves out and what short reads
# cannot see. Every item there is a claim, so every top-level bullet must cite
# its source: a PMID, a DOI or a URL.
#
# Checks:
#   1. every bullet of the page cites a source, and the page has at least
#      MIN_BULLETS bullets (an empty or reshaped page would pass otherwise);
#   2. a planted copy with one bullet that cites nothing fails;
#   3. a planted copy whose only "source" is a link to another page of the
#      docs fails (a relative link is not a source).
#
# Run: bash tests/test_not_reported_sources.sh
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
PAGE="${REPO}/docs/not-reported.md"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
MIN_BULLETS=25

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# check_page FILE: prints each bullet without a source; returns 1 when there
# is one, or when the page has fewer than MIN_BULLETS bullets.
check_page() {
  awk -v min="$MIN_BULLETS" '
    /^```/ { fence = !fence; next }
    fence { next }
    /^- / {
      n++
      if ($0 !~ /PMID [0-9]+/ && $0 !~ /(DOI|doi)[: ]+10\.[0-9]+\// && $0 !~ /https?:\/\//) {
        printf "line %d has no PMID, DOI or URL: %s\n", NR, substr($0, 1, 100)
        bad++
      }
    }
    END {
      if (n < min) { printf "%d bullets, expected at least %d\n", n, min; bad++ }
      exit bad ? 1 : 0
    }' "$1"
}

if out=$(check_page "$PAGE"); then
  pass "docs/not-reported.md: $(grep -c '^- ' "$PAGE") bullets, each with a PMID, DOI or URL"
else
  fail "docs/not-reported.md:"
  printf '%s\n' "$out"
fi

cp "$PAGE" "${WORK}/page.md"
printf '\n- **A planted item.** Short reads cannot see it, trust us.\n' >> "${WORK}/page.md"
if out=$(check_page "${WORK}/page.md"); then
  fail "control: a planted bullet with no source passed"
else
  pass "control: a planted bullet with no source fails: ${out}"
fi

cp "$PAGE" "${WORK}/page.md"
printf '\n- **Another planted item.** See [paralogs](35-paralogs.md).\n' >> "${WORK}/page.md"
if out=$(check_page "${WORK}/page.md"); then
  fail "control: a planted bullet citing only a docs page passed"
else
  pass "control: a planted bullet citing only a docs page fails: ${out}"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "RESULT: ${FAILS} FAILED"
  exit 1
fi
echo "RESULT: ALL PASS"
