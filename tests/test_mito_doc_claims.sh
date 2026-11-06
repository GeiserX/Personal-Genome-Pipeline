#!/usr/bin/env bash
# docs/12-mito-haplogroup.md must not claim that a mitochondrial haplogroup
# says anything about longevity, lifespan or metabolic traits. Those links come
# from small case-control studies that population structure in mtDNA misleads,
# and the pipeline reports no health reading for a haplogroup
# (docs/not-reported.md).
#
# Checks:
#   1. the page exists, is about the haplogroup, and has none of the claims;
#   2. a planted copy with a longevity claim fails;
#   3. a planted copy with a metabolic-trait claim fails.
#
# Run: bash tests/test_mito_doc_claims.sh
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
PAGE="${REPO}/docs/12-mito-haplogroup.md"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

CLAIMS='longevity|long-lived|lifespan|life span|metabolic trait'

# check_page FILE: prints each line with a claim; returns 1 when there is one,
# or when the file is not the haplogroup page at all.
check_page() {
  local f=$1
  if ! grep -q 'haplogroup' "$f" 2>/dev/null; then
    echo "${f}: missing, or not the haplogroup page"
    return 1
  fi
  if grep -niE "$CLAIMS" "$f"; then
    return 1
  fi
  return 0
}

if out=$(check_page "$PAGE"); then
  pass "docs/12-mito-haplogroup.md makes no longevity, lifespan or metabolic claim"
else
  fail "docs/12-mito-haplogroup.md: ${out}"
fi

cp "$PAGE" "${WORK}/page.md"
printf '\nSome haplogroups have known population-level health associations (e.g., longevity).\n' >> "${WORK}/page.md"
if out=$(check_page "${WORK}/page.md"); then
  fail "control: a planted longevity claim passed"
else
  pass "control: a planted longevity claim fails: ${out}"
fi

cp "$PAGE" "${WORK}/page.md"
printf '\nHaplogroup J is linked to metabolic traits.\n' >> "${WORK}/page.md"
if out=$(check_page "${WORK}/page.md"); then
  fail "control: a planted metabolic-trait claim passed"
else
  pass "control: a planted metabolic-trait claim fails: ${out}"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "RESULT: ${FAILS} FAILED"
  exit 1
fi
echo "RESULT: ALL PASS"
