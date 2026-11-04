#!/usr/bin/env bash
# Every version of T1K and Parascopy named in a page under docs/ (for example
# "**T1K** v1.0.10" or "Parascopy 1.19.0") must be the version pinned in
# versions.env (T1K_IMAGE, PARASCOPY_IMAGE). A page that still names an old
# version after a bump says the wrong thing about the tool that runs.
#
# Checks:
#   1. docs/*.md agree with versions.env, and each tool is named with a
#      version at least once (a check with nothing to compare would always pass);
#   2. a planted copy of the docs with an old T1K version fails;
#   3. a planted copy with an old Parascopy version fails;
#   4. a planted copy with no versioned mention of T1K fails.
#
# Run: bash tests/test_doc_versions.sh
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# The version part of an image tag: quay.io/biocontainers/t1k:1.0.10--h5814d7d_0 -> 1.0.10
pinned() {
  local var=$1 line
  line=$(grep -E "^${var}=" "${REPO}/versions.env" | head -n 1)
  line=${line#*=}
  line=${line%%#*}
  line=${line//\"/}
  line=${line##*:}
  printf '%s\n' "${line%%--*}"
}

# check_docs DIR: prints one line per mismatch, returns 1 on any mismatch or
# on a tool that no page names with a version.
check_docs() {
  local dir=$1 rc=0 tool var want hits n h f found
  for tool in T1K Parascopy; do
    case "$tool" in
      T1K) var=T1K_IMAGE ;;
      Parascopy) var=PARASCOPY_IMAGE ;;
    esac
    want=$(pinned "$var")
    if [ -z "$want" ]; then
      echo "${var} not found in versions.env"
      rc=1
      continue
    fi
    n=0
    # A tool name, optional Markdown emphasis or code marks, an optional "v"
    # or "version ", then a dotted version. One grep over every page.
    hits=$(grep -noE "\b${tool}[*\`]*[[:space:]]+(v|version[[:space:]]+)?[0-9]+\.[0-9]+(\.[0-9]+)?" "$dir"/*.md || true)
    if [ -n "$hits" ]; then
      while IFS= read -r h; do
        n=$((n + 1))
        found=${h##*[!0-9.]}
        if [ "$found" != "$want" ]; then
          f=${h%%:*}
          h=${h#*:}
          echo "$(basename "$f"):${h%%:*}: ${tool} ${found}, versions.env pins ${want}"
          rc=1
        fi
      done <<< "$hits"
    fi
    if [ "$n" -eq 0 ]; then
      echo "no page names ${tool} with a version"
      rc=1
    fi
  done
  return "$rc"
}

if out=$(check_docs "${REPO}/docs"); then
  pass "docs/*.md name T1K $(pinned T1K_IMAGE) and Parascopy $(pinned PARASCOPY_IMAGE), as versions.env pins"
else
  fail "docs/*.md disagree with versions.env:"
  printf '%s\n' "$out"
fi

plant() {
  rm -rf "${WORK}/docs"
  cp -R "${REPO}/docs" "${WORK}/docs"
}

plant
sed -i.bak -E 's/(\*\*T1K\*\* )v[0-9.]+/\1v1.0.9/' "${WORK}/docs/08-hla-typing.md"
if out=$(check_docs "${WORK}/docs"); then
  fail "control: a planted T1K v1.0.9 in 08-hla-typing.md passed"
else
  pass "control: a planted T1K v1.0.9 fails: ${out}"
fi

plant
sed -i.bak -E 's/(\*\*Parascopy\*\* )[0-9.]+/\11.18.0/' "${WORK}/docs/35-paralogs.md"
if out=$(check_docs "${WORK}/docs"); then
  fail "control: a planted Parascopy 1.18.0 in 35-paralogs.md passed"
else
  pass "control: a planted Parascopy 1.18.0 fails: ${out}"
fi

plant
sed -i.bak -E 's/(\*\*T1K\*\*) v[0-9.]+/\1/' "${WORK}/docs/08-hla-typing.md"
if out=$(check_docs "${WORK}/docs"); then
  fail "control: docs that never name a T1K version passed"
else
  pass "control: docs that never name a T1K version fail: ${out}"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "RESULT: ${FAILS} FAILED"
  exit 1
fi
echo "RESULT: ALL PASS"
