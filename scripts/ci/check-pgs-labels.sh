#!/usr/bin/env bash
# check-pgs-labels.sh: compare the disease label of every PGS ID in
# scripts/25-prs.sh with the trait the PGS Catalog reports for that ID.
#
# Two scores once carried the wrong disease for months (PGS000020 was printed
# as inflammatory bowel disease and is type 2 diabetes; PGS000738 was printed
# as schizophrenia and is vitiligo). This check asks the PGS Catalog REST API
# (https://www.pgscatalog.org/rest/score/<ID>) for `trait_reported` and fails
# when a label names something else.
#
# A label matches when every word of it appears in trait_reported, after
# lower-casing and dropping apostrophes and punctuation. So "Type 2 diabetes
# (T2D)" and the older "type_2_diabetes" both match, and
# "inflammatory_bowel_disease" does not match "Type 2 diabetes (T2D)".
#
# Usage:
#   scripts/ci/check-pgs-labels.sh [FILE]   check FILE (default scripts/25-prs.sh)
#   scripts/ci/check-pgs-labels.sh --self-test
#
# Reads both map formats the script has used: "PGS000018|Coronary artery
# disease" entries and the older ["coronary_artery_disease"]="PGS000018".
# Prints a Markdown table on stdout. Exit 0 when every label matches, 1 when a
# label is wrong, an ID is unknown, the API answers empty or not at all, or the
# file holds no PGS ID (an empty check is a failure, never a pass).
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
API=${PGS_API:-https://www.pgscatalog.org/rest/score}

# Prints "ID<TAB>label" for every entry of the map in $1.
read_map() {
  local f=$1
  {
    # Current format: "PGS000018|Coronary artery disease"
    sed -nE 's/^[[:space:]]*"(PGS[0-9]{6})\|([^"]*)".*/\1\t\2/p' "$f"
    # Older format: ["coronary_artery_disease"]="PGS000018"
    sed -nE 's/^[[:space:]]*\["([^"]+)"\]="(PGS[0-9]{6})".*/\2\t\1/p' "$f"
  }
}

# Prints the trait_reported of PGS ID $1, or nothing (with exit 1) when the
# API fails, times out, or answers without that ID. The API answers {} with
# HTTP 200 for an unknown ID, so the id field is checked too.
trait_of() {
  local id=$1 body
  body=$(curl -fsS --connect-timeout 10 --max-time 30 --retry "${PGS_RETRY:-3}" --retry-delay 5 --retry-all-errors \
    -H 'Accept: application/json' "${API}/${id}") || return 1
  ID="$id" python3 -c '
import json, os, sys
try:
    d = json.loads(sys.stdin.read())
except ValueError:
    sys.exit(1)
if not isinstance(d, dict) or d.get("id") != os.environ["ID"]:
    sys.exit(1)
t = (d.get("trait_reported") or "").strip()
if not t:
    sys.exit(1)
print(t)
' <<<"$body"
}

# Exit 0 when every word of label $1 appears in trait $2.
label_matches() {
  python3 - "$1" "$2" <<'EOF'
import re, sys
def words(s):
    s = s.lower().replace("’", "").replace("'", "")
    return set(re.sub(r"[^a-z0-9]+", " ", s).split())
label, trait = words(sys.argv[1]), words(sys.argv[2])
sys.exit(0 if label and label <= trait else 1)
EOF
}

check_file() {
  local f=$1 n=0 bad=0 id label trait
  if [ ! -f "$f" ]; then
    echo "ERROR: ${f} not found" >&2
    return 1
  fi
  echo "| PGS ID | Label in $(basename "$f") | PGS Catalog trait_reported | Result |"
  echo "|---|---|---|---|"
  while IFS=$'\t' read -r id label; do
    n=$((n + 1))
    if ! trait=$(trait_of "$id"); then
      echo "| ${id} | ${label} | (no answer) | ERROR: lookup failed or empty |"
      bad=$((bad + 1))
    elif label_matches "$label" "$trait"; then
      echo "| ${id} | ${label} | ${trait} | ok |"
    else
      echo "| ${id} | ${label} | ${trait} | MISMATCH |"
      bad=$((bad + 1))
    fi
  done < <(read_map "$f")
  echo ""
  if [ "$n" -eq 0 ]; then
    echo "ERROR: no PGS ID found in ${f}; the map format changed or the file is wrong." >&2
    return 1
  fi
  echo "Checked ${n} PGS IDs, ${bad} problem(s)."
  [ "$bad" -eq 0 ]
}

self_test() {
  local rc out fail=0
  d=$(mktemp -d)
  trap 'rm -rf "$d"' EXIT

  # Control 1: a wrong label must fail and a right one must pass.
  printf '%s\n' 'PGS_SCORES=(' '  "PGS000018|Coronary artery disease"' \
    '  "PGS000020|Inflammatory bowel disease"' ')' > "${d}/new.sh"
  rc=0; out=$(check_file "${d}/new.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q '| PGS000020 .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000018 .*| ok |' <<<"$out"; then
    echo "SELF-TEST FAIL: a wrong label was not caught, or a right one was not passed:"; echo "$out"; fail=1
  fi

  # Control 2: the older associative-array format is read too.
  printf '%s\n' 'PGS_IDS=(' '  ["schizophrenia"]="PGS000738"' ')' > "${d}/old.sh"
  rc=0; out=$(check_file "${d}/old.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q '| PGS000738 | schizophrenia | Vitiligo | MISMATCH |' <<<"$out"; then
    echo "SELF-TEST FAIL: the older map format was not read or not judged:"; echo "$out"; fail=1
  fi

  # Control 3: an unknown ID (the API answers {} with HTTP 200) is an error.
  printf '%s\n' '  "PGS999999|Anything"' > "${d}/unknown.sh"
  rc=0; out=$(check_file "${d}/unknown.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q 'ERROR: lookup failed or empty' <<<"$out"; then
    echo "SELF-TEST FAIL: an empty API answer was not an error:"; echo "$out"; fail=1
  fi

  # Control 4: an API that does not answer is an error, not a pass.
  rc=0; out=$(PGS_API=https://192.0.2.1/rest/score PGS_RETRY=0 check_file "${d}/new.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "SELF-TEST FAIL: an unreachable API passed:"; echo "$out"; fail=1
  fi

  # Control 5: a file with no PGS ID is an error.
  echo '# nothing here' > "${d}/empty.sh"
  rc=0; out=$(check_file "${d}/empty.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "SELF-TEST FAIL: a file with no PGS ID passed:"; echo "$out"; fail=1
  fi

  if [ "$fail" -eq 0 ]; then
    echo "SELF-TEST OK: wrong label, old format, empty answer, no answer and empty map are all caught."
  fi
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -h|--help) sed -n '2,25p' "$0" ;;
  *) check_file "${1:-${ROOT}/scripts/25-prs.sh}" ;;
esac
