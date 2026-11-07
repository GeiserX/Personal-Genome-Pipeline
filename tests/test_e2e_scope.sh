#!/usr/bin/env bash
# scripts/ci/e2e-scope.sh picks the e2e cases of a pull request, and the
# required e2e check trusts it: a case it leaves out never runs. Each check
# here is run against a planted fault first: a case with no paths row, a
# second row, a row for no case, a base row with globs, a glob that matches
# no file, a glob that matches a file outside the e2e scope. Then the
# picks: nothing for a path outside the scope, every case for shared wiring, the smoke cases plus the matching ones for a helper,
# the case whose file or row changed, the cases those need, never a base
# case. Last, the real table must pass the check.
set -euo pipefail
export LC_ALL=C

REPO=$(cd "$(dirname "$0")/.." && pwd)
S="${REPO}/scripts/ci/e2e-scope.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
ok()  { echo "ok:   $*"; }
bad() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
TAB=$'\t'

# Two base cases, the two smoke cases (the second needs the first), a case
# for one helper, a case that needs it, and a case for one module.
export E2E_CASE_DIR="${WORK}/cases" E2E_PATHS="${WORK}/paths.tsv" E2E_TIMES="${WORK}/times.tsv"
mkdir -p "$E2E_CASE_DIR"
CASES="01-base 02-base nextflow-from-fastq-1-bash nextflow-from-fastq-2-nextflow traits-1 traits-2 lpa-1"
for c in $CASES lib; do : > "${E2E_CASE_DIR}/${c}.sh"; done
printf '%s\n' "01-base${TAB}base${TAB}10${TAB}-" "02-base${TAB}base${TAB}10${TAB}-" \
  "lpa-1${TAB}shard${TAB}10${TAB}-" \
  "nextflow-from-fastq-1-bash${TAB}shard${TAB}10${TAB}-" \
  "nextflow-from-fastq-2-nextflow${TAB}shard${TAB}10${TAB}nextflow-from-fastq-1-bash" \
  "traits-1${TAB}shard${TAB}10${TAB}-" "traits-2${TAB}shard${TAB}10${TAB}traits-1" > "$E2E_TIMES"
paths() {
  printf '# a comment\n\n'
  printf '%s\n' "01-base${TAB}-" "02-base${TAB}-" "lpa-1${TAB}modules/local/lpa/* scripts/46-lpa.sh" \
    "nextflow-from-fastq-1-bash${TAB}scripts/01b-fastp-qc.sh" "nextflow-from-fastq-2-nextflow${TAB}-" \
    "traits-1${TAB}-" "traits-2${TAB}bin/traits_report.py"
}
paths > "$E2E_PATHS"

# --- check: the planted faults ------------------------------------------------
if out=$("$S" check); then ok "check passes the good table: ${out}"
else bad "check failed the good table: ${out}"; fi
broken() {  # broken DESC WANT FILTER [ROW...]: the table through FILTER plus ROWs must fail check, naming WANT
  local desc=$1 want=$2 filter=$3; shift 3
  { paths | grep -v -- "$filter"; [ "$#" -gt 0 ] && printf '%s\n' "$@"; } > "${WORK}/broken.tsv"
  if out=$(E2E_PATHS="${WORK}/broken.tsv" "$S" check 2>&1); then bad "check accepted ${desc}"
  elif grep -q -- "$want" <<< "$out"; then ok "check refuses ${desc}: $(grep -- "$want" <<< "$out")"
  else bad "check refused ${desc} without naming it: ${out}"; fi
}
broken "a case with no row" "case traits-1 has no row" "^traits-1${TAB}"
broken "a second row for a case" "traits-2 has more than one row" "^NONE" "traits-2${TAB}main.nf"
broken "a row for a case that is not a file" "row for gone" "^NONE" "gone${TAB}main.nf"
broken "a base row with globs" "base case 01-base" "^01-base${TAB}" "01-base${TAB}main.nf"
broken "a glob that matches no file" "bin/no_such_helper.py matches no tracked file" "^traits-2${TAB}" \
  "traits-2${TAB}bin/traits_report.py bin/no_such_helper.py"
broken "a module glob that matches no file" "modules/local/nope/\\* matches no tracked file" "^lpa-1${TAB}" \
  "lpa-1${TAB}modules/local/nope/*"
broken "a row of one column" "has 1 columns" "^traits-2${TAB}" "traits-2"
broken "a glob that matches a file outside the e2e scope" "docs/testing.md matches docs/testing.md, which is outside the e2e scope" \
  "^traits-2${TAB}" "traits-2${TAB}bin/traits_report.py docs/testing.md"
broken "a wildcard glob that reaches outside the e2e scope" "docs/\\* matches docs/" "^traits-2${TAB}" \
  "traits-2${TAB}bin/traits_report.py docs/*"
mkdir -p "${WORK}/nosmoke"
cp "${E2E_CASE_DIR}"/*.sh "${WORK}/nosmoke/"
rm "${WORK}/nosmoke/nextflow-from-fastq-2-nextflow.sh"
paths | grep -v "^nextflow-from-fastq-2-nextflow${TAB}" > "${WORK}/nosmoke.tsv"
if out=$(E2E_CASE_DIR="${WORK}/nosmoke" E2E_PATHS="${WORK}/nosmoke.tsv" "$S" check 2>&1); then bad "check accepted a smoke case with no file"
elif grep -q "smoke case nextflow-from-fastq-2-nextflow is not a file" <<< "$out"; then ok "check refuses a smoke case with no file"
else bad "check refused a missing smoke case without naming it: ${out}"; fi

# --- select -------------------------------------------------------------------
pick() {  # pick CHANGED-PATHS... : the picked cases on one line
  printf '%s\n' "$@" > "${WORK}/changed.txt"
  "$S" select "${WORK}/changed.txt" "${ROWS:-/dev/null}" 2> "${WORK}/why.txt" | paste -sd' ' -
}
expect() {  # expect DESC WANT CHANGED-PATHS...
  local desc=$1 want=$2 got; shift 2
  got=$(pick "$@")
  if [ "$got" = "$want" ]; then ok "${desc}: ${got:-nothing}"
  else bad "${desc}: got '${got}', want '${want}' ($(paste -sd' ' "${WORK}/why.txt"))"; fi
}
SMOKE="nextflow-from-fastq-1-bash nextflow-from-fastq-2-nextflow"
expect "a docs-only change runs nothing" "" docs/testing.md README.md
expect "a change outside the scope runs nothing" "" tests/test_cpic_parser.py .github/workflows/lint.yml
expect "main.nf runs every shard case" "lpa-1 ${SMOKE} traits-1 traits-2" main.nf docs/testing.md
expect "a workflow file runs every shard case" "lpa-1 ${SMOKE} traits-1 traits-2" workflows/pgx.nf
expect "the e2e runner runs every shard case" "lpa-1 ${SMOKE} traits-1 traits-2" tests/e2e/lib.sh
expect "one helper runs the smoke cases, its case and the case that needs it" "${SMOKE} traits-1 traits-2" bin/traits_report.py
expect "a file in one module runs that module's case" "lpa-1 ${SMOKE}" modules/local/lpa/main.nf
expect "an in-scope path no case runs runs the smoke cases" "${SMOKE}" scripts/ci/check-name.sh
expect "a committed fixture a case reads is in the scope" "${SMOKE}" tests/fixtures/consumer_export/sites.tsv.gz
expect "a fixture only the unit tests read runs nothing" "" tests/fixtures/pharmcat/report-3.4.0.json
expect "a case's own file runs it" "lpa-1 ${SMOKE}" tests/e2e/lpa-1.sh
expect "a base step runs the smoke cases only (the base leg runs every base case)" "${SMOKE}" tests/e2e/01-base.sh
printf '%s\n' traits-2 > "${WORK}/rows.txt"
ROWS="${WORK}/rows.txt" expect "a changed table row runs its case and what it needs" "${SMOKE} traits-1 traits-2" docs/x.md
expect "a deleted case file runs the smoke cases" "${SMOKE}" tests/e2e/gone-1.sh
if grep -q 'traits-1: needed by traits-2' "${WORK}/why.txt" 2> /dev/null || { pick bin/traits_report.py > /dev/null; grep -q 'traits-2: runs bin/traits_report.py' "${WORK}/why.txt"; }; then
  ok "select says why it picked each case"
else bad "select gave no reasons: $(cat "${WORK}/why.txt")"; fi
paths | grep -v "^traits-1${TAB}" > "${WORK}/norow.tsv"
printf '%s\n' bin/traits_report.py > "${WORK}/changed.txt"
if E2E_PATHS="${WORK}/norow.tsv" "$S" select "${WORK}/changed.txt" > /dev/null 2>&1; then bad "select ran with a case that has no row"
else ok "select refuses a table with a case that has no row"; fi

# --- rows: the cases of changed table rows ------------------------------------
cat > "${WORK}/diff.txt" << EOF
diff --git a/scripts/ci/e2e-case-paths.tsv b/scripts/ci/e2e-case-paths.tsv
--- a/scripts/ci/e2e-case-paths.tsv
+++ b/scripts/ci/e2e-case-paths.tsv
@@ -3 +3,2 @@
-# an old comment
+traits-2${TAB}bin/traits_report.py
+lpa-1${TAB}-
diff --git a/scripts/ci/e2e-case-times.tsv b/scripts/ci/e2e-case-times.tsv
--- a/scripts/ci/e2e-case-times.tsv
+++ b/scripts/ci/e2e-case-times.tsv
@@ -9 +9 @@
-traits-2${TAB}shard${TAB}5${TAB}-
+traits-2${TAB}shard${TAB}10${TAB}traits-1
EOF
if [ "$("$S" rows < "${WORK}/diff.txt" | paste -sd' ' -)" = "traits-2 lpa-1" ]; then ok "rows names each changed row's case once, never a comment or a file header"
else bad "rows gave: $("$S" rows < "${WORK}/diff.txt" | paste -sd' ' -)"; fi

# --- the real table and the workflow -------------------------------------------
unset E2E_CASE_DIR E2E_PATHS E2E_TIMES
if out=$("$S" check 2>&1); then ok "real table: ${out}"
else bad "real table: ${out}"; fi
have=$(awk '/^  push:/ {p = 1; next} p && /^  [a-z_]+:/ {exit}
            p && /^      - / {sub(/^      - /, ""); gsub(/\047/, ""); print}' "${REPO}/.github/workflows/e2e.yml" | sort)
if [ "$("$S" scope)" = "$have" ]; then ok "the push paths of e2e.yml are the scope"
else bad "push paths and scope differ: $(diff <("$S" scope) <(echo "$have") || true)"; fi
printf '%s\n' bin/traits_report.py > "${WORK}/changed.txt"
"$S" select "${WORK}/changed.txt" > "${WORK}/real.txt" 2> /dev/null
if [ -s "${WORK}/real.txt" ] && E2E_SELECT="${WORK}/real.txt" "${REPO}/scripts/ci/e2e-shards.sh" plan 2 > /dev/null; then
  ok "a real pick splits: $(paste -sd' ' "${WORK}/real.txt")"
else bad "a real pick did not split: $(paste -sd' ' "${WORK}/real.txt")"; fi

printf '%s\n' tests/fixtures/consumer_export/sites.tsv.gz > "${WORK}/changed.txt"
if "$S" select "${WORK}/changed.txt" 2> /dev/null | grep -qx consumer-export-1; then ok "the real table runs consumer-export-1 when its fixture changes"
else bad "the real table does not run consumer-export-1 when its fixture changes"; fi

[ "$FAILS" -eq 0 ]
