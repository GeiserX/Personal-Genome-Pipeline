#!/usr/bin/env bash
# scripts/ci/e2e-shards.sh splits the e2e cases over the parallel legs of
# e2e.yml, and the job e2e trusts its cover and merge checks to turn red when a
# case runs twice, never runs or fails. Each check here is run against a
# planted fault first: a duplicate and a missing case in a split, a stale or
# malformed table row, a result missing, failed or reported twice, a trace that
# saw nothing. Then the real table must split the real cases once each.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
S="${REPO}/scripts/ci/e2e-shards.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
ok()  { echo "ok:   $*"; }
bad() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
TAB=$'\t'

# A case directory of nine cases: two base cases, a chain c1 <- c2 <- c3, a
# pair d1 <- d2, and three cases on their own, one of them (z9) with no row.
export E2E_CASE_DIR="${WORK}/cases" E2E_TIMES="${WORK}/times.tsv"
mkdir -p "$E2E_CASE_DIR"
for c in 01-base 02-base c1 c2 c3 d1 d2 e1 z9 lib; do : > "${E2E_CASE_DIR}/${c}.sh"; done
table() {
  printf '# a comment\n\n'
  printf '%s\n' "01-base${TAB}base${TAB}10${TAB}-" "02-base${TAB}base${TAB}20${TAB}-" \
    "c1${TAB}shard${TAB}100${TAB}-" "c2${TAB}shard${TAB}100${TAB}c1" "c3${TAB}shard${TAB}100${TAB}c2" \
    "d1${TAB}shard${TAB}150${TAB}-" "d2${TAB}shard${TAB}100${TAB}d1" "e1${TAB}shard${TAB}50${TAB}-"
}
table > "$E2E_TIMES"

# --- base and plan ----------------------------------------------------------
if [ "$("$S" base | paste -sd' ' -)" = "01-base 02-base" ]; then ok "base lists the two base cases"
else bad "base gave: $("$S" base | paste -sd' ' -)"; fi

"$S" plan 2 > "${WORK}/plan.tsv"
# Heaviest group first, ties by name, each to the least loaded shard (the
# lower number on a tie): c (300) to 1, z9 (no row: 300) to 2, d (250) to 1,
# e1 (50) to 2.
want="1${TAB}c1 1${TAB}c2 1${TAB}c3 1${TAB}d1 1${TAB}d2 2${TAB}e1 2${TAB}z9"
if [ "$(paste -sd' ' "${WORK}/plan.tsv")" = "$want" ]; then ok "plan 2 keeps each chain in one shard, in name order, heaviest group first"
else bad "plan 2 gave: $(paste -sd' ' "${WORK}/plan.tsv")"; fi

{ "$S" base; cut -f2 "${WORK}/plan.tsv"; } > "${WORK}/all.txt"
if "$S" cover "${WORK}/all.txt" > /dev/null; then ok "base plus plan lists every case once"
else bad "cover rejected the real split: $("$S" cover "${WORK}/all.txt" || true)"; fi

if [ "$("$S" plan 1 | cut -f1 | sort -u)" = 1 ] && [ "$("$S" plan 1 | wc -l | tr -d ' ')" = 7 ]; then
  ok "plan 1 puts every shard case in shard 1"
else bad "plan 1 gave: $("$S" plan 1 | paste -sd' ' -)"; fi
if "$S" plan 0 > /dev/null 2>&1; then bad "plan 0 was accepted"; else ok "plan 0 is refused"; fi

# --- count: as many shards as the critical path asks for ---------------------
# Groups c 300, z9 300, d 250, e1 50: 900 seconds, the heaviest 300.
cnt() { E2E_SHARD_FLOOR=${2:-900} "$S" count "$1"; }
if [ "$(cnt 6)" = 1 ]; then ok "count: 900 seconds over a 900-second floor is one shard"
else bad "count 6 gave $(cnt 6)"; fi
if [ "$(cnt 6 100)" = 3 ]; then ok "count: 900 seconds over a 300-second critical path is three shards"
else bad "count 6 with floor 100 gave $(cnt 6 100)"; fi
if [ "$(cnt 2 100)" = 2 ]; then ok "count never goes over MAX"
else bad "count 2 with floor 100 gave $(cnt 2 100)"; fi
if [ "$(cnt 6 400)" = 2 ]; then ok "count adds a shard for a leftover of at least the floor (900 = 2 x 400 + 100: two)"
else bad "count 6 with floor 400 gave $(cnt 6 400)"; fi
if "$S" count 0 > /dev/null 2>&1; then bad "count 0 was accepted"; else ok "count 0 is refused"; fi

# --- E2E_SELECT: only the picked cases --------------------------------------
printf '%s\n' c1 c2 c3 e1 > "${WORK}/select.txt"
want="1${TAB}c1 1${TAB}c2 1${TAB}c3 2${TAB}e1"
if [ "$(E2E_SELECT="${WORK}/select.txt" "$S" plan 2 | paste -sd' ' -)" = "$want" ]; then ok "plan with E2E_SELECT splits the picked cases only"
else bad "plan with E2E_SELECT gave: $(E2E_SELECT="${WORK}/select.txt" "$S" plan 2 | paste -sd' ' -)"; fi
{ "$S" base; E2E_SELECT="${WORK}/select.txt" "$S" plan 2 | cut -f2; } > "${WORK}/picked.txt"
if E2E_SELECT="${WORK}/select.txt" "$S" cover "${WORK}/picked.txt" > /dev/null; then ok "cover with E2E_SELECT wants the base and picked cases"
else bad "cover with E2E_SELECT rejected the picked split"; fi
if "$S" cover "${WORK}/picked.txt" > /dev/null; then bad "cover without E2E_SELECT passed a split that leaves cases out"
else ok "cover without E2E_SELECT wants every case"; fi
printf '%s\n' c1 c3 > "${WORK}/select-hole.txt"
if out=$(E2E_SELECT="${WORK}/select-hole.txt" "$S" plan 2 2>&1); then bad "plan accepted a pick without c2, which c3 needs"
elif grep -q 'c3 but not c2' <<< "$out"; then ok "a pick without a need is refused: ${out#ERROR: }"
else bad "a missing need not named: ${out}"; fi
printf '%s\n' 01-base e1 > "${WORK}/select-base.txt"
if E2E_SELECT="${WORK}/select-base.txt" "$S" plan 2 > /dev/null 2>&1; then bad "plan accepted a pick naming a base case"
else ok "a pick naming a base case is refused"; fi
printf '%s\n' ghost > "${WORK}/select-ghost.txt"
if E2E_SELECT="${WORK}/select-ghost.txt" "$S" plan 2 > /dev/null 2>&1; then bad "plan accepted a pick naming no case"
else ok "a pick naming no case is refused"; fi

# --- cover: the planted duplicate and the planted missing case ---------------
{ cat "${WORK}/all.txt"; echo c2; } > "${WORK}/dup.txt"
if out=$("$S" cover "${WORK}/dup.txt"); then bad "a split listing c2 twice passed cover"
elif grep -q 'more than once: c2' <<< "$out"; then ok "planted duplicate reported: ${out}"
else bad "duplicate not named: ${out}"; fi
grep -vx d2 "${WORK}/all.txt" > "${WORK}/missing.txt"
if out=$("$S" cover "${WORK}/missing.txt"); then bad "a split without d2 passed cover"
elif grep -q 'Not listed: d2' <<< "$out"; then ok "planted missing case reported: ${out}"
else bad "missing case not named: ${out}"; fi
{ cat "${WORK}/all.txt"; echo ghost; } > "${WORK}/extra.txt"
if out=$("$S" cover "${WORK}/extra.txt"); then bad "a split naming a case that does not exist passed cover"
else ok "planted unknown case reported: ${out}"; fi

# --- the table: a broken row stops the split ---------------------------------
broken() {  # broken DESC ROW...: the table plus ROWs must make plan fail
  local desc=$1; shift
  { table; printf '%s\n' "$@"; } > "${WORK}/broken.tsv"
  if out=$(E2E_TIMES="${WORK}/broken.tsv" "$S" plan 2 2>&1); then bad "plan accepted ${desc}"
  else ok "plan refuses ${desc}: $(tail -n 1 <<< "$out")"; fi
}
broken "a second row for a case" "e1${TAB}shard${TAB}5${TAB}-"
broken "a row for a case that is not a file" "gone${TAB}shard${TAB}5${TAB}-"
: > "${E2E_CASE_DIR}/a0.sh"   # a case that runs first, for this row only
broken "a need on a case that runs later" "z9${TAB}shard${TAB}5${TAB}-" "a0${TAB}shard${TAB}5${TAB}z9"
rm -f "${E2E_CASE_DIR}/a0.sh"
broken "a need on a base case" "z9${TAB}shard${TAB}5${TAB}01-base"
broken "a base case with a need" "z9${TAB}base${TAB}5${TAB}c1"
broken "a row of three columns" "z9${TAB}shard${TAB}5"
broken "an unknown stage" "z9${TAB}later${TAB}5${TAB}-"

# --- check: every case has a row ---------------------------------------------
# z9 has no row in the planted table: check names it, while plan still runs it.
if out=$("$S" check 2>&1); then bad "check passed a table with no row for z9"
elif grep -qx 'z9' <<< "$out"; then ok "check names the case with no row: $(tr '\n' ' ' <<< "$out")"
else bad "check did not name z9: ${out}"; fi
table | grep -v "^e1${TAB}" > "${WORK}/two-missing.tsv"
if out=$(E2E_TIMES="${WORK}/two-missing.tsv" "$S" check 2>&1); then bad "check passed a table with no row for e1 and z9"
elif grep -qx 'e1' <<< "$out" && grep -qx 'z9' <<< "$out"; then ok "check names both cases with no row"
else bad "check did not name both e1 and z9: ${out}"; fi
if E2E_TIMES="${WORK}/no-such.tsv" "$S" check > /dev/null 2>&1; then bad "check passed with no table"
else ok "check refuses a table that is not there"; fi
mkdir -p "${WORK}/empty"
if E2E_CASE_DIR="${WORK}/empty" E2E_TIMES=/dev/null "$S" check > /dev/null 2>&1; then bad "check passed with no case and an empty table"
else ok "check refuses a case directory with no case"; fi
{ table; printf '%s\n' "z9${TAB}shard${TAB}5${TAB}-"; } > "${WORK}/whole.tsv"
if out=$(E2E_TIMES="${WORK}/whole.tsv" "$S" check 2>&1); then ok "check passes once z9 has its row: ${out}"
else bad "check failed a table with a row for every case: ${out}"; fi
{ table; printf '%s\n' "z9${TAB}shard${TAB}5${TAB}-" "z9${TAB}shard${TAB}6${TAB}-"; } > "${WORK}/twice.tsv"
if E2E_TIMES="${WORK}/twice.tsv" "$S" check > /dev/null 2>&1; then bad "check passed a table with two rows for z9"
else ok "check refuses a second row too"; fi

# --- merge: one pass per case, or red ---------------------------------------
res() { printf '%s\t%s\t%s\t%s\n' "$@"; }
while read -r c; do res "$c" pass 1 "shard 1"; done < "${WORK}/all.txt" > "${WORK}/r-all.tsv"
if out=$("$S" merge "${WORK}/all.txt" "${WORK}/r-all.tsv"); then
  if [ "$(grep -c '| pass |' <<< "$out")" = 9 ]; then ok "merge passes when every case passed once"
  else bad "merge table lacks rows: ${out}"; fi
else bad "merge failed on an all-pass run: ${out}"; fi
grep -v "^d2${TAB}" "${WORK}/r-all.tsv" > "${WORK}/r-missing.tsv"
if out=$("$S" merge "${WORK}/all.txt" "${WORK}/r-missing.tsv"); then bad "merge passed with d2 never run"
elif grep -q '| d2 | FAIL (did not run)' <<< "$out"; then ok "a case that did not run is red"
else bad "missing case not shown: ${out}"; fi
{ cat "${WORK}/r-all.tsv"; res c1 pass 1 "shard 2"; } > "${WORK}/r-twice.tsv"
if out=$("$S" merge "${WORK}/all.txt" "${WORK}/r-twice.tsv"); then bad "merge passed with c1 run twice"
elif grep -q '| c1 | FAIL (ran 2 times)' <<< "$out"; then ok "a case that ran twice is red"
else bad "twice-run case not shown: ${out}"; fi
sed "s/^e1${TAB}pass/e1${TAB}FAIL (exit 1)/" "${WORK}/r-all.tsv" > "${WORK}/r-fail.tsv"
if "$S" merge "${WORK}/all.txt" "${WORK}/r-fail.tsv" > /dev/null; then bad "merge passed with e1 failed"
else ok "a failed case is red"; fi
{ cat "${WORK}/r-all.tsv"; res stray pass 1 "shard 1"; } > "${WORK}/r-stray.tsv"
if "$S" merge "${WORK}/all.txt" "${WORK}/r-stray.tsv" > /dev/null; then bad "merge passed with a case nobody planned"
else ok "an unplanned result is red"; fi
if "$S" merge "${WORK}/all.txt" /dev/null > /dev/null; then bad "merge passed with no results at all"
else ok "no results at all is red"; fi

# --- times: the traced needs land in the table ------------------------------
printf '%s\n' "01-base${TAB}11${TAB}setup${TAB}3" "02-base${TAB}21${TAB}01-base${TAB}2" \
  "e1${TAB}55${TAB}setup,c1${TAB}4" "z9${TAB}33${TAB}02-base,d2${TAB}1" > "${WORK}/trace.tsv"
if "$S" times "${WORK}/trace.tsv" > "${WORK}/new.tsv"; then
  if grep -qx "e1${TAB}shard${TAB}55${TAB}c1" "${WORK}/new.tsv" && grep -qx "z9${TAB}shard${TAB}33${TAB}d2" "${WORK}/new.tsv" \
     && grep -qx "c3${TAB}shard${TAB}100${TAB}c2" "${WORK}/new.tsv" && grep -qx "02-base${TAB}base${TAB}21${TAB}-" "${WORK}/new.tsv" \
     && grep -qx '# a comment' "${WORK}/new.tsv" \
     && E2E_TIMES="${WORK}/new.tsv" "$S" plan 3 > /dev/null; then
    ok "times keeps stages and needs, adds traced needs, drops base and setup, and its table splits"
  else bad "times gave: $(cat "${WORK}/new.tsv")"; fi
else bad "times failed on a good trace"; fi
printf '%s\n' "c1${TAB}5${TAB}-${TAB}0" "c2${TAB}5${TAB}-${TAB}0" > "${WORK}/blind.tsv"
if "$S" times "${WORK}/blind.tsv" > /dev/null 2>&1; then bad "times accepted a trace that saw no read"
else ok "a trace that saw no read is refused"; fi
printf '%s\n' "c1${TAB}5${TAB}setup${TAB}7" "c2${TAB}5${TAB}setup${TAB}3" > "${WORK}/setup-only.tsv"
if "$S" times "${WORK}/setup-only.tsv" > /dev/null 2>&1; then bad "times accepted a trace whose reads were all of layout files"
else ok "a trace with reads of layout files only is refused"; fi

# --- the real table splits the real cases once each --------------------------
unset E2E_CASE_DIR E2E_TIMES
n=$(awk '/^  E2E_SHARDS:/ {print $2}' "${REPO}/.github/workflows/e2e.yml")
[[ "$n" =~ ^[0-9]+$ ]] || { bad "no E2E_SHARDS in e2e.yml"; n=6; }
if { "$S" base; "$S" plan "$n" | cut -f2; } > "${WORK}/real.txt" && out=$("$S" cover "${WORK}/real.txt"); then
  ok "real table, ${n} shards: ${out}"
else bad "real split: ${out:-plan failed}"; fi
if out=$("$S" check 2>&1); then ok "real table: ${out}"
else bad "real table: ${out}"; fi

[ "$FAILS" -eq 0 ]
