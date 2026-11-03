#!/usr/bin/env bash
# check-pgs-labels.sh: compare the disease label of every PGS ID in
# assets/pgs_scores.tsv (the list step 25 and the PRS process score) with the
# trait the PGS Catalog reports for that ID.
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
# The score list also names, per score, the evaluations whose effect the
# reports print: ppm_id, ancestry and effect, three columns of entries
# separated by ";" and in the same order (PPM000591;PPM000594 with EUR;AFR
# and HR 1.53;HR 1.27). Each PPM id is asked of
# https://www.pgscatalog.org/rest/performance/<PPM>: it must exist, belong to
# that score, have been measured in samples of that ancestry (the catalog's
# ancestry_broad, mapped to the 1000 Genomes group pgsc_calc names), and
# report that odds ratio (OR) or hazard ratio (HR). A wrong PPM id, ancestry
# or effect is a failure, as a wrong label is.
#
# Usage:
#   scripts/ci/check-pgs-labels.sh [FILE]   check FILE (default assets/pgs_scores.tsv)
#   scripts/ci/check-pgs-labels.sh --self-test
#
# Reads the score list (a "PGS000018<TAB>Coronary artery disease" line per
# score; a third column, applies_to, is not read) and the two map formats step
# 25 used before it: "PGS000018|Coronary artery disease" entries and the older
# ["coronary_artery_disease"]="PGS000018".
# Prints a Markdown table on stdout (and one of the evaluations, when the file
# has them). Exit 0 when every label and evaluation matches, 1 when a label or
# an evaluation is wrong, an ID is unknown, the API answers empty or not at all, or the
# file holds no PGS ID (an empty check is a failure, never a pass).
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# Prints "ID<TAB>label" for every entry of the map in $1.
read_map() {
  local f=$1
  {
    # Current format, assets/pgs_scores.tsv: "PGS000018<TAB>Coronary artery disease"
    awk -F'\t' '$1 ~ /^PGS[0-9]{6}$/ && NF >= 2 && $2 != "" { print $1 "\t" $2 }' "$f"
    # Step 25's older map: "PGS000018|Coronary artery disease"
    sed -nE 's/^[[:space:]]*"(PGS[0-9]{6})\|([^"]*)".*/\1\t\2/p' "$f"
    # Older format: ["coronary_artery_disease"]="PGS000018"
    sed -nE 's/^[[:space:]]*\["([^"]+)"\]="(PGS[0-9]{6})".*/\2\t\1/p' "$f"
  }
}

# Prints the trait_reported of PGS ID $1, or nothing (with exit 1) when the
# API fails, times out, or answers without that ID. The API answers {} with
# HTTP 200 for an unknown ID, so the id field is checked too. PGS_API is read
# on every call, so a caller can point a single check at another host.
trait_of() {
  local id=$1 body api=${PGS_API:-https://www.pgscatalog.org/rest/score}
  body=$(curl -fsS --connect-timeout 10 --max-time 30 --retry "${PGS_RETRY:-3}" --retry-delay 5 --retry-all-errors \
    -H 'Accept: application/json' "${api}/${id}") || return 1
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

# Prints "ID<TAB>ppm_id<TAB>ancestry<TAB>effect" for every row of a score
# list ($1) with a ppm_id, the columns found by the header's names.
read_metrics() {
  awk -F'\t' '
    !h && $1 == "pgs_id" { for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next }
    h && c["ppm_id"] && $1 ~ /^PGS[0-9]{6}$/ && $c["ppm_id"] != "" {
      print $1 "\t" $c["ppm_id"] "\t" (c["ancestry"] ? $c["ancestry"] : "") "\t" (c["effect"] ? $c["effect"] : "")
    }' "$1"
}

# Exit 0 when evaluation $2 of the PGS Catalog belongs to score $1, was
# measured in samples of ancestry $3 and reports effect $4 ("OR 1.61").
# Prints what the catalog says, or why it does not match, on stdout. An
# API that fails or answers without that id is exit 2.
metric_matches() {
  local id=$1 ppm=$2 anc=$3 eff=$4 body api=${PGS_PERF_API:-https://www.pgscatalog.org/rest/performance}
  body=$(curl -fsS --connect-timeout 10 --max-time 30 --retry "${PGS_RETRY:-3}" --retry-delay 5 --retry-all-errors \
    -H 'Accept: application/json' "${api}/${ppm}") || { echo "(no answer)"; return 2; }
  ID="$id" PPM="$ppm" ANC="$anc" EFF="$eff" python3 -c '
import json, os, re, sys
# The PGS Catalog ancestry_broad terms of each 1000 Genomes group pgsc_calc
# reports (its ancestry_distribution uses the same codes).
GROUPS = {
    "EUR": {"European"},
    "AFR": {"African unspecified", "African American or Afro-Caribbean", "Sub-Saharan African"},
    "AMR": {"Hispanic or Latin American"},
    "EAS": {"East Asian"},
    "SAS": {"South Asian"},
}
def group(term):
    for g, terms in GROUPS.items():
        if term in terms:
            return g
    return None
try:
    d = json.loads(sys.stdin.read())
except ValueError:
    print("(not JSON)"); sys.exit(2)
if not isinstance(d, dict) or d.get("id") != os.environ["PPM"]:
    print("(no such evaluation)"); sys.exit(2)
samples = (d.get("sampleset") or {}).get("samples") or []
terms = [s.get("ancestry_broad") or "" for s in samples]
found = {group(t) for t in terms}
effects = [(e.get("name_short"), e.get("estimate")) for e in (d.get("performance_metrics") or {}).get("effect_sizes") or []]
said = "%s, %s, %s" % (d.get("associated_pgs_id"), "; ".join(terms) or "no samples",
                       "; ".join("%s %s" % e for e in effects) or "no effect size")
print(said)
m = re.fullmatch(r"(OR|HR) ([0-9]+(?:\.[0-9]+)?)", os.environ["EFF"])
ok = (d.get("associated_pgs_id") == os.environ["ID"]
      and os.environ["ANC"] in GROUPS and found == {os.environ["ANC"]}
      and m is not None
      and any(n == m.group(1) and isinstance(v, (int, float)) and abs(v - float(m.group(2))) < 1e-9 for n, v in effects))
sys.exit(0 if ok else 1)
' <<<"$body"
}

# The evaluations table of file $1; exit 1 when one is wrong. A score's
# three lists must have the same number of entries with none empty, and an
# ancestry may appear once per score.
check_metrics() {
  local f=$1 bad=0 n=0 id ppms ancs effs i said rc
  local -a P A E
  [ -n "$(read_metrics "$f")" ] || return 0
  echo ""
  echo "| PGS ID | PPM ID | Ancestry | Effect | PGS Catalog (score, ancestry, effect) | Result |"
  echo "|---|---|---|---|---|---|"
  while IFS=$'\t' read -r id ppms ancs effs; do
    # read -a drops a trailing empty entry, so "PPM000591;" would pass as one.
    if [[ "${ppms};" == *';;'* || "${ancs};" == *';;'* || "${effs};" == *';;'* \
          || $ppms == ';'* || $ancs == ';'* || $effs == ';'* ]]; then
      echo "| ${id} | ${ppms} | ${ancs} | ${effs} | | ERROR: an empty entry (a ; at the start or end, or two in a row) |"
      bad=$((bad + 1)); continue
    fi
    IFS=';' read -r -a P <<<"$ppms"
    IFS=';' read -r -a A <<<"$ancs"
    IFS=';' read -r -a E <<<"$effs"
    if [ "${#P[@]}" -ne "${#A[@]}" ] || [ "${#P[@]}" -ne "${#E[@]}" ]; then
      echo "| ${id} | ${ppms} | ${ancs} | ${effs} | | ERROR: ppm_id, ancestry and effect have different numbers of entries |"
      bad=$((bad + 1)); continue
    fi
    if [ "$(printf '%s\n' "${A[@]}" | sort | uniq -d)" != "" ]; then
      echo "| ${id} | ${ppms} | ${ancs} | ${effs} | | ERROR: an ancestry is listed twice |"
      bad=$((bad + 1)); continue
    fi
    for i in "${!P[@]}"; do
      n=$((n + 1)); rc=0
      said=$(metric_matches "$id" "${P[$i]}" "${A[$i]}" "${E[$i]}") || rc=$?
      case $rc in
        0) echo "| ${id} | ${P[$i]} | ${A[$i]} | ${E[$i]} | ${said} | ok |" ;;
        1) echo "| ${id} | ${P[$i]} | ${A[$i]} | ${E[$i]} | ${said} | MISMATCH |"; bad=$((bad + 1)) ;;
        *) echo "| ${id} | ${P[$i]} | ${A[$i]} | ${E[$i]} | ${said} | ERROR: lookup failed or empty |"; bad=$((bad + 1)) ;;
      esac
    done
  done < <(read_metrics "$f")
  echo ""
  echo "Checked ${n} evaluations, ${bad} problem(s)."
  [ "$bad" -eq 0 ]
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
  check_metrics "$f" || bad=$((bad + 1))
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

  # Control 1b: the same in the score list format (assets/pgs_scores.tsv).
  printf 'pgs_id\ttrait_reported\nPGS000018\tCoronary artery disease\nPGS000020\tInflammatory bowel disease\n' > "${d}/list.tsv"
  rc=0; out=$(check_file "${d}/list.tsv" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q '| PGS000020 .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000018 .*| ok |' <<<"$out"; then
    echo "SELF-TEST FAIL: the score list format was not read or not judged:"; echo "$out"; fail=1
  fi

  # Control 1c: with the applies_to column after the label, the label is
  # still the second column (a wrong one fails, a right one passes).
  printf 'pgs_id\ttrait_reported\tapplies_to\nPGS000662\tProstate cancer\tmale\nPGS000020\tInflammatory bowel disease\t\n' > "${d}/sex.tsv"
  rc=0; out=$(check_file "${d}/sex.tsv" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q '| PGS000020 | Inflammatory bowel disease | .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000662 | Prostate cancer | Prostate cancer | ok |' <<<"$out"; then
    echo "SELF-TEST FAIL: the label was not read as the second column beside applies_to:"; echo "$out"; fail=1
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

  # Control 4: an API that does not answer is an error, not a pass. The file
  # holds only a right label, so the failure can only come from the host
  # (192.0.2.1 is a documentation address that never answers).
  printf '%s\n' '  "PGS000018|Coronary artery disease"' > "${d}/right.sh"
  rc=0; out=$(PGS_API=https://192.0.2.1/rest/score PGS_RETRY=0 check_file "${d}/right.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ] || ! grep -q '| PGS000018 | Coronary artery disease | (no answer) | ERROR: lookup failed or empty |' <<<"$out"; then
    echo "SELF-TEST FAIL: an unreachable API passed:"; echo "$out"; fail=1
  fi

  # Control 5: a file with no PGS ID is an error.
  echo '# nothing here' > "${d}/empty.sh"
  rc=0; out=$(check_file "${d}/empty.sh" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "SELF-TEST FAIL: a file with no PGS ID passed:"; echo "$out"; fail=1
  fi

  # Control 6: the evaluations. Right ones pass (two of one score, in
  # order); a PPM id of another score, a wrong ancestry, a wrong effect, an
  # unknown PPM id, lists of different lengths and lists ending in ";" each
  # fail. The labels are
  # all right, so the failure can only come from the evaluations.
  printf 'pgs_id\ttrait_reported\tapplies_to\tppm_id\tancestry\teffect\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n' \
    $'PGS000018\tCoronary artery disease\t\tPPM000591;PPM000594\tEUR;AFR\tHR 1.53;HR 1.27' \
    $'PGS000004\tBreast cancer\tfemale\tPPM000591\tEUR\tHR 1.53' \
    $'PGS000016\tAtrial fibrillation\t\tPPM014988\tAFR\tHR 1.4' \
    $'PGS000334\tLate-onset Alzheimer’s disease\t\tPPM022981\tEUR\tHR 1.99' \
    $'PGS000662\tProstate cancer\tmale\tPPM999999\tEUR\tOR 2.04' \
    $'PGS000014\tType 2 diabetes (T2D)\t\tPPM021437;PPM021589\tEUR\tOR 1.71' \
    $'PGS000017\tInflammatory bowel disease\t\tPPM000591;\tEUR;\tHR 1.53;' > "${d}/metrics.tsv"
  rc=0; out=$(check_file "${d}/metrics.tsv" 2>&1) || rc=$?
  echo "Control 6, the planted evaluations (three must be a mismatch, three an error, two must pass):"
  grep '^| PGS0' <<<"$out" | grep 'PPM' || true
  if [ "$rc" -eq 0 ] || ! grep -q 'Checked 7 PGS IDs, 0 problem(s).' <<<"$out" \
     || ! grep -q '| PGS000018 | PPM000591 | EUR | HR 1.53 | .*| ok |' <<<"$out" \
     || ! grep -q '| PGS000018 | PPM000594 | AFR | HR 1.27 | .*| ok |' <<<"$out" \
     || ! grep -q '| PGS000004 | PPM000591 | EUR | HR 1.53 | .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000016 | PPM014988 | AFR | HR 1.4 | .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000334 | PPM022981 | EUR | HR 1.99 | .*| MISMATCH |' <<<"$out" \
     || ! grep -q '| PGS000662 | PPM999999 | EUR | OR 2.04 | .*| ERROR: lookup failed or empty |' <<<"$out" \
     || ! grep -q '| PGS000014 | .*| ERROR: ppm_id, ancestry and effect have different numbers of entries |' <<<"$out" \
     || ! grep -q '| PGS000017 | .*| ERROR: an empty entry (a ; at the start or end, or two in a row) |' <<<"$out" \
     || ! grep -q 'Checked 6 evaluations, 6 problem(s).' <<<"$out"; then
    echo "SELF-TEST FAIL: an evaluation was judged wrongly:"; echo "$out"; fail=1
  fi

  if [ "$fail" -eq 0 ]; then
    echo "SELF-TEST OK: wrong label (list and map formats), old format, empty answer, no answer, empty map, and a wrong PPM id, ancestry or effect or an empty entry are all caught."
  fi
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -h|--help) awk 'NR > 1 && /^set -euo/ {exit} NR > 1' "$0" ;;
  *) check_file "${1:-${ROOT}/assets/pgs_scores.tsv}" ;;
esac
