#!/usr/bin/env bash
# e2e-shards.sh — split the e2e cases of tests/e2e over parallel jobs.
#
# Usage: scripts/ci/e2e-shards.sh base
#        scripts/ci/e2e-shards.sh plan N
#        scripts/ci/e2e-shards.sh cover FILE
#        scripts/ci/e2e-shards.sh merge EXPECTED RESULTS...
#        scripts/ci/e2e-shards.sh times TRACE...
#
#   base     prints the base cases, one per line, in name order. The job
#            "e2e base" runs them once, and every shard starts from what
#            they wrote.
#   plan N   prints "SHARD<TAB>CASE" for every other case, shards 1 to N, in
#            name order within a shard. A case shares a shard with every case
#            it needs (the needs column of e2e-case-times.tsv, followed
#            through), so a group of cases that read each other's outputs is
#            never split. Groups go to the least loaded shard, heaviest first.
#   cover FILE  exits 1 unless FILE lists every case exactly once (it names
#            each duplicate, missing and unknown case).
#   merge EXPECTED RESULTS...  prints one markdown table of every case of
#            EXPECTED (a list of case names) from the result files the jobs
#            wrote ("CASE<TAB>RESULT<TAB>SECONDS<TAB>JOB" lines), and exits 1
#            unless each case has exactly one row and that row is a pass.
#   times TRACE...  prints a new e2e-case-times.tsv from the trace files of a
#            traced run ("CASE<TAB>SECONDS<TAB>WRITERS<TAB>READS" lines,
#            scripts/ci/e2e-run.sh with E2E_TRACE=true). Keeps each case's
#            stage and needs, adds the traced needs. Exits 1 when the trace
#            saw no case read another case's files, which means it saw nothing.
#
# Env: E2E_CASE_DIR  the case directory (default tests/e2e)
#      E2E_TIMES     the table (default scripts/ci/e2e-case-times.tsv)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
CASE_DIR=${E2E_CASE_DIR:-${REPO}/tests/e2e}
TIMES=${E2E_TIMES:-${REPO}/scripts/ci/e2e-case-times.tsv}
DEFAULT_SECONDS=300
export LC_ALL=C

die() { echo "ERROR: $*" >&2; exit 1; }

# cases: every case name, in the order e2e-run.sh runs them.
cases() {
  local f
  for f in "$CASE_DIR"/*.sh; do
    [ -f "$f" ] || continue
    f=${f##*/}
    [ "$f" = lib.sh ] || echo "${f%.sh}"
  done | sort
}

# rows: the table without comments and blank lines.
rows() { grep -vE '^[[:space:]]*(#|$)' "$TIMES" || true; }

# check_table: a malformed, duplicate or stale row, or a need that cannot be
# met, stops the split instead of silently moving a case.
check_table() {
  local bad
  bad=$(rows | awk -F'\t' -v list="$(cases | paste -sd' ' -)" '
    BEGIN { n = split(list, c, " "); for (i = 1; i <= n; i++) { known[c[i]] = 1; pos[c[i]] = i } }
    {
      if (NF != 4) { print "row " NR " has " NF " columns, not 4: " $0; next }
      if (seen[$1]++) print "case " $1 " has more than one row"
      if (!($1 in known)) print "row for " $1 ", which is not a file in the case directory"
      if ($2 != "base" && $2 != "shard") print "case " $1 ": stage is " $2 ", not base or shard"
      if ($3 !~ /^[0-9]+$/) print "case " $1 ": seconds is " $3 ", not a whole number"
      stage[$1] = $2; needs[$1] = $4
    }
    END {
      for (k in needs) {
        if (needs[k] == "-") continue
        if (stage[k] == "base") { print "base case " k " needs " needs[k] ": a base case needs nothing"; continue }
        m = split(needs[k], d, ",")
        for (j = 1; j <= m; j++) {
          if (!(d[j] in known)) print "case " k " needs " d[j] ", which is not a case"
          else if (stage[d[j]] == "base") print "case " k " needs base case " d[j] ": leave base cases out of needs"
          else if (pos[d[j]] >= pos[k]) print "case " k " needs " d[j] ", which does not run before it"
        }
      }
    }')
  [ -z "$bad" ] || die "$(printf '%s\n' "in ${TIMES}:" "$bad")"
}

base() {
  check_table
  awk -F'\t' 'FILENAME == ARGV[1] { if ($2 == "base") b[$1] = 1; next } ($1 in b)' <(rows) <(cases)
}

plan() {
  local n=$1
  [[ "$n" =~ ^[1-9][0-9]*$ ]] || die "plan: N must be a positive whole number, got '${n}'"
  check_table
  # Union-find over needs: one line per group, "WEIGHT<TAB>FIRST<TAB>MEMBERS".
  awk -F'\t' -v def="$DEFAULT_SECONDS" '
    function find(x) { while (up[x] != x) x = up[x]; return x }
    FILENAME == ARGV[1] { stage[$1] = $2; secs[$1] = $3; needs[$1] = $4; next }
    {
      c = $1
      if (stage[c] == "base") next
      order[++n] = c; up[c] = c
      w[c] = (c in secs) ? secs[c] : def
    }
    END {
      for (i = 1; i <= n; i++) {
        c = order[i]
        if (!(c in needs) || needs[c] == "-") continue
        m = split(needs[c], d, ",")
        for (j = 1; j <= m; j++) { a = find(c); b = find(d[j]); if (a != b) up[a] = b }
      }
      for (i = 1; i <= n; i++) {
        c = order[i]; r = find(c)
        weight[r] += w[c]
        members[r] = (r in members) ? members[r] " " c : c
        if (!(r in first)) first[r] = c
      }
      for (r in members) printf "%d\t%s\t%s\n", weight[r], first[r], members[r]
    }' <(rows) <(cases) \
  | sort -t "$(printf '\t')" -k1,1nr -k2,2 \
  | awk -F'\t' -v n="$n" '
    BEGIN { for (s = 1; s <= n; s++) load[s] = 0 }
    {
      best = 1
      for (s = 2; s <= n; s++) if (load[s] < load[best]) best = s
      load[best] += $1
      m = split($3, c, " ")
      for (j = 1; j <= m; j++) printf "%d\t%s\n", best, c[j]
    }' \
  | sort -t "$(printf '\t')" -k1,1n -k2,2
}

cover() {
  local file=$1 dup missing extra rc=0
  [ -f "$file" ] || die "cover: no file ${file}"
  dup=$(sort "$file" | uniq -d)
  missing=$(comm -23 <(cases) <(sort -u "$file"))
  extra=$(comm -13 <(cases) <(sort -u "$file"))
  if [ -n "$dup" ]; then echo "Listed more than once: $(paste -sd' ' - <<< "$dup")"; rc=1; fi
  if [ -n "$missing" ]; then echo "Not listed: $(paste -sd' ' - <<< "$missing")"; rc=1; fi
  if [ -n "$extra" ]; then echo "Not a case: $(paste -sd' ' - <<< "$extra")"; rc=1; fi
  [ "$rc" -eq 0 ] && echo "Every one of the $(cases | wc -l | tr -d ' ') cases is listed once."
  return "$rc"
}

merge() {
  local expected=$1; shift
  [ -f "$expected" ] || die "merge: no file ${expected}"
  awk -F'\t' '
    FILENAME == ARGV[1] { if (NF) want[++n] = $1; next }
    NF >= 3 {
      rows[$1]++
      res[$1] = $2; secs[$1] = $3; job[$1] = ($4 == "" ? "?" : $4)
    }
    END {
      for (i = 1; i <= n; i++) wanted[want[i]] = 1
      bad = 0
      print "| Case | Result | Time | Job |"
      print "|---|---|---|---|"
      for (i = 1; i <= n; i++) {
        c = want[i]
        if (!(c in rows)) { r = "FAIL (did not run)"; t = "-"; j = "-"; bad++ }
        else if (rows[c] > 1) { r = "FAIL (ran " rows[c] " times)"; t = secs[c] "s"; j = job[c]; bad++ }
        else { r = res[c]; t = secs[c] "s"; j = job[c]; if (r != "pass") bad++ }
        printf "| %s | %s | %s | %s |\n", c, r, t, j
      }
      for (c in rows) if (!(c in wanted)) { printf "| %s | FAIL (not a planned case) | %ss | %s |\n", c, secs[c], job[c]; bad++ }
      print ""
      printf "%d of %d cases failed or did not run once.\n", bad, n
      exit (bad > 0)
    }' "$expected" "$@"
}

times() {
  check_table
  [ "$#" -gt 0 ] || die "times: no trace file given"
  # Header: the comment lines of the current table, kept as they are.
  grep -E '^[[:space:]]*#' "$TIMES" || true
  awk -F'\t' -v def="$DEFAULT_SECONDS" '
    FILENAME == ARGV[1] { order[++n] = $1; next }
    FILENAME == ARGV[2] { stage[$1] = $2; secs[$1] = $3; needs[$1] = $4; next }
    {
      if (NF < 4) next
      tsecs[$1] = $2; reads += $4
      if ($3 != "-" && $3 != "") { tw[$1] = $3; edges++ }
    }
    END {
      if (reads == 0 || edges == 0) {
        print "ERROR: the trace saw " reads + 0 " reads of earlier files and " edges + 0 " cases reading another case: atime tracing did not work" > "/dev/stderr"
        exit 1
      }
      for (i = 1; i <= n; i++) {
        c = order[i]
        st = (c in stage) ? stage[c] : "shard"
        s = (c in tsecs) ? tsecs[c] : ((c in secs) ? secs[c] : def)
        delete have
        out = ""
        if (st != "base") {
          k = 0
          if ((c in needs) && needs[c] != "-") k = split(needs[c], d, ",")
          if (c in tw) { m = split(tw[c], e, ","); for (j = 1; j <= m; j++) d[++k] = e[j] }
          for (j = 1; j <= k; j++) {
            x = d[j]
            if (x == "" || x == c || (x in have) || stage[x] == "base" || x == "setup") continue
            have[x] = 1; out = (out == "" ? x : out "," x)
          }
        }
        printf "%s\t%s\t%d\t%s\n", c, st, s, (out == "" ? "-" : out)
      }
    }' <(cases) <(rows) "$@"
}

cmd=${1:-}
[ "$#" -gt 0 ] && shift
case "$cmd" in
  base)  base ;;
  plan)  plan "${1:?Usage: $0 plan N}" ;;
  cover) cover "${1:?Usage: $0 cover FILE}" ;;
  merge) [ "$#" -ge 1 ] || die "Usage: $0 merge EXPECTED RESULTS..."; merge "$@" ;;
  times) times "$@" ;;
  *) sed -n '2,32p' "$0" >&2; exit 2 ;;
esac
