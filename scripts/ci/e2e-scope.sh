#!/usr/bin/env bash
# e2e-scope.sh — pick the e2e cases a pull request can affect.
#
# Usage: scripts/ci/e2e-scope.sh check
#        scripts/ci/e2e-scope.sh scope
#        scripts/ci/e2e-scope.sh rows < DIFF
#        scripts/ci/e2e-scope.sh select CHANGED [ROWS]
#
#   check    exits 1 unless every case of tests/e2e has exactly one row in
#            e2e-case-paths.tsv, every row names a case, a base case's row is
#            "-" and every glob matches a tracked file, and only files inside
#            the e2e scope. It names each fault.
#   scope    prints the e2e scope as the push paths of e2e.yml are written
#            (a directory as DIR/**, then the single files), sorted.
#   rows     reads `git diff -U0` output of the two case tables
#            (e2e-case-paths.tsv, e2e-case-times.tsv) on stdin and prints the
#            case of every added or removed row, once each.
#   select CHANGED [ROWS]  prints the shard cases to run, one per line, in
#            name order, from CHANGED (the paths the pull request changed, one
#            per line) and ROWS (the cases whose table rows it changed, from
#            `rows`). Nothing when no path is in the e2e scope. Every shard
#            case when a path is shared wiring (FULL below). Otherwise the
#            smoke cases, every case whose own file or table row changed,
#            every case with a glob that matches a changed path, and the
#            cases those need (the needs column of e2e-case-times.tsv,
#            followed through). Base cases are never printed: whenever a case
#            runs, the base leg runs all of them. The reason for each pick
#            goes to stderr.
#
# A glob is matched as in bash [[ PATH == GLOB ]]: * matches any characters,
# slashes too, so modules/local/lpa/* covers the whole module.
#
# Env: E2E_CASE_DIR  the case directory (default tests/e2e)
#      E2E_PATHS     the paths table (default scripts/ci/e2e-case-paths.tsv)
#      E2E_TIMES     the times table (default scripts/ci/e2e-case-times.tsv)
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
CASE_DIR=${E2E_CASE_DIR:-${REPO}/tests/e2e}
PATHS=${E2E_PATHS:-${REPO}/scripts/ci/e2e-case-paths.tsv}
TIMES=${E2E_TIMES:-${REPO}/scripts/ci/e2e-case-times.tsv}
export LC_ALL=C

# The paths the e2e cases exercise. A pull request that changes none of them
# runs no case. The push paths of e2e.yml are the same list (e2e-plan checks).
SCOPE_DIRS=(scripts modules workflows bin conf assets tests/e2e .github/actions)
SCOPE_FILES=(tests/fixtures/VERSION main.nf nextflow.config nextflow_schema.json versions.env .github/workflows/e2e.yml
  # The committed fixtures the cases read (their rows name them too).
  tests/fixtures/annotsv/class45.sh tests/fixtures/comphet/make_fixture.py tests/fixtures/consumer_export/sites.tsv.gz
  tests/fixtures/families/stub_families.tsv tests/fixtures/pharmcat/refblocks.g.vcf tests/fixtures/sv/cnvpytor.vcf.in
  tests/fixtures/sv/delly.vcf.in tests/fixtures/sv/manta.vcf.in tests/fixtures/vcf/combined_allele_clinvar.vcf
  tests/fixtures/vcf/combined_allele_sample.vcf tests/fixtures/vcf/grch37_no_contigs.vcf)

# Shared wiring: every case runs through these, so a change to one runs every
# case. The pipeline's entry and config, the shared library and launcher, the
# step table, the alignment and calling steps whose BAM and VCF every case
# starts from, the fixture, and the e2e runner and planner themselves.
FULL=(main.nf nextflow.config nextflow_schema.json versions.env assets/steps.tsv
  'workflows/*' 'conf/*' 'scripts/lib/*' scripts/run-all.sh
  scripts/02-alignment.sh scripts/03-deepvariant.sh
  tests/fixtures/VERSION scripts/ci/build-fixture.sh
  scripts/ci/e2e-run.sh scripts/ci/e2e-shards.sh scripts/ci/e2e-scope.sh
  tests/e2e/lib.sh 'tests/e2e/bin/*' .github/workflows/e2e.yml '.github/actions/*')

# Run on every pull request that runs any case: the bash steps and the
# Nextflow pipeline from the fixture's FASTQ, and the parity check of the two.
SMOKE=(nextflow-from-fastq-1-bash nextflow-from-fastq-2-nextflow)

die() { echo "ERROR: $*" >&2; exit 1; }

cases() {
  local f
  for f in "$CASE_DIR"/*.sh; do
    [ -f "$f" ] || continue
    f=${f##*/}
    [ "$f" = lib.sh ] || echo "${f%.sh}"
  done | sort
}
rows_of() { grep -vE '^[[:space:]]*(#|$)' "$1" || true; }
base_cases() { rows_of "$TIMES" | awk -F'\t' '$2 == "base" {print $1}' | sort; }

tracked() {
  if git -C "$REPO" rev-parse --is-inside-work-tree > /dev/null 2>&1; then
    git -C "$REPO" ls-files
  else
    (cd "$REPO" && find . -type f ! -path './.git/*' | sed 's|^\./||')
  fi
}

check() {
  local bad="" all base c globs g f hit gl out
  [ -f "$PATHS" ] || die "no table ${PATHS}"
  all=$(cases); base=$(base_cases)
  for c in "${SMOKE[@]}"; do
    grep -qx -- "$c" <<< "$all" || bad+="smoke case ${c} is not a file in the case directory"$'\n'
  done
  bad+=$(rows_of "$PATHS" | awk -F'\t' -v list="$(paste -sd' ' - <<< "$all")" -v blist="$(paste -sd' ' - <<< "$base")" '
    BEGIN { n = split(list, c, " "); for (i = 1; i <= n; i++) known[c[i]] = 1
            m = split(blist, b, " "); for (i = 1; i <= m; i++) isbase[b[i]] = 1 }
    {
      if (NF != 2 || $2 == "") { print "row " NR " has " NF " columns, not 2: " $0; next }
      if (seen[$1]++) print "case " $1 " has more than one row"
      if (!($1 in known)) print "row for " $1 ", which is not a file in the case directory"
      if (($1 in isbase) && $2 != "-") print "base case " $1 ": its row must be -, the base leg runs every base case"
    }
    END { for (k in known) if (!(k in seen)) print "case " k " has no row: add one with the globs of the files it runs" }' | sort)
  local -A have=()
  local scope_globs=() p
  for p in "${SCOPE_DIRS[@]}"; do scope_globs+=("${p}/*"); done
  scope_globs+=("${SCOPE_FILES[@]}")
  mapfile -t files < <(tracked)
  for f in "${files[@]}"; do have[$f]=1; done
  while IFS=$'\t' read -r c globs; do
    [ "$globs" = "-" ] && continue
    read -ra gl <<< "$globs"
    for g in "${gl[@]}"; do
      hit=0
      if [[ $g != *[*?[]* ]]; then
        [ -n "${have[$g]:-}" ] && hit=1
      else
        for f in "${files[@]}"; do
          # shellcheck disable=SC2053  # g is a glob on purpose
          if [[ $f == $g ]]; then hit=1; break; fi
        done
      fi
      [ "$hit" = 1 ] || bad+=$'\n'"case ${c}: ${g} matches no tracked file"
      # A file outside the e2e scope never starts a run, so a glob that
      # reaches one would let a change to it merge untested.
      out=""
      if [[ $g != *[*?[]* ]]; then
        in_list "$g" "${scope_globs[@]}" || out=$g
      else
        for f in "${files[@]}"; do
          # shellcheck disable=SC2053  # g is a glob on purpose
          if [[ $f == $g ]] && ! in_list "$f" "${scope_globs[@]}"; then out=$f; break; fi
        done
      fi
      [ -z "$out" ] || bad+=$'\n'"case ${c}: ${g} matches ${out}, which is outside the e2e scope (SCOPE_DIRS, SCOPE_FILES): a change to it would run no case"
    done
  done < <(rows_of "$PATHS" | awk -F'\t' 'NF == 2')
  bad=$(sed '/^$/d' <<< "$bad")
  [ -z "$bad" ] || die "$(printf '%s\n' "in ${PATHS}:" "$bad")"
  echo "Every one of the $(wc -l <<< "$all" | tr -d ' ') cases has one paths row, and every glob matches a file, inside the e2e scope."
}

scope() { { printf '%s/**\n' "${SCOPE_DIRS[@]}"; printf '%s\n' "${SCOPE_FILES[@]}"; } | sort; }

rows() {
  # "+case<TAB>..." or "-case<TAB>...", never the "+++"/"---" file headers.
  awk '/^(\+\+\+|---) / {next} /^[+-][^#[:space:]]/ {s = substr($0, 2); sub(/\t.*/, "", s); if (!seen[s]++) print s}'
}

# in_list PATH GLOB...: whether PATH matches one of the globs.
in_list() {
  local x=$1 g; shift
  # shellcheck disable=SC2053  # g is a glob on purpose
  for g in "$@"; do [[ $x == $g ]] && return 0; done
  return 1
}

select_cases() {
  local changed=$1 rowsf=${2:-/dev/null} p c globs g full="" inscope="" picked out gl
  [ -f "$changed" ] || die "select: no file ${changed}"
  [ -e "$rowsf" ] || die "select: no file ${rowsf}"
  check > /dev/null
  local scope_globs=()
  for p in "${SCOPE_DIRS[@]}"; do scope_globs+=("${p}/*"); done
  scope_globs+=("${SCOPE_FILES[@]}")
  while IFS= read -r p; do
    [ -n "$p" ] || continue
    in_list "$p" "${scope_globs[@]}" && inscope=1
    if in_list "$p" "${FULL[@]}"; then full=$p; break; fi
  done < "$changed"
  [ -s "$rowsf" ] && inscope=1
  if [ -z "$inscope" ]; then
    echo "No changed path is in the e2e scope: no case runs." >&2
    return 0
  fi
  local base; base=$(base_cases)
  if [ -n "$full" ]; then
    echo "${full} is shared wiring: every case runs." >&2
    comm -23 <(cases) <(echo "$base")
    return 0
  fi
  picked=$(
    for c in "${SMOKE[@]}"; do printf '%s\tsmoke\n' "$c"; done
    while IFS= read -r p; do
      [[ $p == tests/e2e/*.sh ]] && { c=${p#tests/e2e/}; printf '%s\tits file changed\n' "${c%.sh}"; }
    done < "$changed"
    while IFS= read -r c; do [ -n "$c" ] && printf '%s\tits table row changed\n' "$c"; done < "$rowsf"
    while IFS=$'\t' read -r c globs; do
      [ "$globs" = "-" ] && continue
      read -ra gl <<< "$globs"
      while IFS= read -r p; do
        [ -n "$p" ] || continue
        for g in "${gl[@]}"; do
          # shellcheck disable=SC2053  # g is a glob on purpose
          if [[ $p == $g ]]; then printf '%s\truns %s\n' "$c" "$p"; continue 3; fi
        done
      done < "$changed"
    done < <(rows_of "$PATHS")
  )
  # Keep only real cases outside the base stage, then add their needs.
  out=$(awk -F'\t' -v list="$(cases | paste -sd' ' -)" -v blist="$(paste -sd' ' - <<< "$base")" '
    BEGIN { n = split(list, cs, " "); for (i = 1; i <= n; i++) known[cs[i]] = 1
            m = split(blist, b, " "); for (i = 1; i <= m; i++) isbase[b[i]] = 1 }
    FILENAME == ARGV[1] { if ($0 !~ /^[[:space:]]*(#|$)/ && $4 != "-") needs[$1] = $4; next }
    ($1 in known) && !($1 in isbase) && !($1 in why) { why[$1] = $2; order[++k] = $1 }
    END {
      for (i = 1; i <= k; i++) {
        c = order[i]
        if (!(c in needs)) continue
        m = split(needs[c], d, ",")
        for (j = 1; j <= m; j++) if (!(d[j] in why)) { why[d[j]] = "needed by " c; order[++k] = d[j] }
      }
      for (i = 1; i <= k; i++) print order[i] "\t" why[order[i]]
    }' "$TIMES" <(printf '%s\n' "$picked") | sort)
  awk -F'\t' '{print "  " $1 ": " $2}' <<< "$out" >&2
  echo "$(grep -c . <<< "$out") of $(comm -23 <(cases) <(echo "$base") | wc -l | tr -d ' ') shard cases run, and the base cases." >&2
  cut -f1 <<< "$out"
}

cmd=${1:-}
[ "$#" -gt 0 ] && shift
case "$cmd" in
  check)  check ;;
  scope)  scope ;;
  rows)   rows ;;
  select) select_cases "${1:?Usage: $0 select CHANGED [ROWS]}" "${2:-}" ;;
  *) sed -n '2,36p' "$0" >&2; exit 2 ;;
esac
