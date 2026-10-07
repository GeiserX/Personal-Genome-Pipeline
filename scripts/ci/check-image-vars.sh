#!/usr/bin/env bash
# check-image-vars.sh: fail when a script under scripts/ expands a *_IMAGE
# variable that neither versions.env nor the script itself assigns.
#
# Under `set -u` such a script dies at run time with "X_IMAGE: unbound
# variable". ShellCheck alone did not see it: SC2154 skips uppercase names
# unless check-unassigned-uppercase is enabled (see .shellcheckrc).
#
# Usage:
#   scripts/ci/check-image-vars.sh              check this repository
#   scripts/ci/check-image-vars.sh <repo_root>  check another tree
#   scripts/ci/check-image-vars.sh --self-test  prove the check can fail
#
# A reference counts as defined when versions.env sets it, when the same
# script assigns it (NAME=...), or when it carries a default (${NAME:-x}).
# Text in single quotes and in comments is neither a use nor an assignment.
set -euo pipefail

# Names versions.env defines. Sourced in a clean shell under `set -u`, so a
# variable exported by the caller cannot stand in for a missing line.
defined_names() {
  # shellcheck disable=SC2016  # $1 expands in the inner shell
  env -i bash --noprofile --norc -c \
    'set -euo pipefail; . "$1"; compgen -v' _ "$1/versions.env" \
    | grep -E '_IMAGE$' || true
}

# scan FILE: print "R<TAB>NAME<TAB>LINE" for every *_IMAGE expansion the shell
# would perform, and "A<TAB>NAME<TAB>LINE" for every NAME= assignment. It
# follows shell quoting: nothing inside single quotes counts (a container's
# own `bash -c '...'` body, a python -c string), comments do not count, a
# quoted here-document body is skipped and an unquoted one only expands, and
# $( ... ) opens a fresh quoting context. A `\$X_IMAGE` is not an expansion.
scan() {
  awk -v q="'" '
    function ref(s, i,    rest, tok, after, name) {
      rest = substr(s, i)
      if (!match(rest, /^[$][{]?[A-Z][A-Z0-9_]*_IMAGE/)) return 1
      tok = substr(rest, 1, RLENGTH)
      after = substr(rest, RLENGTH + 1)
      if (after ~ /^[A-Za-z0-9_]/) return RLENGTH                          # a longer name (FOO_IMAGES)
      if (substr(tok, 2, 1) == "{" && after ~ /^:?[-=+]/) return RLENGTH   # ${X_IMAGE:-default}
      name = tok
      sub(/^[$][{]?/, "", name)
      print "R\t" name "\t" NR
      return RLENGTH
    }
    function heredoc(s, j,    rest, strip, quoted, d) {
      strip = 0
      if (substr(s, j, 1) == "-") { strip = 1; j++ }
      while (substr(s, j, 1) ~ /[ \t]/) j++
      rest = substr(s, j)
      if (match(rest, "^[" q "\"][A-Za-z_][A-Za-z0-9_]*[" q "\"]")) {
        quoted = 1; d = substr(rest, 2, RLENGTH - 2)
      } else if (match(rest, /^\\?[A-Za-z_][A-Za-z0-9_]*/)) {
        d = substr(rest, 1, RLENGTH); quoted = (substr(d, 1, 1) == "\\"); sub(/^\\/, "", d)
      } else {
        return j
      }
      nh++; hdelim[nh] = d; hquoted[nh] = quoted; hstrip[nh] = strip
      return j + RLENGTH
    }
    BEGIN { sq = 0; ansi = 0; dq = 0; depth = 0; nh = 0; inhd = 0 }
    inhd {
      cmp = $0
      if (hstrip[1]) sub(/^\t+/, "", cmp)
      if (cmp == hdelim[1]) {
        for (k = 1; k < nh; k++) { hdelim[k] = hdelim[k + 1]; hquoted[k] = hquoted[k + 1]; hstrip[k] = hstrip[k + 1] }
        nh--
        inhd = (nh > 0)
        next
      }
      if (!hquoted[1]) {
        i = 1; n = length($0)
        while (i <= n) {
          c = substr($0, i, 1)
          if (c == "\\") { i += 2; continue }
          if (c == "$") { i += ref($0, i); continue }
          i++
        }
      }
      next
    }
    {
      line = $0; n = length(line); i = 1; prev = " "
      while (i <= n) {
        c = substr(line, i, 1)
        if (sq) {
          if (ansi && c == "\\") { i += 2; continue }
          if (c == q) { sq = 0; ansi = 0; prev = c }
          i++
          continue
        }
        if (c == "\\") { i += 2; prev = "x"; continue }
        if (c == q && !dq) { sq = 1; i++; continue }
        if (c == "\"") { dq = !dq; i++; prev = c; continue }
        if (c == "$") {
          if (substr(line, i + 1, 1) == q && !dq) { sq = 1; ansi = 1; i += 2; continue }
          if (substr(line, i, 3) == "$((") { if (depth) parens[depth] += 2; i += 3; prev = "("; continue }
          if (substr(line, i, 2) == "$(") {
            depth++; saved[depth] = dq; parens[depth] = 0; dq = 0; i += 2; prev = "("; continue
          }
          i += ref(line, i); prev = "x"; continue
        }
        if (!dq) {
          if (c == "#" && prev ~ /[ \t;|&()]/) break
          if (depth && c == "(") parens[depth]++
          if (depth && c == ")") {
            if (parens[depth] > 0) parens[depth]--
            else { dq = saved[depth]; depth--; i++; prev = c; continue }
          }
          if (substr(line, i, 2) == "<<" && substr(line, i, 3) != "<<<") {
            i = heredoc(line, i + 2); prev = "x"; continue
          }
          if (prev !~ /[A-Za-z0-9_$]/ && match(substr(line, i), /^[A-Za-z_][A-Za-z0-9_]*=/)) {
            print "A\t" substr(line, i, RLENGTH - 1) "\t" NR
            i += RLENGTH; prev = "="; continue
          }
        }
        prev = c; i++
      }
      if (nh > 0) inhd = 1
    }
  ' "$1"
}

check() {
  local root=$1 defined rel kind name lineno found assigned scripts=0 refs=0 bad=0
  local -a files
  if [ ! -f "${root}/versions.env" ]; then
    echo "ERROR: ${root}/versions.env not found" >&2
    return 2
  fi
  defined=$(defined_names "$root")

  while IFS= read -r f; do files+=("$f"); done < <(
    find "${root}/scripts" -name '*.sh' -not -path '*/scripts/ci/*' | LC_ALL=C sort)
  if [ "${#files[@]}" -eq 0 ]; then
    echo "ERROR: no scripts found under ${root}/scripts" >&2
    return 2
  fi

  local report=""
  for f in "${files[@]}"; do
    scripts=$((scripts + 1))
    rel=${f#"${root}"/}
    found=$(scan "$f")
    # Names this script assigns itself (NAME=..., export NAME=..., || NAME=...),
    # read from a variable: no pipe, so SIGPIPE cannot drop a match.
    assigned=$(awk -F'\t' '$1 == "A" { print $2 }' <<<"$found")
    while IFS=$'\t' read -r kind name lineno; do
      [ "$kind" = R ] || continue
      refs=$((refs + 1))
      if grep -qxF "$name" <<<"$defined"; then continue; fi
      if grep -qxF "$name" <<<"$assigned"; then continue; fi
      report+="${name}"$'\t'"${rel}:${lineno}"$'\n'
      bad=$((bad + 1))
    done <<<"$found"
  done

  if [ "$bad" -gt 0 ]; then
    echo "FAIL: ${bad} use(s) of an image variable that versions.env does not define:"
    printf '%s' "$report" | awk -F'\t' '
      !($1 in seen) { seen[$1] = 1; order[++n] = $1 }
      { hits[$1] = hits[$1] "    " $2 "\n" }
      END { for (i = 1; i <= n; i++) printf "  %s\n%s", order[i], hits[order[i]] }'
    echo "Define each name in versions.env, or assign it in the script that uses it."
    return 1
  fi
  echo "OK: ${refs} image variable uses in ${scripts} scripts are all defined."
}

self_test() {
  local out rc fail=0
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  mkdir -p "${tmp}/bad/scripts" "${tmp}/good/scripts"
  printf 'FOO_IMAGE="example/foo:1.0"\n' > "${tmp}/bad/versions.env"
  # shellcheck disable=SC2016  # the dollar signs are the test input
  printf '%s\n' \
    '#!/usr/bin/env bash' \
    'docker run --rm "${FOO_IMAGE}" true' \
    'docker run --rm "${NOPE_IMAGE}" true' \
    'LOCAL_IMAGE="example/local:1.0"; docker run --rm "$LOCAL_IMAGE" true' \
    'docker run --rm "${DEFAULTED_IMAGE:-example/d:1.0}" true' \
    'echo "${MISSING_IMAGES[@]}"' \
    '# docker run --rm "$COMMENTED_IMAGE" true' \
    "printf '%s\\n' '\$QUOTED_IMAGE'" \
    "docker run --rm \"\${FOO_IMAGE}\" bash -c 'echo \$INNER_IMAGE'" \
    'docker run --rm "$TRAILING_IMAGE" true # TRAILING_IMAGE=example/t:1.0' \
    "cat <<'EOF'" \
    "it's \$QUOTED_HD_IMAGE" \
    'EOF' \
    'docker run --rm "${AFTER_HD_IMAGE}" true' \
    > "${tmp}/bad/scripts/a.sh"
  # A script far bigger than a pipe buffer that assigns its own image on line
  # 2: the assignment lookup must not lose to SIGPIPE under pipefail.
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'BIG_IMAGE="example/big:1.0"' \
    'docker run --rm "${BIG_IMAGE}" true' > "${tmp}/bad/scripts/big.sh"
  awk 'BEGIN { for (i = 0; i < 200000; i++) print ": filler line " i }' >> "${tmp}/bad/scripts/big.sh"
  cp "${tmp}/bad/versions.env" "${tmp}/good/versions.env"
  cp "${tmp}/bad/scripts/big.sh" "${tmp}/good/scripts/big.sh"
  grep -vE 'NOPE_IMAGE|TRAILING_IMAGE|AFTER_HD_IMAGE' "${tmp}/bad/scripts/a.sh" > "${tmp}/good/scripts/a.sh"

  # Must be reported: an undefined use, a use whose only "assignment" is in a
  # comment, and a use after a quoted here-document holding an apostrophe.
  rc=0; out=$(check "${tmp}/bad" 2>&1) || rc=$?
  [ "$rc" -eq 1 ] || { echo "self-test: planted tree exited ${rc}, expected 1: FAIL"; fail=1; }
  for want in 'NOPE_IMAGE scripts/a.sh:3' 'TRAILING_IMAGE scripts/a.sh:10' 'AFTER_HD_IMAGE scripts/a.sh:14'; do
    n=${want%% *}
    if grep -qx "  ${n}" <<<"$out" && grep -qx "    ${want#* }" <<<"$out"; then
      echo "self-test: planted ${n} is reported at ${want#* }: PASS"
    else
      echo "self-test: planted ${n} was NOT reported at ${want#* }: FAIL"
      fail=1
    fi
  done
  # Must not be reported: defined, defaulted, a longer name, a comment, text
  # in single quotes or in a quoted here-document.
  for n in FOO_IMAGE LOCAL_IMAGE DEFAULTED_IMAGE MISSING_IMAGES COMMENTED_IMAGE BIG_IMAGE \
           QUOTED_IMAGE INNER_IMAGE QUOTED_HD_IMAGE; do
    if grep -q "^  ${n}\$" <<<"$out"; then
      echo "self-test: ${n} reported although it is defined or not a use: FAIL"
      fail=1
    fi
  done
  [ "$fail" -eq 0 ] || printf '%s\n' "$out"

  rc=0; out=$(check "${tmp}/good" 2>&1) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "self-test: clean tree passes (exit 0): PASS"
  else
    echo "self-test: clean tree failed (exit ${rc}): FAIL"
    printf '%s\n' "$out"
    fail=1
  fi
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -h|--help) sed -n '2,15p' "$0" ;;
  *) check "${1:-$(cd "$(dirname "$0")/../.." && pwd)}" ;;
esac
