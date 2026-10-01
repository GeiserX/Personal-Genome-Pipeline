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
set -euo pipefail

# Names versions.env defines. Sourced in a clean shell under `set -u`, so a
# variable exported by the caller cannot stand in for a missing line.
defined_names() {
  # shellcheck disable=SC2016  # $1 expands in the inner shell
  env -i bash --noprofile --norc -c \
    'set -euo pipefail; . "$1"; compgen -v' _ "$1/versions.env" \
    | grep -E '_IMAGE$' || true
}

# Print "NAME<TAB>LINE" for every *_IMAGE expansion in one script.
references() {
  awk '
    /^[[:space:]]*#/ { next }
    {
      line = $0
      while (match(line, /[$][{]?[A-Z][A-Z0-9_]*_IMAGE/)) {
        tok  = substr(line, RSTART, RLENGTH)
        pre  = (RSTART > 1) ? substr(line, RSTART - 1, 1) : ""
        line = substr(line, RSTART + RLENGTH)
        if (line ~ /^[A-Za-z0-9_]/) continue          # a longer name (FOO_IMAGES)
        if (pre == "\\") continue                    # \$X_IMAGE: expanded later, elsewhere
        braced = (substr(tok, 2, 1) == "{")
        if (braced && line ~ /^:?[-=+]/) continue    # ${X_IMAGE:-default}
        name = tok
        sub(/^[$][{]?/, "", name)
        print name "\t" NR
      }
    }
  ' "$1"
}

check() {
  local root=$1 defined rel name lineno scripts=0 refs=0 bad=0
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
    while IFS=$'\t' read -r name lineno; do
      [ -n "$name" ] || continue
      refs=$((refs + 1))
      if grep -qxF "$name" <<<"$defined"; then continue; fi
      # Assigned in this script (NAME=..., export NAME=..., || NAME=...)?
      if grep -v '^[[:space:]]*#' "$f" | grep -qE "(^|[^A-Za-z0-9_\$])${name}="; then continue; fi
      report+="${name}"$'\t'"${rel}:${lineno}"$'\n'
      bad=$((bad + 1))
    done < <(references "$f")
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
    > "${tmp}/bad/scripts/a.sh"
  # A script far bigger than a pipe buffer that assigns its own image on line
  # 2: the assignment lookup must not lose to SIGPIPE under pipefail.
  # shellcheck disable=SC2016
  printf '%s\n' '#!/usr/bin/env bash' 'BIG_IMAGE="example/big:1.0"' \
    'docker run --rm "${BIG_IMAGE}" true' > "${tmp}/bad/scripts/big.sh"
  awk 'BEGIN { for (i = 0; i < 200000; i++) print ": filler line " i }' >> "${tmp}/bad/scripts/big.sh"
  cp "${tmp}/bad/versions.env" "${tmp}/good/versions.env"
  cp "${tmp}/bad/scripts/big.sh" "${tmp}/good/scripts/big.sh"
  grep -v NOPE_IMAGE "${tmp}/bad/scripts/a.sh" > "${tmp}/good/scripts/a.sh"

  rc=0; out=$(check "${tmp}/bad" 2>&1) || rc=$?
  if [ "$rc" -eq 1 ] && grep -q 'NOPE_IMAGE' <<<"$out" && grep -q 'scripts/a.sh:3' <<<"$out"; then
    echo "self-test: planted \${NOPE_IMAGE} is reported (exit 1): PASS"
  else
    echo "self-test: planted \${NOPE_IMAGE} was NOT reported (exit ${rc}): FAIL"
    fail=1
  fi
  for n in FOO_IMAGE LOCAL_IMAGE DEFAULTED_IMAGE MISSING_IMAGES COMMENTED_IMAGE BIG_IMAGE; do
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
