#!/usr/bin/env bash
# changed-images.sh: which *_IMAGE variables of versions.env the Container Test
# workflow must run, one name per line, in versions.env order.
#
# Usage:
#   scripts/ci/changed-images.sh --all     every *_IMAGE variable
#   scripts/ci/changed-images.sh BASE      the variables to test for the change
#                                          from commit BASE to the working tree
#   scripts/ci/changed-images.sh --self-test
#
# A variable is listed for a change when
#   - its value differs from BASE, or it is new;
#   - a coupled data variable changed that one of its rows in
#     tests/smoke/commands.tsv names with `also=` (PYPGX_BUNDLE_VERSION is
#     tested by the PYPGX_IMAGE row, for example);
#   - or the test itself changed (tests/smoke/, image-smoke.sh, this script or
#     the workflow): then every variable is listed, so a change to a row or to
#     the runner is proved on every image before it merges.
# Values are compared after the shell reads the file, so a comment or a
# quoting change is not a change.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
HARNESS_RE='^(tests/smoke/|scripts/ci/image-smoke\.sh$|scripts/ci/changed-images\.sh$|\.github/workflows/container-test\.yml$)'

# values FILE: NAME=value for every NAME="..." line, as the shell reads them.
values() {
  # shellcheck disable=SC2016  # expanded by the inner bash
  env -i bash --noprofile --norc -c '
    set -a; . "$1" >/dev/null 2>&1; set +a
    for n in $(grep -oE "^[A-Z][A-Z0-9_]*=" "$1" | tr -d =); do printf "%s=%s\n" "$n" "${!n}"; done' _ "$1"
}

image_vars() { grep -oE '^[A-Z][A-Z0-9_]*_IMAGE=' "$1" | tr -d = | awk '!seen[$0]++'; }

changed() {
  local root=$1 base=$2 old new files var n
  git -C "$root" cat-file -e "${base}^{commit}" 2>/dev/null || { echo "ERROR: ${base} is not a commit" >&2; return 2; }
  files=$(git -C "$root" diff --name-only "$base" --)
  files+=$'\n'$(git -C "$root" ls-files --others --exclude-standard)
  if grep -Eq "$HARNESS_RE" <<< "$files"; then
    echo "the image test itself changed: every image" >&2
    image_vars "${root}/versions.env"
    return 0
  fi
  old=$(mktemp) new=$(mktemp)
  git -C "$root" show "${base}:versions.env" > "$old" 2>/dev/null || : > "$old"
  local -a diff=()
  mapfile -t diff < <(comm -13 <(values "$old" | sort) <(values "${root}/versions.env" | sort) | cut -d= -f1)
  rm -f "$old" "$new"
  [ "${#diff[@]}" -gt 0 ] && echo "changed in versions.env: ${diff[*]}" >&2
  declare -A pick=()
  for n in "${diff[@]}"; do
    if [[ "$n" == *_IMAGE ]]; then
      pick[$n]=1
    else
      # Rows that name this data variable with also=.
      while read -r var; do
        [ -n "$var" ] && pick[$var]=1
      done < <(awk -F'\t' -v n="$n" '!/^#/ && NF >= 2 {k = split($2, o, ","); for (i = 1; i <= k; i++) if (o[i] == "also=" n) print $1}' \
                 "${root}/tests/smoke/commands.tsv")
    fi
  done
  while read -r var; do
    [ -n "${pick[$var]:-}" ] && echo "$var"
  done < <(image_vars "${root}/versions.env")
  return 0
}

self_test() {
  local tmp fails=0 got
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' RETURN
  mkdir -p "${tmp}/tests/smoke" "${tmp}/scripts/ci"
  printf '%s\n' 'A_IMAGE="a:1"' 'B_IMAGE="b:1"  # a comment' 'C_IMAGE="c:1"' 'C_DATA="7"' > "${tmp}/versions.env"
  printf 'A_IMAGE\t-\t-\ttrue\ttrue\nB_IMAGE\t-\t-\ttrue\ttrue\nC_IMAGE\talso=C_DATA\t-\ttrue\ttrue\n' > "${tmp}/tests/smoke/commands.tsv"
  git -C "$tmp" init -q
  git -C "$tmp" -c core.hooksPath=/dev/null -c user.name=t -c user.email=t@t add -A
  git -C "$tmp" -c core.hooksPath=/dev/null -c user.name=t -c user.email=t@t commit -qm base
  expect() {  # expect DESCRIPTION WANT
    got=$(changed "$tmp" HEAD 2>/dev/null | paste -sd ' ' -)
    if [ "$got" = "$2" ]; then echo "[PASS] $1 (${got:-nothing})"
    else echo "[FAIL] $1 (got '${got}', want '$2')"; fails=$((fails + 1)); fi
  }
  expect "no change lists nothing" ""
  sed -i.bak 's/b:1"  # a comment/b:1"/' "${tmp}/versions.env"
  expect "a comment change lists nothing" ""
  sed -i.bak 's/b:1/b:2/' "${tmp}/versions.env"
  expect "a bumped image is listed alone" "B_IMAGE"
  sed -i.bak 's/C_DATA="7"/C_DATA="8"/' "${tmp}/versions.env"
  expect "a coupled data variable lists the image of its also= row" "B_IMAGE C_IMAGE"
  printf 'D_IMAGE="d:1"\n' >> "${tmp}/versions.env"
  expect "a new image is listed" "B_IMAGE C_IMAGE D_IMAGE"
  git -C "$tmp" checkout -q -- versions.env
  echo "# note" >> "${tmp}/tests/smoke/commands.tsv"
  expect "a change to the table lists every image" "A_IMAGE B_IMAGE C_IMAGE"
  git -C "$tmp" checkout -q -- tests/smoke/commands.tsv
  touch "${tmp}/scripts/ci/image-smoke.sh"
  expect "a new file of the runner lists every image" "A_IMAGE B_IMAGE C_IMAGE"
  if [ "$fails" -gt 0 ]; then echo "self-test: ${fails} case(s) failed" >&2; return 1; fi
  echo "self-test: all cases passed"
}

case "${1:-}" in
  --all) image_vars "${REPO}/versions.env" ;;
  --self-test) self_test ;;
  ""|-h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 2 ;;
  *) changed "$REPO" "$1" ;;
esac
