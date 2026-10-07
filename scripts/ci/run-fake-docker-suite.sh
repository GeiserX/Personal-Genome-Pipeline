#!/usr/bin/env bash
# run-fake-docker-suite.sh: run the pipeline's shell entry points against a
# fake docker, wget, curl and df (scripts/ci/fake-docker/), one case file at a
# time from tests/fake-docker/*.sh.
#
# Usage: scripts/ci/run-fake-docker-suite.sh [case-name ...]
#        scripts/ci/run-fake-docker-suite.sh --self-test
#
# Each case runs under `env -i` in its own temp directory (see
# scripts/ci/fake-docker/lib.sh for what it receives). A case fails when it
# exits non-zero, runs longer than CASE_TIMEOUT seconds (default 300), or when
# the fake docker logged an empty image name, a missing image or an unknown
# docker option. Adding a case is adding a file; this runner never changes.
# --self-test runs five planted cases and checks the runner fails the four
# that must fail (a non-zero exit, a quoted empty image, an unquoted empty
# image variable and an unknown --name=value option, each docker error
# swallowed) and passes the clean one.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
FAKE_BIN="${ROOT}/scripts/ci/fake-docker"
CASE_DIR=${FAKE_DOCKER_CASE_DIR:-${ROOT}/tests/fake-docker}
CASE_TIMEOUT=${CASE_TIMEOUT:-300}

if [ "${1:-}" = "--self-test" ]; then
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  printf '%s\n' '#!/usr/bin/env bash' 'docker run --rm example/tool:1.0 true' > "${tmp}/planted-clean.sh"
  printf '%s\n' '#!/usr/bin/env bash' 'docker run --rm "" true || true' > "${tmp}/planted-empty-image.sh"
  printf '%s\n' '#!/usr/bin/env bash' 'exit 3' > "${tmp}/planted-exit.sh"
  # shellcheck disable=SC2016  # the dollar signs are the test input
  printf '%s\n' '#!/usr/bin/env bash' 'IMG=""' 'docker run --rm $IMG tool --flag || true' > "${tmp}/planted-unquoted-empty.sh"
  printf '%s\n' '#!/usr/bin/env bash' 'docker run --rm --imaginary=1 example/tool:1.0 true || true' > "${tmp}/planted-unknown-eq-option.sh"
  rc=0
  out=$(FAKE_DOCKER_CASE_DIR="$tmp" GITHUB_ACTIONS='' "$0" 2>&1) || rc=$?
  ok=true
  [ "$rc" -eq 1 ] || { echo "self-test: runner exited ${rc}, expected 1"; ok=false; }
  for want in 'PASS planted-clean' 'FAIL planted-empty-image' 'bad docker call: EMPTY_IMAGE run' \
              'FAIL planted-exit' 'FAIL planted-unquoted-empty' 'bad docker call: NO_IMAGE run untagged=tool' \
              'FAIL planted-unknown-eq-option' 'bad docker call: UNKNOWN_OPTION run --imaginary=1' \
              '1 passed, 4 failed'; do
    grep -qF "$want" <<<"$out" || { echo "self-test: missing '${want}'"; ok=false; }
  done
  if $ok; then
    echo "self-test: a swallowed empty image (quoted or not), an unknown --name=value option and a failing case are reported, a clean case passes: PASS"
    exit 0
  fi
  printf '%s\n' "$out"
  echo "self-test: FAIL"
  exit 1
fi

cases=()
if [ $# -gt 0 ]; then
  for n in "$@"; do cases+=("${CASE_DIR}/${n%.sh}.sh"); done
else
  for c in "${CASE_DIR}"/*.sh; do
    [ -f "$c" ] && cases+=("$c")
  done
fi
if [ "${#cases[@]}" -eq 0 ]; then
  echo "ERROR: no case files in ${CASE_DIR}; a suite with nothing to run proves nothing." >&2
  exit 1
fi

group() { [ -n "${GITHUB_ACTIONS:-}" ] && echo "::group::$*" || echo "=== $*"; }
endgroup() { [ -n "${GITHUB_ACTIONS:-}" ] && echo "::endgroup::" || true; }

timeout_cmd=()
if command -v timeout >/dev/null 2>&1; then
  timeout_cmd=("$(command -v timeout)" "$CASE_TIMEOUT")
fi

passed=() failed=()
for c in "${cases[@]}"; do
  name=$(basename "$c" .sh)
  if [ ! -f "$c" ]; then
    echo "FAIL ${name}: no such case file" >&2
    failed+=("$name")
    continue
  fi
  work=$(mktemp -d "${TMPDIR:-/tmp}/fake-docker.${name}.XXXXXX")
  mkdir -p "${work}/home"
  : > "${work}/docker.log"

  rc=0
  env -i \
    HOME="${work}/home" \
    PATH="${FAKE_BIN}:/usr/local/bin:/usr/bin:/bin" \
    LANG=C.UTF-8 TERM=dumb TMPDIR="${work}" \
    REPO_ROOT="$ROOT" CASE_WORK="$work" FAKE_DOCKER_LOG="${work}/docker.log" \
    "${timeout_cmd[@]}" bash "$c" > "${work}/case.out" 2>&1 || rc=$?

  bad_calls=$(grep -E '^(EMPTY_IMAGE|NO_IMAGE|UNKNOWN_OPTION)' "${work}/docker.log" || true)

  group "${name} (exit ${rc})"
  cat "${work}/case.out"
  echo "--- docker/wget/curl calls ($(wc -l < "${work}/docker.log" | tr -d ' ')):"
  cat "${work}/docker.log"
  endgroup

  if [ "$rc" -eq 0 ] && [ -z "$bad_calls" ]; then
    echo "PASS ${name}"
    passed+=("$name")
  else
    echo "FAIL ${name}"
    if [ "$rc" -eq 124 ]; then echo "  timed out after ${CASE_TIMEOUT}s"; fi
    if [ "$rc" -ne 0 ]; then
      echo "  exit ${rc}"
      grep -E 'ASSERT FAIL|unbound variable' "${work}/case.out" | sed 's/^/  /' | head -20 || true
    fi
    if [ -n "$bad_calls" ]; then printf '%s\n' "$bad_calls" | sed 's/^/  bad docker call: /'; fi
    failed+=("$name")
  fi
  rm -rf "$work"
done

echo ""
echo "Fake-docker suite: ${#passed[@]} passed, ${#failed[@]} failed."
if [ "${#failed[@]}" -gt 0 ]; then
  printf '  failed: %s\n' "${failed[@]}"
  exit 1
fi
