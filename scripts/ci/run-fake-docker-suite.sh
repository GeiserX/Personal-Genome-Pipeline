#!/usr/bin/env bash
# run-fake-docker-suite.sh: run the pipeline's shell entry points against a
# fake docker, wget, curl and df (scripts/ci/fake-docker/), one case file at a
# time from tests/fake-docker/*.sh.
#
# Usage: scripts/ci/run-fake-docker-suite.sh [case-name ...]
#
# Each case runs under `env -i` in its own temp directory (see
# scripts/ci/fake-docker/lib.sh for what it receives). A case fails when it
# exits non-zero, runs longer than CASE_TIMEOUT seconds (default 300), or when
# the fake docker logged an empty image name, a missing image or an unknown
# docker option. Adding a case is adding a file; this runner never changes.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
FAKE_BIN="${ROOT}/scripts/ci/fake-docker"
CASE_TIMEOUT=${CASE_TIMEOUT:-300}

cases=()
if [ $# -gt 0 ]; then
  for n in "$@"; do cases+=("${ROOT}/tests/fake-docker/${n%.sh}.sh"); done
else
  for c in "${ROOT}"/tests/fake-docker/*.sh; do
    [ -f "$c" ] && cases+=("$c")
  done
fi
if [ "${#cases[@]}" -eq 0 ]; then
  echo "ERROR: no case files in tests/fake-docker/; a suite with nothing to run proves nothing." >&2
  exit 1
fi

group() { [ -n "${GITHUB_ACTIONS:-}" ] && echo "::group::$*" || echo "=== $*"; }
endgroup() { [ -n "${GITHUB_ACTIONS:-}" ] && echo "::endgroup::" || true; }

timeout_cmd=()
command -v timeout >/dev/null 2>&1 && timeout_cmd=(timeout "$CASE_TIMEOUT")

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
    [ "$rc" -eq 124 ] && echo "  timed out after ${CASE_TIMEOUT}s"
    [ "$rc" -ne 0 ] && grep -E 'ASSERT FAIL|unbound variable' "${work}/case.out" | sed 's/^/  /' | head -20
    [ -n "$bad_calls" ] && printf '%s\n' "$bad_calls" | sed 's/^/  bad docker call: /'
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
