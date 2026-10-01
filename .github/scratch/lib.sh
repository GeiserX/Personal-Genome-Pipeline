# shellcheck shell=bash
# Shared helpers for the scratch checks (sourced). Temporary: deleted before merge.
set -uo pipefail

REPO="${GITHUB_WORKSPACE:-$(pwd)}"
NEW="$REPO"
OLD=/tmp/old-main          # origin/main scripts + versions.env, for red-first runs
SCRATCH="$REPO/.github/scratch"
RESULTS=/tmp/scratch-results.tsv
LOGS=/tmp/scratch-logs
mkdir -p "$LOGS"
: > "$RESULTS"
BAD=0

prelude() {
  # Old tree from origin/main for the red-first controls
  git -C "$REPO" fetch -q origin main
  rm -rf "$OLD" && mkdir -p "$OLD"
  git -C "$REPO" archive origin/main scripts versions.env | tar -x -C "$OLD"
  echo "old tree: $(git -C "$REPO" rev-parse origin/main)"

  # docker shim: the scripts ask for up to 8 CPUs, a hosted runner has 4.
  # It clamps --cpus to nproc and can add DOCKER_RUN_EXTRA (e.g. --network none).
  mkdir -p /tmp/shim
  cat > /tmp/shim/docker <<'EOF'
#!/usr/bin/env bash
max=$(nproc)
args=()
while [ $# -gt 0 ]; do
  if [ "$1" = "--cpus" ]; then
    v="$2"; shift 2
    if awk "BEGIN{exit !($v > $max)}"; then v=$max; fi
    args+=(--cpus "$v"); continue
  fi
  args+=("$1"); shift
done
if [ "${args[0]:-}" = run ] && [ -n "${DOCKER_RUN_EXTRA:-}" ]; then
  # shellcheck disable=SC2206
  args=(run $DOCKER_RUN_EXTRA "${args[@]:1}")
fi
exec /usr/bin/docker "${args[@]}"
EOF
  chmod +x /tmp/shim/docker
  export PATH="/tmp/shim:$PATH"
  sudo apt-get -qq update >/dev/null && sudo apt-get -qq install -y tabix samtools bcftools >/dev/null
  echo "host tools: $(samtools --version | head -1), $(bcftools --version | head -1)"
}

# expect_ok <label> cmd...     the command must exit 0
expect_ok() {
  local label="$1"; shift
  local log="$LOGS/${label// /_}.log" rc=0
  "$@" > "$log" 2>&1 || rc=$?
  if [ "$rc" -eq 0 ]; then
    printf '%s\t%s\t%s\n' "$label" "exit 0 (expected 0)" "PASS" >> "$RESULTS"
  else
    printf '%s\t%s\t%s\n' "$label" "exit $rc (expected 0)" "UNEXPECTED" >> "$RESULTS"; BAD=1
  fi
  echo "::group::$label (exit $rc)"; tail -n 60 "$log"; echo "::endgroup::"
  return 0
}

# expect_fail <label> cmd...   the command must exit non-zero (red-first control)
expect_fail() {
  local label="$1"; shift
  local log="$LOGS/${label// /_}.log" rc=0
  "$@" > "$log" 2>&1 || rc=$?
  if [ "$rc" -ne 0 ]; then
    printf '%s\t%s\t%s\n' "$label" "exit $rc (expected non-zero)" "PASS" >> "$RESULTS"
  else
    printf '%s\t%s\t%s\n' "$label" "exit 0 (expected non-zero)" "UNEXPECTED" >> "$RESULTS"; BAD=1
  fi
  echo "::group::$label (exit $rc)"; tail -n 60 "$log"; echo "::endgroup::"
  return 0
}

# observe <label> cmd...      records the exit code without judging it
observe() {
  local label="$1"; shift
  local log="$LOGS/${label// /_}.log" rc=0
  "$@" > "$log" 2>&1 || rc=$?
  printf '%s\t%s\t%s\n' "$label" "exit $rc" "INFO" >> "$RESULTS"
  echo "::group::$label (exit $rc)"; tail -n 60 "$log"; echo "::endgroup::"
  return 0
}

# check <label> <observed> <expected-description> <0|1 pass>
check() {
  if [ "$4" = 1 ]; then
    printf '%s\t%s\t%s\n' "$1" "$2 ($3)" "PASS" >> "$RESULTS"
  else
    printf '%s\t%s\t%s\n' "$1" "$2 ($3)" "UNEXPECTED" >> "$RESULTS"; BAD=1
  fi
  echo "CHECK $1: $2 [expected $3]"
}

finish() {
  {
    echo "| check | observed | verdict |"
    echo "|---|---|---|"
    awk -F'\t' '{print "| " $1 " | " $2 " | " $3 " |"}' "$RESULTS"
  } >> "${GITHUB_STEP_SUMMARY:-/dev/stdout}"
  cat "$RESULTS"
  exit "$BAD"
}
