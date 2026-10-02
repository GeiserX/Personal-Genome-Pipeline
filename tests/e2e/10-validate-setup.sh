#!/usr/bin/env bash
# validate-setup.sh runs to its summary. It still reports [FAIL] items here
# (the VEP cache and most images are absent on the runner), so the check is
# that it gets through every section, not its exit code.
. "$(dirname "$0")/lib.sh"

run_step validate-setup.sh "$SAMPLE"
OUT=$(cat "$STEP_LOG")
check "no 'unbound variable' error" lacks 'unbound variable' "$OUT"
check "gets past the Docker image list" has '=== Sample Data' "$OUT"
check "reaches its summary" has '=== Summary ===' "$OUT"

finish
