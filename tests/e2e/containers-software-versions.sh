#!/usr/bin/env bash
# The Nextflow case's pipeline_info/software_versions.yml names every process
# that ran, each with the tag of the image conf/containers.config gives it
# (generated from versions.env). Reads the run of 60-nextflow.sh.
. "$(dirname "$0")/lib.sh"

INFO="${GENOME_DIR}/nf-results/pipeline_info"
TRACE=$(find "$INFO" -name 'trace_*.txt' 2>/dev/null | LC_ALL=C sort | awk 'END {print}')
# Process names of the tasks that completed: "PGX:PHARMCAT (HG002)" -> PHARMCAT.
RAN=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
                  $c["status"] == "COMPLETED" || $c["status"] == "CACHED" {n = $c["name"]; sub(/ \(.*$/, "", n); sub(/.*:/, "", n); print n}' \
  "${TRACE:-/dev/null}" 2>/dev/null | LC_ALL=C sort -u)
N_RAN=$(grep -c . <<< "$RAN" || true)
echo "processes that completed: $(tr '\n' ' ' <<< "$RAN")"
check_ge "processes the Nextflow case completed" "$N_RAN" 8
check "software_versions.yml exists" test -s "${INFO}/software_versions.yml"

# shellcheck disable=SC2086  # one argument per process name
OUT=$("${REPO}/scripts/ci/gen-containers-config.sh" --check-versions "${INFO}/software_versions.yml" $RAN 2>&1)
echo "$OUT"
check_eq "processes whose version line is their image tag" "$(grep -c '^OK ' <<< "$OUT" || true)" "$N_RAN"
check_eq "processes with a wrong, missing or unknown version" "$(grep -c '^FAIL' <<< "$OUT" || true)" 0

finish
