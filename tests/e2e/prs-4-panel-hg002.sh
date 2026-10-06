#!/usr/bin/env bash
# HG002 on the 1000 Genomes panel case prs-3 installed: the checks of case
# prs-2 (ancestry table, percentile and its group, both reports) with the
# panel users install instead of the synthetic one. HG002 is of European
# ancestry, but the fixture holds a few megabases, so the label is printed,
# not checked. Runs when prs-3 did (monthly and dispatched runs, or
# PGSC_MEASURE=1); case prs-2 removes the panel afterwards.
. "$(dirname "$0")/lib.sh"
. "${REPO}/versions.env"

PANEL="${GENOME_DIR}/reference/pgsc_calc/${PGSC_PANEL}.tar.zst"
if [ ! -s "$PANEL" ]; then
  if [ "${GITHUB_EVENT_NAME:-}" = pull_request ] && [ "${PGSC_MEASURE:-}" != 1 ]; then
    echo "pull request run: case prs-3 did not install the 1000 Genomes panel, so there is nothing to project onto."
    finish
  fi
  fail "case prs-3 left no panel at ${PANEL}"
  finish
fi
t0=$(date +%s)
PANEL_NAME="$PGSC_PANEL" bash "$(dirname "$0")/prs-2-ancestry.sh"
check_eq "case prs-2 on ${PGSC_PANEL} passes" "$?" 0
printf 'fixture_run_seconds_with_%s\t%s\n' "$PGSC_PANEL" "$(( $(date +%s) - t0 ))" >> "${E2E_WORK}/logs/pgsc-measure.tsv"
finish
