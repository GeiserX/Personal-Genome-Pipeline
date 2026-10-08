#!/usr/bin/env bash
# benchmark-variants.sh --giab: the step-03 calls against both GIAB HG002
# truth sets, v4.2.1 and v5.0q, each downloaded and checked by md5. --regions
# with the chr20 slice keeps each set's benchmark regions inside the slice.
#
# v4.2.1 comes from GIAB's S3 mirror and must pass. v5.0q has no mirror, only
# NCBI: when every try of one of its downloads gets HTTP 404, the leg prints
# one SKIPPED line (also written to the run's summary) and does not fail the
# case. Any other failure of either leg (a wrong md5, a network error, hap.py,
# the recall or the precision) fails it. E2E_GIAB_SETS runs only the sets it
# names, for a check of one leg by hand.
. "$(dirname "$0")/lib.sh"

# NCBI's GIAB folder has answered 404 for several minutes at a time and then
# come back (2026-10-07); fetch waits out up to 14 minutes per file.
export FETCH_TRIES=${FETCH_TRIES:-15} FETCH_WAIT=${FETCH_WAIT:-60}

# upstream_404: true when the step stopped only because one download got HTTP
# 404 on every try: one ERROR line, the fetch's "could not download", and no
# curl error but 404. Prints the number of tries.
upstream_404() {
  local last done_tries tries
  [ "$STEP_RC" != 0 ] || return 1
  [ "$(grep -c '^ERROR:' "$STEP_LOG")" = 1 ] || return 1
  grep -q '^ERROR: could not download ' "$STEP_LOG" || return 1
  last=$(grep -Eo 'Download attempt [0-9]+/[0-9]+ failed' "$STEP_LOG" | tail -n 1)
  done_tries=$(sed -E 's|.* ([0-9]+)/[0-9]+ .*|\1|' <<<"$last")
  tries=$(sed -E 's|.*/([0-9]+) .*|\1|' <<<"$last")
  [ -n "$tries" ] && [ "$done_tries" = "$tries" ] || return 1
  [ "$(grep -Eo 'curl: \([0-9]+\)[^[:cntrl:]]*' "$STEP_LOG" | grep -vc 'returned error: 404')" = 0 ] || return 1
  [ "$(grep -Eo 'curl: \([0-9]+\)[^[:cntrl:]]*' "$STEP_LOG" | grep -c 'returned error: 404')" -ge "$tries" ] || return 1
  echo "$tries"
}

SLICE="${GENOME_DIR}/giab-slice-chr20.bed"
printf 'chr20\t10000000\t10500000\n' > "$SLICE"
for set in ${E2E_GIAB_SETS:-v4.2.1 v5.0q}; do
  run_step benchmark-variants.sh "$SAMPLE" --giab "$set" --regions "$SLICE"
  if [ "$set" = v5.0q ] && tries=$(upstream_404); then
    line="SKIPPED v5.0q: GIAB truth set unavailable upstream (HTTP 404 after ${tries} tries)"
    echo "::warning::${line}"
    printf '#### Skipped\n\n- %s: %s\n\n' "$CASE_NAME" "$line" >> "$E2E_NOTES"
    continue
  fi
  check_step_exit "benchmark-variants.sh --giab ${set}"
  check "the log names the ${set} truth set" has "Truth set: GIAB HG002 ${set} " "$(cat "$STEP_LOG")"
  check "the ${set} regions are cut to the slice" has "Benchmark regions of ${set} inside .*: [0-9]+ intervals" "$(cat "$STEP_LOG")"
  TSV="${GENOME_DIR}/${SAMPLE}/benchmark/comparison.tsv"
  cat "$TSV" 2>/dev/null
  ROW=$(awk -F'\t' '$1 == "DeepVariant"' "$TSV" 2>/dev/null)
  TP=$(cut -f2 <<<"$ROW") PREC=$(cut -f5 <<<"$ROW") REC=$(cut -f6 <<<"$ROW")
  check_ge "DeepVariant SNP true positives against ${set}" "${TP:-0}" 100
  check "DeepVariant SNP recall against ${set} >= 0.9 (${REC:-none})" awk -v x="${REC:-0}" 'BEGIN {exit !(x + 0 >= 0.9)}'
  check "DeepVariant SNP precision against ${set} >= 0.9 (${PREC:-none})" awk -v x="${PREC:-0}" 'BEGIN {exit !(x + 0 >= 0.9)}'
done

finish
