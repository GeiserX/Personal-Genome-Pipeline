#!/usr/bin/env bash
# benchmark-variants.sh --giab: the step-03 calls against both GIAB HG002
# truth sets, v4.2.1 and v5.0q, each downloaded and checked by md5. --regions
# with the chr20 slice keeps each set's benchmark regions inside the slice.
. "$(dirname "$0")/lib.sh"

# NCBI's GIAB folder has answered 404 for several minutes at a time and then
# come back (2026-10-07); fetch waits out up to 14 minutes per file.
export FETCH_TRIES=15 FETCH_WAIT=60

SLICE="${GENOME_DIR}/giab-slice-chr20.bed"
printf 'chr20\t10000000\t10500000\n' > "$SLICE"
for set in v4.2.1 v5.0q; do
  run_step benchmark-variants.sh "$SAMPLE" --giab "$set" --regions "$SLICE"
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
