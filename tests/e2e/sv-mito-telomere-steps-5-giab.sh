#!/usr/bin/env bash
# benchmark-variants.sh --giab: the step-03 calls against both GIAB HG002
# truth sets, v4.2.1 and v5.0q. --regions with the chr20 slice keeps each
# set's benchmark regions inside the slice.
#
# The fixture carries both sets' files exactly as GIAB publishes them (since
# fixture-v6), so the case copies them into GENOME_DIR/giab/ and the step
# checks them by md5 there, as it checks a file a user downloaded by hand. No
# download happens: each set's source is pointed at a host that never
# resolves, so a step that tried one would fail the case. NCBI's GIAB tree,
# the only home of v5.0q, has been 404 for a day at a time, and this keeps
# the case independent of it. E2E_GIAB_SETS runs only the sets it names, for
# a check of one leg by hand.
. "$(dirname "$0")/lib.sh"

# .invalid never resolves (RFC 2606): any fetch fails on its first try.
export GIAB_V421_URL=https://giab-v421.blocked.invalid/ GIAB_V5Q_URL=https://ftp-trace.ncbi.nlm.nih.gov.invalid/
export FETCH_TRIES=1 FETCH_WAIT=0

giab_files() {
  case "$1" in
    v4.2.1) echo HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
                 HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed ;;
    v5.0q)  echo HG002_GRCh38_v5.0q_smvar.vcf.gz HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi \
                 HG002_GRCh38_v5.0q_smvar.benchmark.bed ;;
  esac
}

SLICE="${GENOME_DIR}/giab-slice-chr20.bed"
printf 'chr20\t10000000\t10500000\n' > "$SLICE"
mkdir -p "${GENOME_DIR}/giab"
for set in ${E2E_GIAB_SETS:-v4.2.1 v5.0q}; do
  for f in $(giab_files "$set"); do
    check "the fixture holds GIAB ${set}'s ${f}" cp "${FIXTURE_DIR}/${f}" "${GENOME_DIR}/giab/${f}"
  done
  run_step benchmark-variants.sh "$SAMPLE" --giab "$set" --regions "$SLICE"
  check_step_exit "benchmark-variants.sh --giab ${set}"
  check "the step downloaded nothing for ${set}" lacks "Downloading GIAB" "$(cat "$STEP_LOG")"
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
