#!/usr/bin/env bash
# Step 37 (Yleaf) reads the sex step 16 infers. On the fixture's slices
# indexcov reads HG002 as female (case 34), so the step skips HG002 with one
# line. Two copies of the sample (the BAM hard-linked) carry an indexcov .ped
# of a male and of a female sample: the male one runs Yleaf on the chrY slice
# (chrY:2.7-3.0 Mb) and must write a haplogroup or say "insufficient
# markers"; the female one is skipped.
. "$(dirname "$0")/lib.sh"

run_step 37-y-haplogroup.sh "$SAMPLE"
check_step_exit 37-y-haplogroup.sh
check "HG002 is skipped on the sex indexcov infers from the slices" \
  has '^Step 37 skipped: indexcov infers (female|unknown) from the reads' "$(cat "$STEP_LOG")"

PED_HEAD=$'#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY'
copy() {  # copy NAME SEX_CODE: a copy of the sample with an indexcov .ped of that sex
  local t=$1
  mkdir -p "${GENOME_DIR}/${t}/aligned" "${GENOME_DIR}/${t}/indexcov"
  ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" "${GENOME_DIR}/${t}/aligned/${t}_sorted.bam"
  ln -f "${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai" "${GENOME_DIR}/${t}/aligned/${t}_sorted.bam.bai"
  printf '%s\n%s\t%s\t-9\t-9\t%s\t-9\t1.0\t1.0\n' "$PED_HEAD" "$t" "$t" "$2" > "${GENOME_DIR}/${t}/indexcov/indexcov-indexcov.ped"
}

M="${SAMPLE}ym"
copy "$M" 1
run_step 37-y-haplogroup.sh "$M"
check_step_exit 37-y-haplogroup.sh
OUT="${GENOME_DIR}/${M}/y_haplogroup/${M}_y_haplogroup.txt"
cat "$OUT" 2>/dev/null
check_eq "Yleaf's table has a header and one sample row" "$(grep -c . "$OUT" 2>/dev/null || true)" 2
check "the step writes a haplogroup or a clean 'insufficient markers' line" \
  has '^Y haplogroup: ([A-T][^ ]* \([0-9]+ markers|insufficient markers \(Yleaf found [0-9]+ Y markers)' "$(cat "$STEP_LOG")"

F="${SAMPLE}yf"
copy "$F" 2
run_step 37-y-haplogroup.sh "$F"
check_step_exit 37-y-haplogroup.sh
check "a female sample is skipped with one line" has '^Step 37 skipped: indexcov infers female from the reads' "$(cat "$STEP_LOG")"
check "and gets no output" test ! -e "${GENOME_DIR}/${F}/y_haplogroup"

finish
