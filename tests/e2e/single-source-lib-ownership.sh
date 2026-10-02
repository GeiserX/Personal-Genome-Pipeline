#!/usr/bin/env bash
# Steps run their containers as the calling user (run_in, scripts/lib/common.sh),
# so after the base cases nothing under a sample directory belongs to anyone
# else, root included.
#
# Then a sample whose vcf/ an older version left owned by root: step 11 cannot
# write there, fails, and says how to take the sample directory back.
. "$(dirname "$0")/lib.sh"

ME=$(id -u)
for d in "$SAMPLE" "${SAMPLE}cyrius"; do
  find "${GENOME_DIR}/${d}" > "${CASE_TMP}/all-${d}.txt" 2>/dev/null || true
  check_ge "paths under ${d}/ (written by the base cases)" "$(grep -c . "${CASE_TMP}/all-${d}.txt" || true)" 3
  find "${GENOME_DIR}/${d}" ! -user "$ME" > "${CASE_TMP}/foreign-${d}.txt" 2>/dev/null || true
  check_eq "paths under ${d}/ not owned by uid ${ME}" "$(grep -c . "${CASE_TMP}/foreign-${d}.txt" || true)" 0
  head -n 20 "${CASE_TMP}/foreign-${d}.txt"
done

# A root-owned vcf/ the way a run before scripts/lib/common.sh left it.
L="${SAMPLE}legacy"
in_genome "$BCFTOOLS_IMAGE" bash -c "mkdir -p '${L}/vcf' &&
  cp '${SAMPLE}/vcf/${SAMPLE}.vcf.gz' '${L}/vcf/${L}.vcf.gz' &&
  cp '${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi' '${L}/vcf/${L}.vcf.gz.tbi'"
check "${L}/vcf/ is owned by root" test "$(find "${GENOME_DIR}/${L}/vcf" -maxdepth 0 -user 0 | grep -c .)" -eq 1
run_step 11-roh-analysis.sh "$L"
check "step 11 fails when vcf/ belongs to root (exit ${STEP_RC})" test "$STEP_RC" -ne 0
check "step 11 prints the chown that gives the sample directory back" \
  has "sudo chown -R \"${ME}:$(id -g)\" \"${GENOME_DIR}/${L}\"" "$(cat "$STEP_LOG")"
in_genome "$BCFTOOLS_IMAGE" rm -rf "$L"

finish
