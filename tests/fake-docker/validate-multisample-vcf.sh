#!/usr/bin/env bash
# validate-setup.sh and a sample VCF with more than one sample column:
#   - two samples (a joint-called file) fail, naming the count and the
#     samples, with the bcftools command that keeps one;
#   - more than five samples list the first five, then "...";
#   - one sample passes and is named.
# The VCFs are the synthetic ones in tests/fixtures/vcf/. The fake docker
# answers `bcftools view -h` with the header of the sample's VCF through a
# run hook.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
VCF="${GENOME_DIR}/sample1/vcf/sample1.vcf.gz"
FIX="${REPO_ROOT}/tests/fixtures/vcf"

cat > "${CASE_WORK}/vcf-hook" <<'HOOK'
#!/usr/bin/env bash
shift   # the image
if [ "${1:-}" = bcftools ] && [ "${2:-}" = view ] && [ "${3:-}" = -h ]; then
  gzip -dc "${GENOME_DIR}${4#/genome}" | grep '^#'
fi
exit 0
HOOK
chmod +x "${CASE_WORK}/vcf-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/vcf-hook"

# --- two samples ------------------------------------------------------------------
gzip -c "${FIX}/two_samples.vcf" > "$VCF"
run_expect 1 two "${SCRIPTS}/validate-setup.sh" sample1
output_has two '\[FAIL\].*VCF holds 2 samples \(SAMPLE_A, SAMPLE_B\); the pipeline analyses one sample per run'
output_has two 'bcftools view -s <name> -a -c 1 -Oz -o /genome/<name>/vcf/<name>\.vcf\.gz /genome/sample1/vcf/sample1\.vcf\.gz'
output_has two 'bcftools index -t /genome/<name>/vcf/<name>\.vcf\.gz'
output_has two '\[OK\].*VCF genome build: GRCh38'
[ "$(grep -c '^ *\[FAIL\]' "${CASE_WORK}/two.out")" -eq 1 ] || fail "two: expected exactly one [FAIL], the sample count"

# --- seven samples: the first five, then ... ----------------------------------------
awk -F'\t' -v OFS='\t' '/^#CHROM/ { $10 = "S1"; $11 = "S2\tS3\tS4\tS5\tS6\tS7" } { print }' \
  "${FIX}/two_samples.vcf" | grep '^#' | gzip -c > "$VCF"
run_expect 1 seven "${SCRIPTS}/validate-setup.sh" sample1
output_has seven '\[FAIL\].*VCF holds 7 samples \(S1, S2, S3, S4, S5, \.\.\.\)'

# --- one sample ----------------------------------------------------------------------
gzip -c "${FIX}/one_sample.vcf" > "$VCF"
run_expect 0 one "${SCRIPTS}/validate-setup.sh" sample1
output_has one '\[OK\].*VCF holds one sample \(SAMPLE_A\)'
output_lacks one 'samples \('

echo "validate-multisample-vcf: a VCF with two or seven samples fails with the count and the names; one sample passes"
