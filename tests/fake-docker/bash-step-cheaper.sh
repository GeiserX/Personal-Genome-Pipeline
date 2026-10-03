#!/usr/bin/env bash
# Cheaper and quieter steps:
#   - step 14 prepares every chromosome, chrX included, in one container and
#     one bcftools pass each (it started 66 containers in two passes before);
#   - benchmark-variants.sh normalises each caller's VCF once: five callers,
#     five normalisations (it ran two per pair, twenty, before);
#   - MultiQC (script and module) does not ask its server for a new version;
#   - CPSR passes --secondary_findings unless CPSR_SECONDARY_FINDINGS=false.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook

# --- step 14 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 imputation "${SCRIPTS}/14-imputation-prep.sh" sample1
n=$(grep -c '^run ' "$FAKE_DOCKER_LOG" || true)
[ "$n" -eq 1 ] || fail "step 14 started ${n} containers, want 1"
docker_log_has '^run .* chr22 chrX( |$)' "step 14 did not prepare chrX"
docker_log_has "^run .*bcftools view -f PASS,\\\\. -r " "step 14 does not select PASS records per chromosome in one bcftools view"

# --- benchmark: one normalisation per caller --------------------------------------------
S="${GENOME_DIR}/sample1"
mkdir -p "${S}/vcf_gatk" "${S}/vcf_freebayes" "${S}/vcf_strelka2/results/variants" "${S}/vcf_octopus"
for f in vcf_gatk/sample1.vcf.gz vcf_freebayes/sample1.vcf.gz vcf_strelka2/results/variants/variants.vcf.gz vcf_octopus/sample1.vcf.gz; do
  printf 'placeholder\n' > "${S}/${f}"
  printf 'placeholder\n' > "${S}/${f}.tbi"
done
: > "$FAKE_DOCKER_LOG"
run_expect 0 benchmark "${SCRIPTS}/benchmark-variants.sh" sample1
# The normalisation runs in `bash -c`, which the log shows with escaped spaces.
n=$(grep -cE '^run .*bcftools(\\)? norm' "$FAKE_DOCKER_LOG" || true)
[ "$n" -eq 5 ] || fail "benchmark-variants.sh ran ${n} normalisations for 5 callers, want 5"
n=$(grep -c '^run .*bcftools isec ' "$FAKE_DOCKER_LOG" || true)
[ "$n" -eq 10 ] || fail "benchmark-variants.sh ran ${n} comparisons for 5 callers, want 10"
[ -z "$(find "${S}/benchmark" -name '.norm_*')" ] || fail "benchmark-variants.sh left its normalised copies behind"

# --- MultiQC ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 multiqc "${SCRIPTS}/28-multiqc.sh" sample1
docker_log_has '^run image=[^ ]*multiqc.* multiqc .*--no-version-check ' "MultiQC ran without --no-version-check"
grep -q -- '--no-version-check' "${REPO_ROOT}/modules/local/multiqc/main.nf" \
  || fail "the Nextflow MultiQC module runs without --no-version-check"

# --- CPSR secondary findings -------------------------------------------------------------
mkdir -p "${GENOME_DIR}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38" "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data"
echo "species homo_sapiens" > "${GENOME_DIR}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt"
: > "$FAKE_DOCKER_LOG"
run_expect 0 cpsr "${SCRIPTS}/17-cpsr.sh" sample1
docker_log_has '^run .* cpsr .*--secondary_findings ' "CPSR ran without --secondary_findings by default"
: > "$FAKE_DOCKER_LOG"
CPSR_SECONDARY_FINDINGS=false run_expect 0 cpsr-off "${SCRIPTS}/17-cpsr.sh" sample1
if grep -q -- '--secondary_findings' "$FAKE_DOCKER_LOG"; then fail "CPSR_SECONDARY_FINDINGS=false still passed --secondary_findings"; fi
CPSR_SECONDARY_FINDINGS=no run_expect 1 cpsr-bad "${SCRIPTS}/17-cpsr.sh" sample1
output_has cpsr-bad 'CPSR_SECONDARY_FINDINGS must be true or false'
