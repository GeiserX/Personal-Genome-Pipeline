#!/usr/bin/env bash
# Step 13 (VEP) asks for the same annotation as the Nextflow module and leaves
# no stale file for the steps after it:
#   - vep gets --fasta (without it VEP turns HGVS off), --compress_output
#     bgzip, --cache_version of versions.env and --fork THREADS;
#   - the output is <sample>_vep.vcf.gz with its .tbi;
#   - a finished run removes vcfanno's <sample>_annotated.vcf.gz (built from
#     the previous annotation, and preferred by steps 30, 23 and 31) and the
#     uncompressed <sample>_vep.vcf of earlier versions;
#   - a failed VEP run keeps the previous result and the derived files;
#   - an input VCF without its index is refused before VEP starts.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
echo "species homo_sapiens" > "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
V="${GENOME_DIR}/sample1/vep"
mkdir -p "$V"
for f in sample1_annotated.vcf.gz sample1_annotated.vcf.gz.tbi sample1_vep.vcf; do
  printf 'from the previous VEP run\n' > "${V}/${f}"
done
use_output_hook
# FAKE_VEP_FAIL=1 makes vep exit 1; everything else goes to the generic hook.
cat > "${CASE_WORK}/vep-hook" <<'HOOK'
#!/usr/bin/env bash
if [ "${FAKE_VEP_FAIL:-}" = 1 ] && [ "${2:-}" = vep ]; then echo "fake vep: failed" >&2; exit 1; fi
exec "${CASE_WORK:?}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/vep-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/vep-hook"

THREADS=3 run_expect 0 vep "${SCRIPTS}/13-vep-annotation.sh" sample1
docker_log_has '^run image=[^ ]*ensembl-vep.* vep .*--fasta /genome/reference/GRCh38_no_alt_analysis_set\.fasta ' \
  "step 13 did not pass the reference FASTA to vep"
docker_log_has '^run image=[^ ]*ensembl-vep.* --compress_output bgzip ' "step 13 did not ask vep for bgzip output"
docker_log_has "^run image=[^ ]*ensembl-vep.* --cache_version ${VEP_CACHE_RELEASE} " "step 13 did not pass --cache_version ${VEP_CACHE_RELEASE}"
docker_log_has '^run image=[^ ]*ensembl-vep.* --cpus 3 .* --fork 3( |$)' "step 13 did not use THREADS=3 for --cpus and --fork"
[ -s "${V}/sample1_vep.vcf.gz" ] || fail "step 13 wrote no sample1_vep.vcf.gz"
[ -f "${V}/sample1_vep.vcf.gz.tbi" ] || fail "step 13 did not index sample1_vep.vcf.gz"
for f in sample1_annotated.vcf.gz sample1_annotated.vcf.gz.tbi sample1_vep.vcf; do
  [ ! -e "${V}/${f}" ] || fail "step 13 left ${f} from the previous annotation in place"
done

# A failed rerun keeps the previous result and the files built from it.
printf 'built from the result above\n' > "${V}/sample1_annotated.vcf.gz"
BEFORE=$(cksum < "${V}/sample1_vep.vcf.gz")
FAKE_VEP_FAIL=1 run_rc vep-fails "${SCRIPTS}/13-vep-annotation.sh" sample1
[ "$RC" -ne 0 ] || fail "step 13 exited 0 although vep failed"
[ "$(cksum < "${V}/sample1_vep.vcf.gz")" = "$BEFORE" ] || fail "a failed VEP run changed the previous sample1_vep.vcf.gz"
[ -e "${V}/sample1_annotated.vcf.gz" ] || fail "a failed VEP run removed the files built from the previous result"

# No index, no annotation: the VCF may be half written.
rm "${GENOME_DIR}/sample1/vcf/sample1.vcf.gz.tbi"
: > "$FAKE_DOCKER_LOG"
run_rc vep-no-index "${SCRIPTS}/13-vep-annotation.sh" sample1
[ "$RC" -ne 0 ] || fail "step 13 annotated a VCF without its index"
output_has vep-no-index 'sample1\.vcf\.gz\.tbi'
if awk '/^run / && /ensembl-vep/ { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 13 started VEP on a VCF without its index"
fi
