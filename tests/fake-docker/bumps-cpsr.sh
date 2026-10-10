#!/usr/bin/env bash
# Step 17 (CPSR) on PCGR 2.3: the argv the image gets, and the refusal a user
# with only the PCGR 2.2.5 data sees.
#   - the image is PCGR_IMAGE from ghcr.io (Docker Hub's sigven/pcgr stops at
#     2.2.5), the bundle mounted at /mnt/bundle is pcgr_data/PCGR_DATA_BUNDLE,
#     and VEP gets the cache directory that holds PCGR_VEP_CACHE_RELEASE;
#   - --classify_all is not passed: CPSR 2.3 dropped it and stops on it;
#   - with only the 2.2.5 data (VEP 113 cache, bundle 20250314) the step stops
#     before docker and names the cache and the bundle it needs;
#   - a sample called S1 (CPSR 2.3.2 takes 3 to 40 characters) reaches CPSR
#     as S1_cpsr, and its files end up named S1.cpsr.*.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

# The 2.2.5 data only: the step must refuse, before any docker run.
mkdir -p "${GENOME_DIR}/vep_cache/homo_sapiens/113_GRCh38" "${GENOME_DIR}/pcgr_data/20250314/data"
echo "species homo_sapiens" > "${GENOME_DIR}/vep_cache/homo_sapiens/113_GRCh38/info.txt"
: > "$FAKE_DOCKER_LOG"
run_rc old-cache "${SCRIPTS}/17-cpsr.sh" sample1
[ "$RC" -ne 0 ] || fail "step 17 ran with only the VEP 113 cache"
output_has old-cache "VEP release-${PCGR_VEP_CACHE_RELEASE} cache not found"
output_has old-cache "homo_sapiens_vep_${PCGR_VEP_CACHE_RELEASE}_GRCh38\.tar\.gz"
mkdir -p "${GENOME_DIR}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38"
echo "species homo_sapiens" > "${GENOME_DIR}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt"
run_rc old-bundle "${SCRIPTS}/17-cpsr.sh" sample1
[ "$RC" -ne 0 ] || fail "step 17 ran with only the 20250314 bundle"
output_has old-bundle "pcgr_ref_data\.${PCGR_DATA_BUNDLE}\.grch38\.tgz"
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "step 17 started a container without the data it needs"; fi

# The data of versions.env: the run.
mkdir -p "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data"
: > "$FAKE_DOCKER_LOG"
run_expect 0 cpsr "${SCRIPTS}/17-cpsr.sh" sample1
echo "docker run of step 17:"
grep '^run ' "$FAKE_DOCKER_LOG"
case "$PCGR_IMAGE" in
  ghcr.io/sigven/pcgr:*) ;;
  *) fail "PCGR_IMAGE is ${PCGR_IMAGE}, not an image of ghcr.io/sigven/pcgr" ;;
esac
docker_log_has "^run image=${PCGR_IMAGE//./\\.} " "step 17 did not run PCGR_IMAGE (${PCGR_IMAGE})"
docker_log_has "^run image=[^ ]* :: .* -v [^ ]*/pcgr_data/${PCGR_DATA_BUNDLE}:/mnt/bundle " \
  "step 17 did not mount pcgr_data/${PCGR_DATA_BUNDLE} at /mnt/bundle"
docker_log_has "^run image=[^ ]* :: .* -v [^ ]*/vep_cache:/mnt/\.vep .* --vep_dir /mnt/\.vep " \
  "step 17 did not give CPSR the VEP cache directory"
docker_log_has "^run image=[^ ]* :: .* cpsr .*--panel_id 0 .*--secondary_findings .*--force_overwrite" \
  "step 17 did not pass --panel_id 0, --secondary_findings and --force_overwrite"
if grep -q -- '--classify_all' "$FAKE_DOCKER_LOG"; then fail "step 17 passed --classify_all, which CPSR 2.3 refuses"; fi
# A flag line of the module's script block (the header comment names it too).
if grep -qE -- '^[[:space:]]+--classify_all' "${REPO_ROOT}/modules/local/cpsr/main.nf"; then
  fail "the Nextflow CPSR module passes --classify_all, which CPSR 2.3 refuses"
fi
# A two-character sample: CPSR gets S1_cpsr, and the files it writes under
# that name are renamed to S1's. The hook writes them as CPSR would.
seed_sample "$GENOME_DIR" S1
cat > "${CASE_WORK}/cpsr-hook" <<'HOOK'
#!/usr/bin/env bash
# Called as <image> <command...>, with the -v specs in FAKE_DOCKER_VOLUMES
out=$(sed -n 's|^\(.*\):/mnt/outputs$|\1|p' <<<"${FAKE_DOCKER_VOLUMES:-}") sid="" prev=""
for a in "$@"; do [ "$prev" != --sample_id ] || sid=$a; prev=$a; done
[ -z "$out" ] || [ -z "$sid" ] || for e in html classification.tsv.gz; do : > "${out}/${sid}.cpsr.grch38.${e}"; done
HOOK
chmod +x "${CASE_WORK}/cpsr-hook"
: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/cpsr-hook" run_expect 0 short-id "${SCRIPTS}/17-cpsr.sh" S1
docker_log_has "^run image=[^ ]* :: .* cpsr .*--sample_id S1_cpsr " "step 17 did not give CPSR S1_cpsr for sample S1"
for e in html classification.tsv.gz; do
  [ -e "${GENOME_DIR}/S1/cpsr/S1.cpsr.grch38.${e}" ] || fail "step 17 left no S1.cpsr.grch38.${e}: $(ls "${GENOME_DIR}/S1/cpsr")"
  [ ! -e "${GENOME_DIR}/S1/cpsr/S1_cpsr.cpsr.grch38.${e}" ] || fail "step 17 left S1_cpsr.cpsr.grch38.${e}"
done

echo "PASS: step 17 runs ${PCGR_IMAGE} on bundle ${PCGR_DATA_BUNDLE} and the VEP ${PCGR_VEP_CACHE_RELEASE} cache"
