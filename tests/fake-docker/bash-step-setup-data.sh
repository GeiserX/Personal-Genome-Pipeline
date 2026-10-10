#!/usr/bin/env bash
# setup.sh creates the reference's sequence dictionary (Picard and GATK need
# it; chip-to-vcf.sh now checks for it) and tries to install the small pinned
# data files. The fake downloads fail their checksums, which must not stop
# setup: each file gets a [WARN], and validate-setup.sh names what is missing
# and what the step does without it.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook
# CreateSequenceDictionary writes a SAM header; the generic hook would leave an
# empty file, which setup.sh rightly refuses.
cat > "${CASE_WORK}/dict-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
case " ${*:2} " in
  *" CreateSequenceDictionary "*)
    out=$(sed -n 's/.* -O \([^ ]*\) .*/\1/p' <<<" ${*:2} ")
    printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "$(host_path "$out")"
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/dict-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/dict-hook"

# chip-to-vcf.sh stops on a missing .dict before it converts anything.
mkdir -p "${GENOME_DIR}/reference_hg19" "${GENOME_DIR}/liftover" "${GENOME_DIR}/chip1/raw"
for f in reference_hg19/human_g1k_v37.fasta reference_hg19/human_g1k_v37.fasta.fai liftover/hg19ToHg38.over.chain.gz; do
  printf 'placeholder\n' > "${GENOME_DIR}/${f}"
done
printf 'rs1\t1\t100\tAG\n' > "${GENOME_DIR}/chip1/raw/chip1_raw.txt"
: > "$FAKE_DOCKER_LOG"
run_expect 1 chip-no-dict "${SCRIPTS}/chip-to-vcf.sh" chip1
output_has chip-no-dict 'Required file not found: .*GRCh38_no_alt_analysis_set\.dict'
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "chip-to-vcf.sh started a container without the sequence dictionary"; fi

run_expect 0 setup "${SCRIPTS}/setup.sh" "$GENOME_DIR"
docker_log_has '^run image=[^ ]*gatk.* gatk CreateSequenceDictionary -R /genome/reference/GRCh38_no_alt_analysis_set\.fasta ' \
  "setup.sh did not create the sequence dictionary"
[ -f "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.dict" ] || fail "setup.sh left no GRCh38_no_alt_analysis_set.dict"
for what in 'Delly exclude map' 'GRCh38 chromosome bands' 'IPD-IMGT/HLA' 'GENCODE'; do
  output_has setup "\[WARN\] Could not install the ${what}"
done
docker_log_has '^curl https://raw\.githubusercontent\.com/dellytools/delly/[0-9a-f]{40}/excludeTemplates/human\.hg38\.excl\.tsv ' \
  "setup.sh did not fetch Delly's exclude map from a pinned commit"

run_rc validate "${SCRIPTS}/validate-setup.sh" sample1
output_has validate 'Delly exclude map \(step 19 runs without it'
output_has validate 'IPD-IMGT/HLA [0-9.]+ \(step 08 is skipped without it\) not found'
output_has validate 'GRCh38 chromosome bands \(step 10 falls back to hg19 bands\) not found'

# --- setup.sh --vep-cache: the VEP cache before the first run ------------------------
# The fake download serves a real tarball of the release versions.env pins,
# with the CHECKSUMS line Ensembl publishes (BSD sum), so install_vep_cache
# runs as it does for real: fetch, check, unpack, move, delete the tarball.
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
TARBALL="homo_sapiens_vep_${VEP_CACHE_RELEASE}_GRCh38.tar.gz"
mkdir -p "${CASE_WORK}/vep-src/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38" "${CASE_WORK}/vep-served"
echo "species homo_sapiens" > "${CASE_WORK}/vep-src/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
tar -czf "${CASE_WORK}/vep-served/${TARBALL}" -C "${CASE_WORK}/vep-src" homo_sapiens
printf '%s 1 %s\n' "$(sum "${CASE_WORK}/vep-served/${TARBALL}" | awk '{print $1}')" "$TARBALL" > "${CASE_WORK}/vep-served/CHECKSUMS"
[ ! -e "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38" ] || fail "the VEP cache exists before --vep-cache ran"

: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/vep-served" run_expect 0 vep-cache "${SCRIPTS}/setup.sh" --vep-cache "$GENOME_DIR"
docker_log_has "^curl [^ ]*/release-${VEP_CACHE_RELEASE}/variation/indexed_vep_cache/${TARBALL} -> " \
  "setup.sh --vep-cache did not download the release-${VEP_CACHE_RELEASE} cache"
[ -f "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ] \
  || fail "setup.sh --vep-cache did not install the release-${VEP_CACHE_RELEASE} cache"
[ ! -e "${GENOME_DIR}/vep_cache/${TARBALL}" ] || fail "setup.sh --vep-cache kept the cache tarball"
output_has vep-cache "\[OK\] VEP ${VEP_CACHE_RELEASE} cache: "
output_lacks vep-cache '=== Phase'
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "setup.sh --vep-cache started a container"; fi
run_expect 0 vep-cache-again "${SCRIPTS}/setup.sh" --vep-cache "$GENOME_DIR"
output_has vep-cache-again "\[OK\] VEP ${VEP_CACHE_RELEASE} cache already present"

# The full setup points at the flag while the cache is missing.
output_has setup "\./scripts/setup\.sh --vep-cache "

# --- opt-in data installers pull their step's image ------------------------------------
# The data is already in place (the fake downloads fail their checksums), so
# --parascopy-data and --yleaf-data go straight to the image, which is missing.
P="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
mkdir -p "${P}/homology_table" "${P}/models_GRCh38_1KGP/EUR"
printf 'x\n' > "${P}/homology_table/GRCh38.bed.gz"
printf 'x\n' > "${P}/models_GRCh38_1KGP/EUR/SMN1.gz"
Y="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}/data"
mkdir -p "${Y}/hg38" "${Y}/hg_prediction_tables"
printf 'x\n' > "${Y}/hg38/new_positions.txt"
printf '{}\n' > "${Y}/hg_prediction_tables/tree.json"

: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_MISSING_IMAGES='parascopy|yleaf' run_rc validate-optin "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-optin "\[WARN\].*Opt-in step 35: its data is installed but its image is not pulled"
output_has validate-optin "\[WARN\].*Opt-in step 37: its data is installed but its image is not pulled"
output_has validate-optin "setup\.sh --parascopy-data "

: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_MISSING_IMAGES='parascopy|yleaf' run_expect 0 parascopy "${SCRIPTS}/setup.sh" --parascopy-data "$GENOME_DIR"
docker_log_has "^pull :: pull ${PARASCOPY_IMAGE//./\\.} " "setup.sh --parascopy-data did not pull ${PARASCOPY_IMAGE}"
output_has parascopy "\[OK\] ${PARASCOPY_IMAGE//./\\.} \(step 35\) pulled"
: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_MISSING_IMAGES='parascopy|yleaf' run_expect 0 yleaf "${SCRIPTS}/setup.sh" --yleaf-data "$GENOME_DIR"
docker_log_has "^pull :: pull ${YLEAF_IMAGE//./\\.} " "setup.sh --yleaf-data did not pull ${YLEAF_IMAGE}"

# With both images present, validate-setup reports them as ready.
run_rc validate-optin-ok "${SCRIPTS}/validate-setup.sh" sample1
output_has validate-optin-ok "\[OK\].*Opt-in step 35: its data and image "
output_has validate-optin-ok "\[OK\].*Opt-in step 37: its data and image "
