#!/usr/bin/env bash
# Downloads go through fetch (scripts/lib/common.sh): written to <file>.part,
# checked, then renamed.
#   - a dropped connection leaves only the .part file, and the next run
#     downloads the file;
#   - a file whose checksum does not match is removed and the script fails;
#   - the data directory and what setup.sh writes into it are private;
#   - a lock left by a killed run does not block step 26;
#   - step 13 installs the cache of the release versions.env names even when
#     another release is present, and refuses a tarball with a wrong checksum;
#   - step 32 refuses a pypgx bundle that is not the tag versions.env names.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
use_output_hook

mode() { ls -ld "$1" | cut -c1-10; }

# --- 1. a dropped connection -------------------------------------------------
G="${CASE_WORK}/genome"
mkdir -p "$G"
FAKE_DOWNLOAD_PARTIAL='clinvar\.vcf\.gz$' run_rc setup-dropped "${SCRIPTS}/setup.sh" "$G"
[ "$RC" -ne 0 ] || fail "setup.sh exited 0 although the ClinVar download was cut"
[ -f "${G}/clinvar/clinvar.vcf.gz.part" ] || fail "the cut download left no clinvar.vcf.gz.part to resume"
[ ! -e "${G}/clinvar/clinvar.vcf.gz" ] || fail "a cut download was stored as clinvar.vcf.gz"

: > "$FAKE_DOCKER_LOG"
run_expect 0 setup-again "${SCRIPTS}/setup.sh" "$G"
docker_log_has '^curl [^ ]*/clinvar\.vcf\.gz -> ' "the second run did not download ClinVar again"
[ -s "${G}/clinvar/clinvar.vcf.gz" ] || fail "the second run left no clinvar.vcf.gz"
[ ! -e "${G}/clinvar/clinvar.vcf.gz.part" ] || fail "clinvar.vcf.gz.part is still there after a complete download"
if grep -q '^curl [^ ]*Homo_sapiens_assembly38\.fasta -> ' "$FAKE_DOCKER_LOG"; then
  fail "the second run downloaded the reference again although it was complete"
fi

# --- 2. private data directory -------------------------------------------------
[ "$(mode "$G")" = "drwx------" ] || fail "GENOME_DIR is $(mode "$G"), expected drwx------"
for f in reference/Homo_sapiens_assembly38.fasta clinvar/clinvar.vcf.gz; do
  [ "$(mode "${G}/${f}")" = "-rw-------" ] || fail "${f} is $(mode "${G}/${f}"), expected -rw-------"
done

# --- 3. a checksum that does not match ---------------------------------------------
G2="${CASE_WORK}/genome2"
mkdir -p "$G2" "${CASE_WORK}/served"
printf '00000000000000000000000000000000  clinvar.vcf.gz\n' > "${CASE_WORK}/served/clinvar.vcf.gz.md5"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/served" run_rc setup-bad-md5 "${SCRIPTS}/setup.sh" "$G2"
[ "$RC" -ne 0 ] || fail "setup.sh exited 0 although ClinVar did not match its published md5"
output_has setup-bad-md5 'md5 checksum of clinvar\.vcf\.gz is [0-9a-f]+, expected 0{32}'
[ ! -e "${G2}/clinvar/clinvar.vcf.gz" ] || fail "a ClinVar file with a wrong md5 was kept as clinvar.vcf.gz"
[ ! -e "${G2}/clinvar/clinvar.vcf.gz.part" ] || fail "a ClinVar download with a wrong md5 was kept as .part and would be resumed"

G3="${CASE_WORK}/genome3"
mkdir -p "$G3"
REF_FASTA_MD5=00000000000000000000000000000000 run_rc setup-bad-ref "${SCRIPTS}/setup.sh" "$G3"
[ "$RC" -ne 0 ] || fail "setup.sh exited 0 although the reference did not match its recorded md5"
[ ! -e "${G3}/reference/Homo_sapiens_assembly38.fasta" ] || fail "a reference with a wrong md5 was kept"

# --- 4. a lock left by a killed run -------------------------------------------------
export GENOME_DIR="$G"
seed_sample "$G" sample1
( : ) &
dead=$!
wait "$dead"
mkdir -p "${G}/ancestry_ref/.download.lock"
echo "$dead" > "${G}/ancestry_ref/.download.lock/pid"
run_rc ancestry "${SCRIPTS}/26-ancestry.sh" sample1
output_has ancestry 'Removing stale lock'
[ -s "${G}/ancestry_ref/1kg_common_snps.vcf.gz" ] || fail "step 26 did not prepare the panel after taking over the stale lock"
[ -f "${G}/ancestry_ref/1kg_common_snps.vcf.gz.tbi" ] || fail "step 26 left the panel without its index"
[ ! -e "${G}/ancestry_ref/.download.lock" ] || fail "step 26 left its lock behind"
[ ! -e "${G}/ancestry_ref/1kg_common_snps.part.vcf.gz" ] || fail "step 26 left the .part panel behind"

# --- 5. VEP cache of the right release --------------------------------------------------
# Only the CPSR cache is present; step 13 must still install its own release.
OTHER="$PCGR_VEP_CACHE_RELEASE"
[ "$OTHER" != "$VEP_CACHE_RELEASE" ] || fail "versions.env has one VEP release for both steps; this case needs two"
mkdir -p "${G}/vep_cache/homo_sapiens/${OTHER}_GRCh38"
echo "species homo_sapiens" > "${G}/vep_cache/homo_sapiens/${OTHER}_GRCh38/info.txt"
TARBALL="homo_sapiens_vep_${VEP_CACHE_RELEASE}_GRCh38.tar.gz"
mkdir -p "${CASE_WORK}/vep-src/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38" "${CASE_WORK}/vep-served" "${CASE_WORK}/vep-bad"
echo "species homo_sapiens" > "${CASE_WORK}/vep-src/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
tar -czf "${CASE_WORK}/vep-served/${TARBALL}" -C "${CASE_WORK}/vep-src" homo_sapiens
cp "${CASE_WORK}/vep-served/${TARBALL}" "${CASE_WORK}/vep-bad/${TARBALL}"
SUM=$(sum "${CASE_WORK}/vep-served/${TARBALL}" | awk '{print $1}')
printf '11111 1 some_other_file.tar.gz\n%s 1 %s\n' "$SUM" "$TARBALL" > "${CASE_WORK}/vep-served/CHECKSUMS"
printf '%s 1 %s\n' "$(( (10#$SUM + 1) % 65536 ))" "$TARBALL" > "${CASE_WORK}/vep-bad/CHECKSUMS"

: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/vep-bad" run_rc vep-bad-sum "${SCRIPTS}/13-vep-annotation.sh" sample1
[ "$RC" -ne 0 ] || fail "step 13 exited 0 although the cache tarball did not match CHECKSUMS"
[ ! -e "${G}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38" ] || fail "step 13 unpacked a tarball with a wrong checksum"
if awk '/^run / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 13 ran VEP without its cache"
fi

: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/vep-served" run_expect 0 vep "${SCRIPTS}/13-vep-annotation.sh" sample1
docker_log_has "^curl [^ ]*/${TARBALL} -> " "step 13 did not download the release-${VEP_CACHE_RELEASE} cache although only release ${OTHER} was present"
[ -f "${G}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ] || fail "step 13 did not install the release-${VEP_CACHE_RELEASE} cache"
[ ! -e "${G}/vep_cache/${TARBALL}" ] || fail "step 13 kept the cache tarball"
[ -z "$(find "${G}/vep_cache" -maxdepth 1 -name '.extract.*')" ] || fail "step 13 left its extraction directory behind"
docker_log_has "^run image=[^ ]*ensembl-vep.* --cache_version ${VEP_CACHE_RELEASE} " "step 13 did not pass --cache_version ${VEP_CACHE_RELEASE} to vep"

# --- 6. pypgx bundle of another tag ---------------------------------------------------------
BUNDLE="${G}/reference/pypgx-bundle"
mkdir -p "$BUNDLE"
bgit() { git -C "$BUNDLE" -c user.name=case -c user.email=case@example.invalid -c init.defaultBranch=main "$@"; }
bgit init -q
bgit commit -q --allow-empty -m bundle
bgit tag 0.0.1
: > "$FAKE_DOCKER_LOG"
run_rc pypgx-old "${SCRIPTS}/32-pypgx.sh" sample1
[ "$RC" -ne 0 ] || fail "step 32 exited 0 with a pypgx bundle of another tag"
output_has pypgx-old "'0\.0\.1'"
output_has pypgx-old "needs ${PYPGX_BUNDLE_VERSION//./\\.}"
if awk '/^run / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 32 started pypgx with a bundle of another tag"
fi
bgit tag -d 0.0.1 > /dev/null
bgit tag "$PYPGX_BUNDLE_VERSION"
run_rc pypgx-right "${SCRIPTS}/32-pypgx.sh" sample1
docker_log_has '^run image=[^ ]*pypgx' "step 32 did not run pypgx with the bundle tag versions.env names"
