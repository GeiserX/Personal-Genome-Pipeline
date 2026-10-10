#!/usr/bin/env bash
# setup.sh on an empty GENOME_DIR: the reference and ClinVar downloads go
# through the fake wget, the image checks through the fake docker. It must
# finish with exit 0.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

G="${CASE_WORK}/genome"
mkdir -p "$G"

run_expect 0 setup "${SCRIPTS}/setup.sh" "$G"
output_lacks setup 'unbound variable'
output_has setup 'Setup complete!'

for f in reference/GRCh38_no_alt_analysis_set.fasta reference/GRCh38_no_alt_analysis_set.fasta.fai \
         clinvar/clinvar.vcf.gz clinvar/clinvar.vcf.gz.tbi; do
  [ -s "${G}/${f}" ] || fail "setup.sh did not download ${f}"
done
grep -q '^image inspect' "$FAKE_DOCKER_LOG" || fail "setup.sh never checked a docker image"

# --- a failed pull shows Docker's message and is tried again ------------------------
# A docker in front of the fake one: `pull` of the DeepVariant image fails with
# Docker Hub's rate-limit message the first PULL_FAILS times, then goes through.
mkdir -p "${CASE_WORK}/pull-bin"
cat > "${CASE_WORK}/pull-bin/docker" <<'SHIM'
#!/usr/bin/env bash
if [ "${1:-}" = pull ] && [[ "${2:-}" == *deepvariant* ]]; then
  n=$(( $(cat "${CASE_WORK}/pull-count" 2>/dev/null || echo 0) + 1 ))
  echo "$n" > "${CASE_WORK}/pull-count"
  if [ "$n" -le "${PULL_FAILS:-0}" ]; then
    echo "Error response from daemon: toomanyrequests: You have reached your unauthenticated pull rate limit." >&2
    exit 1
  fi
fi
exec "${REPO_ROOT}/scripts/ci/fake-docker/docker" "$@"
SHIM
chmod +x "${CASE_WORK}/pull-bin/docker"

# Twice refused, the third try goes through: setup succeeds and says why it retried.
rm -f "${CASE_WORK}/pull-count"
PATH="${CASE_WORK}/pull-bin:${PATH}" PULL_FAILS=2 FAKE_DOCKER_MISSING_IMAGES=deepvariant \
  run_expect 0 pull-retry "${SCRIPTS}/setup.sh" --pull-only
output_has pull-retry 'toomanyrequests: You have reached your unauthenticated pull rate limit'
output_has pull-retry 'Pull attempt 1/3 failed: [^ ]*deepvariant'
output_has pull-retry 'Pull attempt 2/3 failed: [^ ]*deepvariant'
output_lacks pull-retry 'Pull attempt 3/3'
output_has pull-retry '\[OK\] Docker images: 1 pulled, [0-9]+ already present, 0 failed'
[ "$(cat "${CASE_WORK}/pull-count")" -eq 3 ] || fail "setup.sh pulled the DeepVariant image $(cat "${CASE_WORK}/pull-count") times, expected 3"

# Refused every time: setup fails after FETCH_TRIES tries, with Docker's message.
rm -f "${CASE_WORK}/pull-count"
PATH="${CASE_WORK}/pull-bin:${PATH}" PULL_FAILS=99 FETCH_TRIES=2 FAKE_DOCKER_MISSING_IMAGES=deepvariant \
  run_expect 1 pull-fail "${SCRIPTS}/setup.sh" --pull-only
output_has pull-fail 'toomanyrequests'
output_has pull-fail 'Pull attempt 2/2 failed'
output_has pull-fail "WARNING: Failed to pull [^ ]*deepvariant[^ ]*; Docker's message is above"
output_lacks pull-fail 'Check the image name/tag'
[ "$(cat "${CASE_WORK}/pull-count")" -eq 2 ] || fail "setup.sh pulled $(cat "${CASE_WORK}/pull-count") times with FETCH_TRIES=2"
