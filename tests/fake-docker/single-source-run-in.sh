#!/usr/bin/env bash
# The docker wrapper every script uses (run_in, scripts/lib/common.sh):
#   - no script starts a container any other way;
#   - an analysis step runs with no network, the data directory read-only,
#     only the sample directory writable, as the calling user;
#   - a step that says --net gets the network and still runs as the caller;
#   - a write outside the sample directory fails unless the call says --rw;
#   - the image is the one versions.env names: changing that one line changes
#     what the script runs, and without versions.env the script does not start.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook

# --- every container goes through run_in ---------------------------------------
# A `docker run` (or `"$CONTAINER_ENGINE" run`) outside a printed hint would
# start a container without the defaults below.
if bypass=$(grep -nE '^[^#]*(^|[^[:alnum:]_])(docker|\$\{?CONTAINER_ENGINE\}?"?) run([^[:alnum:]_-]|$)' "${SCRIPTS}"/*.sh | grep -v 'echo '); then
  fail "these lines start a container without run_in: ${bypass}"
fi

# --- an analysis step ---------------------------------------------------------
run_expect 0 roh "${SCRIPTS}/11-roh-analysis.sh" sample1
docker_log_has "^run image=[^ ]*bcftools[^ ]* :: .*--network none .*--user [0-9]+:[0-9]+ -e HOME=/tmp -v [^ ]*/genome:/genome:ro -v [^ ]*/genome/sample1:/genome/sample1 " \
  "step 11 did not run bcftools with --network none, the caller's user, /genome read-only and the sample directory writable"
[ -f "${GENOME_DIR}/sample1/vcf/sample1_roh.txt" ] || fail "step 11 wrote no output (the hook did not run)"

# --- a step that downloads and installs (Cyrius: pip) -------------------------
: > "$FAKE_DOCKER_LOG"
run_rc cyrius "${SCRIPTS}/21-cyrius.sh" sample1
docker_log_has '^run image=[^ ]*python' "step 21 never ran the python image"
if awk '/^run image=[^ ]*python/ && (/--network none/ || !/--user [0-9]+:[0-9]+ /) { bad = 1 } END { exit !bad }' "$FAKE_DOCKER_LOG"; then
  fail "step 21 (pip install) ran without network or not as the calling user; it needs --net and no --root"
fi

# --- read-only data directory ---------------------------------------------------
cat > "${CASE_WORK}/write-ref.sh" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
SAMPLE=sample1
. "${REPO_ROOT}/scripts/lib/common.sh"
# shellcheck disable=SC2086  # RW_OPT is empty or "--rw DIR"
run_in ${RW_OPT:-} "$BCFTOOLS_IMAGE" bcftools view -Oz -o /genome/reference/shared.vcf.gz /genome/sample1/vcf/sample1.vcf.gz
SH
chmod +x "${CASE_WORK}/write-ref.sh"
run_rc write-ro "${CASE_WORK}/write-ref.sh"
[ "$RC" -ne 0 ] || fail "a container wrote into reference/ although the data directory is mounted read-only"
output_has write-ro 'Read-only file system'
[ ! -e "${GENOME_DIR}/reference/shared.vcf.gz" ] || fail "reference/shared.vcf.gz was written through a read-only mount"
RW_OPT="--rw ${GENOME_DIR}/reference" run_expect 0 write-rw "${CASE_WORK}/write-ref.sh"
[ -f "${GENOME_DIR}/reference/shared.vcf.gz" ] || fail "run_in --rw did not make reference/ writable"

# A path outside GENOME_DIR cannot be made writable: containers never see it.
RW_OPT="--rw ${CASE_WORK}/elsewhere" run_rc write-outside "${CASE_WORK}/write-ref.sh"
[ "$RC" -eq 2 ] || fail "run_in --rw accepted a directory outside GENOME_DIR (exit ${RC}, expected 2)"

# --- one line in versions.env decides the image ---------------------------------
COPY="${CASE_WORK}/repo"
mkdir -p "$COPY"
cp -R "${REPO_ROOT}/scripts" "${REPO_ROOT}/versions.env" "$COPY/"
# shellcheck source=../../versions.env
OLD=$(. "${REPO_ROOT}/versions.env" && printf '%s' "$BCFTOOLS_IMAGE")
sed 's|^BCFTOOLS_IMAGE=.*|BCFTOOLS_IMAGE="example/bcftools:9.99"|' "${REPO_ROOT}/versions.env" > "${COPY}/versions.env"
: > "$FAKE_DOCKER_LOG"
rm -f "${GENOME_DIR}/sample1/vcf/sample1_roh.txt"
run_expect 0 bumped "${COPY}/scripts/11-roh-analysis.sh" sample1
docker_log_has '^run image=example/bcftools:9\.99 ' "step 11 did not run the image the edited versions.env names"
if grep -qF "image=${OLD} " "$FAKE_DOCKER_LOG"; then
  fail "step 11 still ran ${OLD} after versions.env was changed"
fi

# --- no versions.env, no run ------------------------------------------------------
rm "${COPY}/versions.env"
: > "$FAKE_DOCKER_LOG"
run_rc no-versions "${COPY}/scripts/11-roh-analysis.sh" sample1
[ "$RC" -ne 0 ] || fail "step 11 ran without versions.env"
if awk '/^run / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 11 started a container without versions.env (a fallback image is hidden somewhere)"
fi
