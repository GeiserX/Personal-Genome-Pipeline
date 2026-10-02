#!/usr/bin/env bash
# setup.sh and validate-setup.sh take their image list from versions.env:
# a new *_IMAGE line is pulled and checked with no other edit, a line marked
# `# optional` is left for the step that uses it, and `setup.sh --pull-only`
# pulls the list without a data directory and without downloading anything.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

COPY="${CASE_WORK}/repo"
mkdir -p "$COPY"
cp -R "${REPO_ROOT}/scripts" "${REPO_ROOT}/versions.env" "$COPY/"
cat >> "${COPY}/versions.env" <<'ENV'
FOO_IMAGE="example/foo:1.0"
BAR_IMAGE="example/bar:2.0"  # optional: a tool no default step runs
ENV

# Nothing is pulled yet, so setup has to pull every image it lists.
export FAKE_DOCKER_MISSING_IMAGES='.'
run_expect 0 pull-only "${COPY}/scripts/setup.sh" --pull-only
docker_log_has '^pull :: pull example/foo:1\.0 ' "setup.sh --pull-only did not pull the image added to versions.env"
if grep -q 'example/bar' "$FAKE_DOCKER_LOG"; then
  fail "setup.sh --pull-only pulled an image marked '# optional'"
fi
if grep -qE '^(curl|wget|run) ' "$FAKE_DOCKER_LOG"; then
  fail "setup.sh --pull-only downloaded data or started a container"
fi
# Every image the list holds is a line of versions.env, and nothing else.
want=$(grep -E '^[A-Z0-9_]+_IMAGE=' "${COPY}/versions.env" | grep -vc '# optional')
got=$(grep -c '^pull :: ' "$FAKE_DOCKER_LOG")
[ "$got" -eq "$want" ] || fail "setup.sh --pull-only pulled ${got} images, versions.env lists ${want} that are not optional"
[ "$want" -ge 25 ] || fail "only ${want} images counted in versions.env; the count is broken"

# validate-setup.sh checks the same list.
unset FAKE_DOCKER_MISSING_IMAGES
export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
: > "$FAKE_DOCKER_LOG"
run_expect 0 validate "${COPY}/scripts/validate-setup.sh"
docker_log_has '^image inspect :: image inspect example/foo:1\.0 ' "validate-setup.sh did not check the image added to versions.env"
output_has validate "All ${want} Docker images are pulled"
if grep -q 'example/bar' "$FAKE_DOCKER_LOG"; then
  fail "validate-setup.sh checked an image marked '# optional'"
fi

# A missing image is reported by name.
: > "$FAKE_DOCKER_LOG"
FAKE_DOCKER_MISSING_IMAGES='example/foo' run_expect 1 validate-missing "${COPY}/scripts/validate-setup.sh"
output_has validate-missing 'docker pull example/foo:1\.0'
