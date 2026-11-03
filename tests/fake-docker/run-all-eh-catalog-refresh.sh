#!/usr/bin/env bash
# run-all.sh with the bundled ExpansionHunter catalog after an image bump. The
# copy it extracts from the image is named after the image tag, so a copy an
# older image left (here under the old untagged name and under another tag)
# is never passed to the pipeline: the run extracts the catalog of this
# image once and passes that file, and a second run reuses it.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true TOOLS=expansion_hunter EH_CATALOG=bundled
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"
TAG=${EXPANSIONHUNTER_IMAGE##*:}
FRESH="${G}/reference/expansionhunter_variant_catalog-${TAG}.json"
# What older images left: the untagged copy and one of another tag.
echo '[{"LocusId": "OLD"}]' > "${G}/reference/expansionhunter_variant_catalog.json"
echo '[{"LocusId": "OLD"}]' > "${G}/reference/expansionhunter_variant_catalog-0.0.0--old_0.json"

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
case "${*:2}" in
  "cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json") echo '[{"LocusId": "HTT"}]' ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

extracts() { grep -F "run image=${EXPANSIONHUNTER_IMAGE} :: " "$FAKE_DOCKER_LOG" | grep -c 'variant_catalog/grch38/variant_catalog.json' || true; }

run_expect 0 run-all "${SCRIPTS}/run-all.sh" sample1 male
output_lacks run-all 'unbound variable'
[ "$(extracts)" -eq 1 ] || fail "the catalog was extracted from ${EXPANSIONHUNTER_IMAGE} $(extracts) times, expected once"
[ "$(cat "$FRESH" 2>/dev/null)" = '[{"LocusId": "HTT"}]' ] || fail "no fresh catalog at ${FRESH}"
docker_log_has "^nextflow :: .* --expansion_catalog ${FRESH} " "nextflow was not given this image's catalog ${FRESH}"

run_expect 0 run-all-again "${SCRIPTS}/run-all.sh" sample1 male
[ "$(extracts)" -eq 1 ] || fail "the second run extracted the catalog again ($(extracts) extractions in all)"
echo "An image bump gives run-all.sh a fresh bundled ExpansionHunter catalog, extracted once per image tag."
