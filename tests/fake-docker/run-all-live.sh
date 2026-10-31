#!/usr/bin/env bash
# run-all.sh and the live site (docs/live-site.md):
#   1. with PGP_LIVE_URL set, nextflow starts with PGP_LIVE_RUN=sample1 and
#      -c <sample>/nextflow/live.config, which loads the pinned nf-weblog
#      plugin and points the weblog at PGP_LIVE_URL/ingest/sample1,
#      the script step 04b (GRIDSS=true) posts 'started' and its result, and
#      the run ends with a pgp_run_finished message carrying the exit code;
#   2. without PGP_LIVE_URL, nextflow sees no PGP_LIVE_RUN and no live.config,
#      and nothing is posted; a URL that is not a plain http(s) URL exits 2;
#   3. with the site refusing every post, the exit code is the same as in 1 and
#      one warning says so.
# Wrappers in front of the fake nextflow and the fake curl log what they see.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
use_output_hook
seed_sample "$G" sample1

mkdir -p "${CASE_WORK}/bin"
cat > "${CASE_WORK}/bin/nextflow" <<'NF'
#!/usr/bin/env bash
echo "nextflow PGP_LIVE_RUN=${PGP_LIVE_RUN:-unset} PGP_LIVE_URL=${PGP_LIVE_URL:-unset}" >> "${CASE_WORK}/live.log"
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
cat > "${CASE_WORK}/bin/curl" <<'CURL'
#!/usr/bin/env bash
for a in "$@"; do case "$a" in */ingest/*) echo "post ${a} $(cat)" >> "${CASE_WORK}/live.log" ;; esac; done
exec "${REPO_ROOT}/scripts/ci/fake-docker/curl" "$@"
CURL
chmod +x "${CASE_WORK}/bin/nextflow" "${CASE_WORK}/bin/curl"
P="${CASE_WORK}/bin:${PATH}"
LOG="${CASE_WORK}/live.log"
URL=http://127.0.0.1:8788

# --- 1. the site is set ----------------------------------------------------------------------
: > "$LOG"
# GRIDSS fails on the fake docker, so run-all.sh exits 1; what counts is what it posted.
run_rc on env PATH="$P" GRIDSS=true PGP_LIVE_URL="$URL" "${SCRIPTS}/run-all.sh" sample1 male
ON_RC=$RC
[ "$ON_RC" -le 1 ] || fail "run-all.sh exited ${ON_RC}, expected 0 or 1"
grep -q "^nextflow PGP_LIVE_RUN=sample1 PGP_LIVE_URL=${URL}$" "$LOG" || fail "nextflow did not see PGP_LIVE_RUN=sample1: $(cat "$LOG")"
grep -q "^post ${URL}/ingest/sample1 {\"event\":\"pgp_step\",\"step\":\"04b\",\"state\":\"started\"," "$LOG" \
  || fail "no 'started' post for step 04b: $(cat "$LOG")"
grep -qE "^post ${URL}/ingest/sample1 \{\"event\":\"pgp_step\",\"step\":\"04b\",\"state\":\"(ok|failed)\"," "$LOG" \
  || fail "no result post for step 04b: $(cat "$LOG")"
grep -q "^post ${URL}/ingest/sample1 {\"event\":\"pgp_run_finished\",\"exit\":${ON_RC}," "$LOG" \
  || fail "no pgp_run_finished post with exit ${ON_RC}: $(cat "$LOG")"
output_has on "Live site: messages go to ${URL}/ingest/sample1"
LC="${G}/sample1/nextflow/live.config"
grep -q '^nextflow :: .* -c '"$LC" "$FAKE_DOCKER_LOG" || fail "nextflow was not given -c ${LC}: $(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1)"
grep -qE "^    id 'nf-weblog@[0-9]+\.[0-9]+\.[0-9]+'$" "$LC" || fail "live.config does not load a pinned nf-weblog: $(cat "$LC")"
grep -q "^    url     = '${URL}/ingest/sample1'$" "$LC" || fail "live.config does not point at the site: $(cat "$LC")"
# Every post is a message the site takes (bin/pgp_live.py's cut())
sed -n 's/^post [^ ]* //p' "$LOG" > "${CASE_WORK}/bodies.jsonl"
python3 -c 'import sys; sys.path.insert(0, sys.argv[1] + "/bin"); import json, pgp_live
[pgp_live.cut(json.loads(l)) for l in open(sys.argv[2])]' "$REPO_ROOT" "${CASE_WORK}/bodies.jsonl" \
  || fail "a post is not a message the site takes: $(cat "${CASE_WORK}/bodies.jsonl")"

# --- 2. no site --------------------------------------------------------------------------------
: > "$LOG"
run_rc off env PATH="$P" GRIDSS=true "${SCRIPTS}/run-all.sh" sample1 male
grep -q '^nextflow PGP_LIVE_RUN=unset PGP_LIVE_URL=unset$' "$LOG" || fail "nextflow saw a live setting without PGP_LIVE_URL: $(cat "$LOG")"
if grep -q '^post ' "$LOG"; then fail "run-all.sh posted without PGP_LIVE_URL: $(cat "$LOG")"; fi
output_lacks off 'Live site'
if grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1 | grep -q 'live\.config'; then fail "nextflow was given live.config without PGP_LIVE_URL"; fi
[ ! -e "${G}/sample1/nextflow/live.config" ] || fail "an earlier run's live.config was left in place"
run_expect 2 badurl env PATH="$P" PGP_LIVE_URL="http://127.0.0.1:8788/x'y" "${SCRIPTS}/run-all.sh" sample1 male
output_has badurl 'PGP_LIVE_URL must be a plain'

# --- 3. the site refuses every post -------------------------------------------------------------
run_rc down env PATH="$P" GRIDSS=true PGP_LIVE_URL="$URL" FAKE_DOWNLOAD_FAIL='/ingest/' "${SCRIPTS}/run-all.sh" sample1 male
[ "$RC" -eq "$ON_RC" ] || fail "with the site down run-all.sh exited ${RC}, with it up ${ON_RC}"
output_has down 'WARNING: the live site at http://127.0.0.1:8788 did not take a message'
echo "run-all.sh posts its steps when PGP_LIVE_URL is set, nothing without it, and a site that is down changes nothing."
