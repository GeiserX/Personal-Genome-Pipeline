#!/usr/bin/env bash
# e2e-run.sh — run every end-to-end case in tests/e2e/ on the HG002 fixture.
#
# Usage: scripts/ci/e2e-run.sh [pattern]
#        scripts/ci/e2e-run.sh --self-test-pull
#   pattern  optional shell glob; runs only the cases whose file name matches
#            (e.g. '2*' or '*cyrius*'). Later cases read earlier cases' outputs,
#            so a partial run is for debugging only.
#   --self-test-pull  checks the pull retry of step 3 against a fake docker
#            (seconds, no network) and exits.
#
# Env: E2E_WORK  work area (default ${RUNNER_TEMP:-/tmp}/e2e-work), about 25 GB
#      GH_TOKEN  token for `gh release download` (the job's GITHUB_TOKEN)
#
# What it does:
#   1. downloads the release named in tests/fixtures/VERSION and checks
#      SHA256SUMS; when a pull request bumps VERSION, waits up to 75 minutes
#      for the build-fixture job to publish it;
#   2. lays out GENOME_DIR the way setup.sh and step 13 would leave it;
#   3. pulls every image the cases use that is not on the machine yet,
#      trying a pull again only when the registry answered with a 5xx;
#   4. runs tests/e2e/*.sh in name order (C locale: numbered cases first, then
#      the <package-key>-*.sh cases later packages add) and keeps going after a
#      failure, so one run lists every broken step;
#   5. prints a case / result / time / log table to the job summary and exits
#      1 if any case failed.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
PATTERN=${1:-*}

export REPO
export SAMPLE=HG002
export THREADS=4
export E2E_WORK="${E2E_WORK:-${RUNNER_TEMP:-/tmp}/e2e-work}"
export FIXTURE_DIR="${E2E_WORK}/fixture"
export GENOME_DIR="${E2E_WORK}/genome"
export E2E_NOTES="${E2E_WORK}/notes.md"
LOG_DIR="${E2E_WORK}/logs"
CASE_TIMEOUT=${CASE_TIMEOUT:-3600}
GH_REPO=${GITHUB_REPOSITORY:-GeiserX/Personal-Genome-Pipeline}
TAG=$(tr -d '[:space:]' < "${REPO}/tests/fixtures/VERSION")
SUMMARY=${GITHUB_STEP_SUMMARY:-/dev/null}

# The docker shim clamps --cpus to this machine's CPU count.
export PATH="${REPO}/tests/e2e/bin:${PATH}"

# --- Image pulls ----------------------------------------------------------------
# A `docker run` whose image is missing pulls it with one manifest request, and
# a registry 5xx on that request ends the step: in run 37980248545 Docker Hub
# answered 500 for hap.py and the GIAB case failed; the rerun passed. So the
# images are pulled before the cases, and a pull that fails with a registry
# 5xx is tried again, up to PULL_TRIES tries in all. Any other error (a wrong
# tag, a rate limit) fails at once. The tools' own `docker run` is never retried.
PULL_TRIES=${PULL_TRIES:-4}
PULL_WAIT=${PULL_WAIT:-20}   # seconds before the second try; doubles after each
REGISTRY_5XX='(HTTP status|status code):? 5[0-9][0-9]|50[0-4] (Internal Server Error|Bad Gateway|Service Unavailable|Gateway Time-?out)'
# pull_image IMAGE: docker pull with the retry above.
pull_image() {
  local image=$1 try=1 wait=$PULL_WAIT out
  while :; do
    if out=$(docker pull -q "$image" 2>&1); then echo "pulled ${image}"; return 0; fi
    echo "$out" >&2
    if ! grep -qE "$REGISTRY_5XX" <<< "$out"; then
      echo "ERROR: pulling ${image} failed, and not with a registry 5xx: not tried again." >&2
      return 1
    fi
    if [ "$try" -ge "$PULL_TRIES" ]; then
      echo "ERROR: pulling ${image}: a registry 5xx on each of ${try} tries." >&2
      return 1
    fi
    echo "Registry 5xx pulling ${image} (try ${try} of ${PULL_TRIES}); trying again in ${wait}s." >&2
    sleep "$wait"
    try=$((try + 1)) wait=$((wait * 2))
  done
}
# Images no case uses: those of the steps in docs/testing.md, "What no e2e case
# runs", plus BWA_IMAGE (the classic index GRIDSS needs, step 04b), PICARD_IMAGE
# (chip-to-vcf) and SAMTOOLS_HTTPS_IMAGE (docs/quick-test.md). Pulling them
# would only cost time and cache space.
UNUSED_IMAGES=" SAMTOOLS_HTTPS_IMAGE BWA_IMAGE STRELKA_IMAGE OCTOPUS_IMAGE CLAIR3_IMAGE TIDDIT_IMAGE SNIFFLES_IMAGE GRIDSS_IMAGE DUPHOLD_IMAGE CNVPYTOR_IMAGE ANNOTSV_IMAGE VEP_IMAGE PCGR_IMAGE PICARD_IMAGE "
# case_images: the NAME_IMAGE="ref" lines of versions.env and of the cases,
# without the images above, once each.
case_images() {
  grep -hoE '^[A-Z0-9_]+_IMAGE="[^"]+"' "${REPO}/versions.env" "${REPO}"/tests/e2e/*.sh \
    | while IFS='=' read -r name ref; do
        [[ "$UNUSED_IMAGES" == *" ${name} "* ]] || echo "${ref//\"/}"
      done | awk '!seen[$0]++'
}

if [ "${1:-}" = --self-test-pull ]; then
  T=$(mktemp -d)
  trap 'rm -rf "$T"' EXIT
  # A fake docker that only pulls: each image name has one behaviour, and every
  # pull is counted in $T/<image>.count.
  cat > "${T}/docker" <<'FAKE'
#!/usr/bin/env bash
[ "$1" = pull ] || { echo "fake docker: only pull is expected, got: $*" >&2; exit 2; }
img=${*: -1}
c="${FAKE_DIR}/${img//[\/:]/_}.count"
echo x >> "$c"
n=$(wc -l < "$c")
case "$img" in
  flaky500)
    [ "$n" -ge 2 ] && exit 0
    echo 'Error response from daemon: Head "https://registry-1.docker.io/v2/library/flaky500/manifests/latest": received unexpected HTTP status: 500 Internal Server Error' >&2 ;;
  down503)
    echo 'Error response from daemon: Head "https://registry-1.docker.io/v2/library/down503/manifests/latest": received unexpected HTTP status: 503 Service Unavailable' >&2 ;;
  notag)
    echo 'Error response from daemon: manifest for notag not found: manifest unknown: manifest unknown' >&2 ;;
  ratelimit)
    echo 'Error response from daemon: toomanyrequests: You have reached your unauthenticated pull rate limit. https://www.docker.com/increase-rate-limit' >&2 ;;
esac
exit 1
FAKE
  chmod +x "${T}/docker"
  export FAKE_DIR="$T"
  PATH="${T}:${PATH}" PULL_WAIT=0
  st_failed=0
  # expect DESC pass|fail PULLS IMAGE
  expect() {
    local rc=0 got n
    pull_image "$4" > /dev/null 2> "${T}/err" || rc=$?
    got=$([ "$rc" -eq 0 ] && echo pass || echo fail)
    n=$({ cat "${T}/${4//[\/:]/_}.count" 2> /dev/null || true; } | wc -l | tr -d ' ')
    if [ "$got" = "$2" ] && [ "$n" = "$3" ]; then
      echo "ok: $1"
    else
      echo "FAIL: $1 (${got} after ${n} pulls; want ${2} after ${3})"
      sed 's/^/  /' "${T}/err"
      st_failed=1
    fi
  }
  expect "a registry 500 on the first pull is tried again and passes" pass 2 flaky500
  expect "a registry 503 on every pull gives up after ${PULL_TRIES} tries" fail "$PULL_TRIES" down503
  expect "a missing tag fails at once" fail 1 notag
  expect "a rate limit (429) fails at once" fail 1 ratelimit
  LIST=$(case_images)
  # shellcheck source=../../versions.env
  HAPPY=$(. "${REPO}/versions.env" && echo "$HAPPY_IMAGE") VEP=$(. "${REPO}/versions.env" && echo "$VEP_IMAGE")
  if grep -qxF "$HAPPY" <<< "$LIST"; then echo "ok: the pull list holds hap.py's image"; else echo "FAIL: the pull list lacks ${HAPPY}"; st_failed=1; fi
  if grep -qxF "$VEP" <<< "$LIST"; then echo "FAIL: the pull list holds ${VEP}, which no case runs"; st_failed=1; else echo "ok: the pull list leaves VEP out"; fi
  exit "$st_failed"
fi

mkdir -p "$E2E_WORK" "$FIXTURE_DIR" "$LOG_DIR"
: > "$E2E_NOTES"

fixture_complete() {
  [ -f "${FIXTURE_DIR}/SHA256SUMS" ] && (cd "$FIXTURE_DIR" && sha256sum --quiet -c SHA256SUMS)
}

# --- 1. Fixture ---------------------------------------------------------------
echo "=== Fixture ${TAG} ==="
deadline=$(( $(date +%s) + 75 * 60 ))
until fixture_complete; do
  if gh release view "$TAG" -R "$GH_REPO" >/dev/null 2>&1; then
    # Inside the condition, so a failed download (assets still uploading, a
    # network error) waits for the next try instead of ending the run.
    if gh release download "$TAG" -R "$GH_REPO" -D "$FIXTURE_DIR" --clobber && fixture_complete; then
      break
    fi
    echo "Release ${TAG} could not be downloaded, is incomplete or its checksums do not match."
  else
    echo "Release ${TAG} does not exist yet; the build-fixture job publishes it."
  fi
  if [ "$(date +%s)" -ge "$deadline" ]; then
    echo "ERROR: fixture ${TAG} not usable after 75 minutes. Run the E2E workflow with job=build-fixture." >&2
    exit 1
  fi
  sleep 60
done
(cd "$FIXTURE_DIR" && sha256sum -c SHA256SUMS)

# --- 2. GENOME_DIR layout -----------------------------------------------------
# The reference under the name scripts/lib/common.sh reads by default.
echo "=== Layout ${GENOME_DIR} ==="
mkdir -p "${GENOME_DIR}/reference" "${GENOME_DIR}/clinvar" "${GENOME_DIR}/annotations" \
  "${GENOME_DIR}/${SAMPLE}/fastq" "${GENOME_DIR}/${SAMPLE}/vep"
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set"
if [ ! -s "${REF}.fasta" ]; then
  gzip -dc "${FIXTURE_DIR}/fixture_ref.fa.gz" > "${REF}.fasta.tmp"
  mv "${REF}.fasta.tmp" "${REF}.fasta"
fi
cp "${FIXTURE_DIR}/fixture_ref.fa.gz.fai" "${REF}.fasta.fai"
cp "${FIXTURE_DIR}/fixture_ref.dict" "${REF}.dict"
for f in clinvar.vcf.gz clinvar_chr.vcf.gz clinvar_pathogenic_chr.vcf.gz; do
  cp "${FIXTURE_DIR}/${f}" "${FIXTURE_DIR}/${f}.tbi" "${GENOME_DIR}/clinvar/"
done
# The synthetic score file stands in for REVEL (step 30 looks for this name).
cp "${FIXTURE_DIR}/revel_synthetic.tsv.gz" "${GENOME_DIR}/annotations/revel_grch38.tsv.gz"
cp "${FIXTURE_DIR}/revel_synthetic.tsv.gz.tbi" "${GENOME_DIR}/annotations/revel_grch38.tsv.gz.tbi"
# Reads for step 02, and the VEP output step 13 would have written.
cp "${FIXTURE_DIR}/${SAMPLE}_R1.fastq.gz" "${FIXTURE_DIR}/${SAMPLE}_R2.fastq.gz" "${GENOME_DIR}/${SAMPLE}/fastq/"
cp "${FIXTURE_DIR}/${SAMPLE}_vep.vcf" "${GENOME_DIR}/${SAMPLE}/vep/"
df -h "$E2E_WORK"

# --- 3. Images ----------------------------------------------------------------
# Only a full run pre-pulls. A case reaches most images through the steps it
# runs, so the images of a few cases cannot be told apart from the rest; a
# partial run lets its cases pull their own, without the retry.
if [ "$PATTERN" = "*" ]; then
  echo "=== Pulling the images the cases use ==="
  mapfile -t IMAGES < <(case_images)
  [ "${#IMAGES[@]}" -gt 0 ] || { echo "ERROR: no NAME_IMAGE line found in versions.env" >&2; exit 1; }
  pulled=0
  for image in "${IMAGES[@]}"; do
    docker image inspect "$image" > /dev/null 2>&1 && continue
    pull_image "$image" || exit 1
    pulled=$((pulled + 1))
  done
  echo "${#IMAGES[@]} images: $(( ${#IMAGES[@]} - pulled )) already here, ${pulled} pulled"
else
  echo "=== Partial run ('${PATTERN}'): no pre-pull, each case pulls what it uses ==="
fi

# --- 4. Cases -----------------------------------------------------------------
mapfile -t CASES < <(cd "${REPO}/tests/e2e" && find . -maxdepth 1 -type f -name '*.sh' ! -name lib.sh -printf '%f\n' | LC_ALL=C sort)
names=() results=() times=()
failed=0
for c in "${CASES[@]}"; do
  # shellcheck disable=SC2053  # PATTERN is a glob on purpose
  [[ "$c" == $PATTERN ]] || continue
  name="${c%.sh}"
  log="${LOG_DIR}/${name}.log"
  start=$(date +%s)
  echo "::group::${name}"
  set +e
  timeout -k 30 "$CASE_TIMEOUT" bash "${REPO}/tests/e2e/${c}" 2>&1 | tee "$log"
  rc=${PIPESTATUS[0]}
  set -e
  echo "::endgroup::"
  case "$rc" in
    0)   result=pass ;;
    124) result="FAIL (timeout ${CASE_TIMEOUT}s)"; failed=$((failed + 1)) ;;
    *)   result="FAIL (exit ${rc})"; failed=$((failed + 1)) ;;
  esac
  names+=("$name"); results+=("$result"); times+=("$(( $(date +%s) - start ))s")
  echo "${name}: ${result}"
done
if [ "${#names[@]}" -eq 0 ]; then
  echo "ERROR: no case matches '${PATTERN}'" >&2
  exit 1
fi
# Tells the workflow every case ran, so the pulled images are complete enough to cache.
[ "$PATTERN" = "*" ] && touch "${E2E_WORK}/all-cases-ran"

# --- 5. Report ----------------------------------------------------------------
{
  echo "### E2E on fixture ${TAG}: ${failed} of ${#names[@]} cases failed"
  echo
  echo "| Case | Result | Time | Log |"
  echo "|---|---|---|---|"
  for i in "${!names[@]}"; do
    echo "| ${names[$i]} | ${results[$i]} | ${times[$i]} | logs/${names[$i]}.log |"
  done
  echo
  echo "Logs are in the e2e-logs artifact of this run."
  for i in "${!names[@]}"; do
    [ "${results[$i]}" = pass ] && continue
    echo
    echo "#### ${names[$i]}"
    echo '```'
    grep -E '^\[FAIL\]' "${LOG_DIR}/${names[$i]}.log" || tail -n 15 "${LOG_DIR}/${names[$i]}.log"
    echo '```'
  done
  if [ -s "$E2E_NOTES" ]; then
    echo
    cat "$E2E_NOTES"
  fi
} | tee -a "$SUMMARY"

[ "$failed" -eq 0 ]
