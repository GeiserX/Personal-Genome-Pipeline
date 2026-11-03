#!/usr/bin/env bash
# e2e-run.sh — run every end-to-end case in tests/e2e/ on the HG002 fixture.
#
# Usage: scripts/ci/e2e-run.sh [pattern]
#        scripts/ci/e2e-run.sh --self-test-pull | --image-groups | --pull-images
#   pattern  optional shell glob; runs only the cases whose file name matches
#            (e.g. '2*' or '*cyrius*'). Later cases read earlier cases' outputs,
#            so a partial run is for debugging only.
#   --self-test-pull  checks the pull retry of step 3 against a fake docker
#            (seconds, no network) and the image cache groups, and exits.
#   --image-groups  prints "GROUP IMAGE" for every image the cases use (the
#            e2e workflow's image cache entries) and exits.
#   --pull-images  pulls those images, as step 3 below does, and exits.
#
# Env: E2E_WORK  work area (default ${RUNNER_TEMP:-/tmp}/e2e-work), about 25 GB
#      GH_TOKEN  token for `gh release download` (the job's GITHUB_TOKEN)
#      E2E_CASES  a file of case names (without .sh): run only those, as a
#                full run (pre-pull, all-cases-ran marker). The e2e workflow
#                gives each of its parallel jobs its list from e2e-shards.sh.
#      E2E_JOB   the job name written beside each result (default local)
#      E2E_BASE_OUT  after the cases, pack every file they wrote or changed
#                (under E2E_WORK and the checkout) into this .tar.zst, and
#                write "E2E_WORK CHECKOUT" to E2E_BASE_OUT.from
#      E2E_BASE_IN   after the layout, unpack this .tar.zst from E2E_BASE_OUT,
#                which must have been packed with the same E2E_WORK and
#                checkout, then delete it to give the cases its disk space
#      E2E_TRACE=true  after each case, record how long it took and which
#                earlier case wrote each file it read (atime and ctime, with
#                the file systems remounted strictatime; needs sudo), in
#                $E2E_WORK/report/trace.tsv for `e2e-shards.sh times`
#
# Every run writes $E2E_WORK/report/results.tsv ("CASE RESULT SECONDS JOB",
# tab-separated) and report/extra.md (failed checks and notes), which the
# e2e workflow merges into one table across its jobs.
#
# What it does:
#   1. downloads the release named in tests/fixtures/VERSION and checks
#      SHA256SUMS; when a pull request bumps VERSION, waits up to 75 minutes
#      for the build-fixture job to publish it;
#   2. lays out GENOME_DIR the way setup.sh and step 13 would leave it;
#   3. pulls every image the cases use that is not on the machine yet,
#      trying a pull again only when the registry answered with a 5xx;
#   4. runs tests/e2e/*.sh in name order (C locale: numbered cases first, then
#      the <package-key>-*.sh cases later packages add), or those of E2E_CASES,
#      and keeps going after a failure, so one run lists every broken step;
#   5. prints a case / result / time / log table to the job summary;
#   6. with E2E_BASE_OUT, packs what the cases wrote; exits 1 if any case failed.
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
# runs", plus BWA_IMAGE (the classic index GRIDSS needs, step 04b) and
# SAMTOOLS_HTTPS_IMAGE (docs/quick-test.md). Pulling them would only cost time
# and cache space. PICARD_IMAGE is not one of them: consumer-export-1 runs
# chip-to-vcf.sh, whose liftover is Picard's.
UNUSED_IMAGES=" SAMTOOLS_HTTPS_IMAGE BWA_IMAGE STRELKA_IMAGE OCTOPUS_IMAGE CLAIR3_IMAGE TIDDIT_IMAGE SNIFFLES_IMAGE GRIDSS_IMAGE DUPHOLD_IMAGE CNVPYTOR_IMAGE ANNOTSV_IMAGE VEP_IMAGE PCGR_IMAGE "
# case_images: the NAME_IMAGE="ref" lines of versions.env and of the cases,
# without the images above, once each.
case_images() {
  grep -hoE '^[A-Z0-9_]+_IMAGE="[^"]+"' "${REPO}/versions.env" "${REPO}"/tests/e2e/*.sh \
    | while IFS='=' read -r name ref; do
        [[ "$UNUSED_IMAGES" == *" ${name} "* ]] || echo "${ref//\"/}"
      done | awk '!seen[$0]++'
}

# --- Image cache groups -------------------------------------------------------
# The e2e workflow caches the images of case_images in IMAGE_GROUPS cache
# entries, so no single entry, and not their sum, reaches the repository's
# 10 GB cache limit. An image's group comes from its repository name without
# the tag or digest, so a version bump moves no image to another group and
# changes the cache key of one group only.
IMAGE_GROUPS=4
# image_group REF: the group (0 to IMAGE_GROUPS-1) of an image reference.
image_group() {
  local repo=${1%@*}
  [[ "${repo##*/}" == *:* ]] && repo=${repo%:*}
  echo $(( $(printf '%s' "$repo" | cksum | cut -d' ' -f1) % IMAGE_GROUPS ))
}
# image_groups: one "GROUP REF" line per image of case_images.
image_groups() {
  local ref
  case_images | while IFS= read -r ref; do echo "$(image_group "$ref") ${ref}"; done
}

if [ "${1:-}" = --image-groups ]; then
  image_groups
  exit 0
fi

# pull_all: pull every image of case_images that is not here yet.
pull_all() {
  local image pulled=0
  mapfile -t IMAGES < <(case_images)
  [ "${#IMAGES[@]}" -gt 0 ] || { echo "ERROR: no NAME_IMAGE line found in versions.env" >&2; return 1; }
  for image in "${IMAGES[@]}"; do
    docker image inspect "$image" > /dev/null 2>&1 && continue
    pull_image "$image" || return 1
    pulled=$((pulled + 1))
  done
  echo "${#IMAGES[@]} images: $(( ${#IMAGES[@]} - pulled )) already here, ${pulled} pulled"
}

# The e2e workflow's shard legs pull while the base leg runs.
if [ "${1:-}" = --pull-images ]; then
  pull_all
  exit
fi

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
  # Cache groups: a version bump keeps an image in its group, and every image
  # of the pull list is in exactly one group.
  for pair in "example/tool:2.31--h118bc1c_0 example/tool:2.32--h0_0" \
              "example/pinned:3.15.0@sha256:16e7e5fc example/pinned:3.16.0@sha256:0123abcd" \
              "registry.example:5000/tool:1.0 registry.example:5000/tool:2.0"; do
    read -r a b <<< "$pair"
    if [ "$(image_group "$a")" = "$(image_group "$b")" ]; then echo "ok: ${a} and ${b} share a cache group"
    else echo "FAIL: ${a} and ${b} fall in different cache groups"; st_failed=1; fi
  done
  GROUPS_LIST=$(image_groups)
  if [ "$(cut -d' ' -f2 <<< "$GROUPS_LIST")" = "$LIST" ] && ! grep -qvE "^[0-9]+ " <<< "$GROUPS_LIST" \
     && [ -z "$(awk -v n="$IMAGE_GROUPS" '$1 >= n' <<< "$GROUPS_LIST")" ]; then
    echo "ok: each of the $(wc -l <<< "$LIST" | tr -d ' ') listed images is in one cache group of ${IMAGE_GROUPS}"
  else
    echo "FAIL: the cache groups do not hold the pull list once each"; st_failed=1
  fi
  exit "$st_failed"
fi

REPORT_DIR="${E2E_WORK}/report"
JOB=${E2E_JOB:-local}
TRACE=${E2E_TRACE:-false}
CASE_LIST=${E2E_CASES:-}
if [ -n "$CASE_LIST" ]; then
  [ -f "$CASE_LIST" ] || { echo "ERROR: E2E_CASES names no file: ${CASE_LIST}" >&2; exit 1; }
  [ "$PATTERN" = "*" ] || { echo "ERROR: give E2E_CASES or a pattern, not both" >&2; exit 1; }
fi
mkdir -p "$E2E_WORK" "$FIXTURE_DIR" "$LOG_DIR" "$REPORT_DIR"
: > "$E2E_NOTES"
: > "${REPORT_DIR}/results.tsv"
: > "${REPORT_DIR}/extra.md"

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
# Files whose ctime is newer than this were written by the cases (or unpacked
# from the base job), not by the layout above.
touch "${E2E_WORK}/.layout-done"

# The files the cases of the base job wrote, at the paths they were written to.
# base-from names the work area and checkout of the packing run; anywhere else
# the paths in the archive would land beside the files the cases read.
if [ -n "${E2E_BASE_IN:-}" ]; then
  echo "=== Unpacking the base cases' outputs ==="
  # Checked before anything is written: an archive from another work area
  # would put its files where no case reads them.
  from=$(zstd -dc "$E2E_BASE_IN" | tar -P -xOf - "${E2E_WORK}/base-from" 2> /dev/null || true)
  if [ "$from" != "${E2E_WORK} ${REPO}" ]; then
    echo "ERROR: ${E2E_BASE_IN} was not packed with E2E_WORK=${E2E_WORK} and the checkout at ${REPO} (it says: ${from:-nothing})" >&2
    exit 1
  fi
  zstd -dc "$E2E_BASE_IN" | sudo tar -P -xpf -
  echo "Unpacked $(zstd -dc "$E2E_BASE_IN" | tar -tf - | wc -l) entries"
  # The archive is about 4 GB on the disk the cases write to. Shard 1 of run
  # 41524201317 filled that disk while the archive was still there.
  rm -f "$E2E_BASE_IN"
fi

# --- 3. Images ----------------------------------------------------------------
# Only a full run pre-pulls. A case reaches most images through the steps it
# runs, so the images of a few cases cannot be told apart from the rest; a
# partial run lets its cases pull their own, without the retry.
if [ "$PATTERN" = "*" ]; then
  echo "=== Pulling the images the cases use ==="
  pull_all
else
  echo "=== Partial run ('${PATTERN}'): no pre-pull, each case pulls what it uses ==="
fi

# --- Trace -------------------------------------------------------------------
# With E2E_TRACE=true each case's reads are told apart by atime, which only
# moves on every read under strictatime, and its writes by ctime, which no
# copy can carry over from the source. trace_scan prints "CTIME ATIME PATH" for
# every file under the work area and the checkout, without the fixture, the
# report and .git. The logs stay in: a case may read an earlier case's log.
TRACE_DIR="${REPORT_DIR}/trace-state"
trace_scan() {
  sudo find "$E2E_WORK" "$REPO" \( -path "$FIXTURE_DIR" -o -path "$REPORT_DIR" -o -path "${REPO}/.git" \) \
    -prune -o -type f -printf '%C@\t%A@\t%p\n' 2> /dev/null || true
}
if [ "$TRACE" = true ]; then
  mkdir -p "$TRACE_DIR"
  for d in "$E2E_WORK" "$REPO"; do sudo mount -o remount,strictatime "$(findmnt -n -o TARGET -T "$d")"; done
  # Positive control: a second read must move atime too (relatime moves it on
  # the first read after a write only), or the trace would miss readers.
  for d in "$E2E_WORK" "$REPO"; do
    probe="${d}/.atime-probe"
    echo x > "$probe"; sleep 1.1; cat "$probe" > /dev/null
    before=$(stat -c %x "$probe"); sleep 1.1; cat "$probe" > /dev/null
    if [ "$(stat -c %x "$probe")" = "$before" ]; then
      echo "ERROR: E2E_TRACE: a second read did not move atime under ${d}" >&2; exit 1
    fi
    rm -f "$probe"
  done
  # Every file there before the first case was written by the layout or the base job.
  trace_scan | awk -F'\t' '{print $3 "\tsetup"}' > "${TRACE_DIR}/writers.tsv"
  : > "${REPORT_DIR}/trace.tsv"
fi
# trace_case NAME START SECONDS: one trace.tsv line for the case that just ran,
# "NAME SECONDS WRITERS READS", WRITERS being the comma-separated cases (or
# "setup") that wrote the files it read, READS how many such files it read.
trace_case() {
  trace_scan > "${TRACE_DIR}/scan.tsv"
  awk -F'\t' -v s="$2" -v k="$1" -v secs="$3" -v out="${TRACE_DIR}/writers.new" '
    FILENAME == ARGV[1] { w[$1] = $2; next }
    {
      p = $3
      if (p in w) {
        if ($2 + 0 >= s + 0) { reads++; if (!(w[p] in seen)) { seen[w[p]] = 1; list = (list == "" ? w[p] : list "," w[p]) } }
        if ($1 + 0 >= s + 0) w[p] = k
      } else w[p] = k
      print p "\t" w[p] > out
    }
    END { printf "%s\t%s\t%s\t%d\n", k, secs, (list == "" ? "-" : list), reads }' \
    "${TRACE_DIR}/writers.tsv" "${TRACE_DIR}/scan.tsv" >> "${REPORT_DIR}/trace.tsv"
  mv "${TRACE_DIR}/writers.new" "${TRACE_DIR}/writers.tsv"
}

# --- 4. Cases -----------------------------------------------------------------
mapfile -t CASES < <(cd "${REPO}/tests/e2e" && find . -maxdepth 1 -type f -name '*.sh' ! -name lib.sh -printf '%f\n' | LC_ALL=C sort)
if [ -n "$CASE_LIST" ]; then
  unknown=$(grep -vxF -f <(printf '%s\n' "${CASES[@]%.sh}") "$CASE_LIST" || true)
  [ -z "$unknown" ] || { echo "ERROR: E2E_CASES lists names that are no case: ${unknown//$'\n'/ }" >&2; exit 1; }
fi
names=() results=() times=()
failed=0
for c in "${CASES[@]}"; do
  name="${c%.sh}"
  if [ -n "$CASE_LIST" ]; then
    grep -qxF "$name" "$CASE_LIST" || continue
  else
    # shellcheck disable=SC2053  # PATTERN is a glob on purpose
    [[ "$c" == $PATTERN ]] || continue
  fi
  log="${LOG_DIR}/${name}.log"
  start=$(date +%s)
  tstart=$(date +%s.%N)
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
  secs=$(( $(date +%s) - start ))
  names+=("$name"); results+=("$result"); times+=("${secs}s")
  printf '%s\t%s\t%s\t%s\n' "$name" "$result" "$secs" "$JOB" >> "${REPORT_DIR}/results.tsv"
  [ "$TRACE" = true ] && trace_case "$name" "$tstart" "$secs"
  echo "${name}: ${result}"
done
if [ "${#names[@]}" -eq 0 ]; then
  echo "ERROR: no case matches '${PATTERN}'${CASE_LIST:+ or is listed in ${CASE_LIST}}" >&2
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
  echo "Logs are in the e2e-logs artifacts of this run."
  {
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
  } | tee "${REPORT_DIR}/extra.md"
} | tee -a "$SUMMARY"

# --- 6. Pack for the shards ---------------------------------------------------
# Everything the cases wrote or changed, at its absolute path, owners, modes
# and times (ns) kept, the base cases' logs too (a case may read one).
# Directories go in as entries alone (--no-recursion), so the reference and
# the fixture copies of the layout stay out.
if [ -n "${E2E_BASE_OUT:-}" ] && [ "$failed" -eq 0 ]; then
  echo "${E2E_WORK} ${REPO}" > "${E2E_WORK}/base-from"
  # Beside the archive too, so a shard learns the work area before it unpacks.
  cp "${E2E_WORK}/base-from" "${E2E_BASE_OUT}.from"
  sudo find "$E2E_WORK" "$REPO" \( -path "$FIXTURE_DIR" -o -path "$REPORT_DIR" \
      -o -path "${REPO}/.git" -o -path "$E2E_NOTES" -o -path "${E2E_WORK}/all-cases-ran" \) -prune \
      -o -cnewer "${E2E_WORK}/.layout-done" -print0 \
    | sudo tar -P --format=posix --no-recursion --null -T - -cf - | zstd -T0 -3 -q -o "$E2E_BASE_OUT"
  echo "Packed the base outputs: $(du -h "$E2E_BASE_OUT" | cut -f1), $(zstd -dc "$E2E_BASE_OUT" | tar -tf - | wc -l) entries. The largest:"
  zstd -dc "$E2E_BASE_OUT" | tar -tvf - | sort -k3,3n | tail -n 5
fi

[ "$failed" -eq 0 ]
