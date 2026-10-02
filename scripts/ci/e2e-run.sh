#!/usr/bin/env bash
# e2e-run.sh — run every end-to-end case in tests/e2e/ on the HG002 fixture.
#
# Usage: scripts/ci/e2e-run.sh [pattern]
#   pattern  optional shell glob; runs only the cases whose file name matches
#            (e.g. '2*' or '*cyrius*'). Later cases read earlier cases' outputs,
#            so a partial run is for debugging only.
#
# Env: E2E_WORK  work area (default ${RUNNER_TEMP:-/tmp}/e2e-work), about 25 GB
#      GH_TOKEN  token for `gh release download` (the job's GITHUB_TOKEN)
#
# What it does:
#   1. downloads the release named in tests/fixtures/VERSION and checks
#      SHA256SUMS; when a pull request bumps VERSION, waits up to 75 minutes
#      for the build-fixture job to publish it;
#   2. lays out GENOME_DIR the way setup.sh and step 13 would leave it;
#   3. runs tests/e2e/*.sh in name order (C locale: numbered cases first, then
#      the <package-key>-*.sh cases later packages add) and keeps going after a
#      failure, so one run lists every broken step;
#   4. prints a case / result / time / log table to the job summary and exits
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
# Today's scripts read reference/Homo_sapiens_assembly38.fasta; the no-alt name
# is there too so a script that switches to it finds the same file.
echo "=== Layout ${GENOME_DIR} ==="
mkdir -p "${GENOME_DIR}/reference" "${GENOME_DIR}/clinvar" "${GENOME_DIR}/annotations" \
  "${GENOME_DIR}/${SAMPLE}/fastq" "${GENOME_DIR}/${SAMPLE}/vep"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38"
if [ ! -s "${REF}.fasta" ]; then
  gzip -dc "${FIXTURE_DIR}/fixture_ref.fa.gz" > "${REF}.fasta.tmp"
  mv "${REF}.fasta.tmp" "${REF}.fasta"
fi
cp "${FIXTURE_DIR}/fixture_ref.fa.gz.fai" "${REF}.fasta.fai"
cp "${FIXTURE_DIR}/fixture_ref.dict" "${REF}.dict"
for ext in fasta fasta.fai dict; do
  ln -f "${REF}.${ext}" "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.${ext}"
done
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

# --- 3. Cases -----------------------------------------------------------------
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

# --- 4. Report ----------------------------------------------------------------
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
