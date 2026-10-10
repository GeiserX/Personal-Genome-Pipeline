#!/usr/bin/env bash
# Delly 2.7.0 genotypes chrX and chrY by sex (`delly sr --sex`). Unless told,
# it infers the sex from coverage (`--sex auto`). The pipeline already knows
# the declared sex, so it hands it over:
#   scripts/19-delly.sh S female  -> delly sr --sex female
#   scripts/19-delly.sh S male    -> delly sr --sex male
#   scripts/19-delly.sh S         -> delly sr --sex auto (Delly infers)
#   scripts/19-delly.sh S other   -> exits non-zero before Delly runs
# and the Nextflow DELLY process passes meta.sex the same way, auto when the
# samplesheet gives none. `--sex none` (no sex-aware genotyping) is never used.
#
# A fake `docker` on PATH logs every call, so the real script runs end to end
# without containers.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/docker" <<'FAKE'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$FAKE_DOCKER_LOG"
FAKE
chmod +x "${WORK}/bin/docker"

GD="${WORK}/genome"
mkdir -p "${GD}/S1/aligned" "${GD}/reference"
: > "${GD}/S1/aligned/S1_sorted.bam"
: > "${GD}/S1/aligned/S1_sorted.bam.bai"
: > "${GD}/reference/GRCh38_no_alt_analysis_set.fasta"
: > "${GD}/reference/GRCh38_no_alt_analysis_set.fasta.fai"
LOG="${WORK}/docker.log"

# run_case [declared sex]; sets OUT, RC and DELLY (the delly sr call, or empty)
run_case() {
  : > "$LOG"
  set +e
  OUT=$(PATH="${WORK}/bin:${PATH}" GENOME_DIR="$GD" FAKE_DOCKER_LOG="$LOG" bash "${REPO}/scripts/19-delly.sh" S1 "$@" 2>&1)
  RC=$?
  set -e
  DELLY=$(grep ' delly sr ' "$LOG" || true)
}

# expect_sex <declared sex or ""> <value Delly must get>
expect_sex() {
  local declared=$1 want=$2 label=${1:-none declared}
  if [ -n "$declared" ]; then run_case "$declared"; else run_case; fi
  if [ "$RC" -ne 0 ]; then
    fail "${label}: the step exited ${RC}: $(tail -3 <<<"$OUT" | tr '\n' '|')"
  elif [ -z "$DELLY" ]; then
    fail "${label}: delly sr was never run"
  elif ! grep -q -- " --sex ${want} " <<<"$DELLY "; then
    fail "${label}: delly sr did not get --sex ${want}: ${DELLY}"
  elif [ "$(grep -o -- ' --sex ' <<<"$DELLY" | wc -l)" -ne 1 ]; then
    fail "${label}: delly sr got --sex more than once: ${DELLY}"
  else
    pass "${label}: delly sr --sex ${want}"
  fi
}

expect_sex female female
expect_sex male male
expect_sex "" auto

run_case unknown
if [ "$RC" -ne 0 ] && [ -z "$DELLY" ] && grep -q "sex must be 'male' or 'female'" <<<"$OUT"; then
  pass "an unknown sex stops the step before Delly runs"
else
  fail "an unknown sex: rc=${RC}, delly call: '${DELLY}', output: $(tail -3 <<<"$OUT" | tr '\n' '|')"
fi

# The Nextflow process: its delly sr command carries --sex from meta.sex,
# auto when the samplesheet gives no sex.
MOD="${REPO}/modules/local/delly/main.nf"
block=$(awk '/^process DELLY \{/ { on = 1 } on && /^}/ { exit } on' "$MOD")
if grep -q "def sex_arg = \"--sex \${meta.sex ?: 'auto'}\"" <<<"$block" \
   && awk '/delly sr/ { on = 1 } on && /\$\{sex_arg\}/ { found = 1 } on && /^ *"""/ { on = 0 } END { exit !found }' <<<"$block"; then
  pass "DELLY passes --sex \${meta.sex ?: 'auto'} to delly sr"
else
  fail "DELLY in ${MOD#"${REPO}"/} does not pass --sex from meta.sex (auto when unset) to delly sr"
fi
if grep -q -- '--sex none' "$MOD" "${REPO}/scripts/19-delly.sh"; then
  fail "--sex none turns off Delly's sex-aware genotyping and must not be used"
fi

if [ "$FAILS" -ne 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "all Delly sex checks passed"
