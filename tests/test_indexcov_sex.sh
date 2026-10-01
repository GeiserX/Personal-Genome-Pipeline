#!/usr/bin/env bash
# scripts/16-indexcov.sh: inferred sex comes from the .ped 'sex' column, and a
# declared sex that disagrees stops the step.
#
# goleft writes "#family_id sample_id paternal_id maternal_id sex phenotype
# CNchrX CNchrY ...". The old parser tested column 6 (phenotype, always -9), so
# every sample was printed as female, and the declared sex was never compared.
#
# A fake `docker` on PATH stands in for goleft and writes a synthetic .ped, so
# the real script runs end to end without containers.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

FAILS=0
fail() { echo "FAIL: $*"; FAILS=$((FAILS + 1)); }
pass() { echo "ok:   $*"; }

# Fake docker: finds --directory /genome/<...> and writes $FAKE_PED there.
mkdir -p "${WORK}/bin"
cat > "${WORK}/bin/docker" <<'FAKE'
#!/usr/bin/env bash
dir=""
while [ $# -gt 0 ]; do
  if [ "$1" = "--directory" ]; then dir=$2; fi
  shift
done
[ -n "$dir" ] || { echo "fake docker: no --directory" >&2; exit 1; }
host="${GENOME_DIR}${dir#/genome}"
mkdir -p "$host"
printf '%b' "$FAKE_PED" > "${host}/$(basename "$host")-indexcov.ped"
FAKE
chmod +x "${WORK}/bin/docker"

GD="${WORK}/genome"
mkdir -p "${GD}/S1/aligned"
: > "${GD}/S1/aligned/S1_sorted.bam"
: > "${GD}/S1/aligned/S1_sorted.bam.bai"

HEADER='#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\tbins.out\tbins.lo\tbins.hi\tbins.in\tslope\tp.out\tPC1\tPC2\tPC3\tPC4\tPC5\n'
MALE_ROW='S1\tS1\t-9\t-9\t1\t-9\t1.01\t0.97\t40\t12\t9\t0.99\t0.98\t0.02\t0\t0\t0\t0\t0\n'
FEMALE_ROW='S1\tS1\t-9\t-9\t2\t-9\t1.98\t0.00\t40\t12\t9\t0.99\t0.98\t0.02\t0\t0\t0\t0\t0\n'

# run_case <ped content> [declared sex]; sets OUT and RC
run_case() {
  local ped=$1; shift
  set +e
  OUT=$(PATH="${WORK}/bin:${PATH}" GENOME_DIR="$GD" FAKE_PED="$ped" bash "${REPO}/scripts/16-indexcov.sh" S1 "$@" 2>&1)
  RC=$?
  set -e
}

run_case "${HEADER}${MALE_ROW}"
if [ "$RC" -eq 0 ] && grep -q 'Predicted sex: male' <<<"$OUT"; then
  pass ".ped sex=1 prints male"
else
  fail ".ped sex=1 did not print male (rc=${RC}): $(grep -i 'sex' <<<"$OUT" | tr '\n' '|')"
fi

run_case "${HEADER}${FEMALE_ROW}"
if [ "$RC" -eq 0 ] && grep -q 'Predicted sex: female' <<<"$OUT"; then
  pass ".ped sex=2 prints female"
else
  fail ".ped sex=2 did not print female (rc=${RC}): $(grep -i 'sex' <<<"$OUT" | tr '\n' '|')"
fi

run_case "${HEADER}${MALE_ROW}" female
if [ "$RC" -ne 0 ] && grep -q 'declared female, inferred male' <<<"$OUT"; then
  pass "declared female vs inferred male exits non-zero naming both"
else
  fail "declared female vs inferred male: rc=${RC}, output: $(tail -3 <<<"$OUT" | tr '\n' '|')"
fi

run_case "${HEADER}${MALE_ROW}" male
if [ "$RC" -eq 0 ] && grep -q 'Sex check: OK' <<<"$OUT"; then
  pass "declared male vs inferred male passes"
else
  fail "declared male vs inferred male: rc=${RC}"
fi

set +e
OUT=$(PATH="${WORK}/bin:${PATH}" GENOME_DIR="$GD" FAKE_PED="${HEADER}${MALE_ROW}" SEX_CHECK=warn \
  bash "${REPO}/scripts/16-indexcov.sh" S1 female 2>&1)
RC=$?
set -e
if [ "$RC" -eq 0 ] && grep -q 'SEX CHECK MISMATCH' <<<"$OUT"; then
  pass "SEX_CHECK=warn prints the mismatch and exits 0"
else
  fail "SEX_CHECK=warn: rc=${RC}"
fi

echo ""
if [ "$FAILS" -gt 0 ]; then
  echo "${FAILS} check(s) failed"
  exit 1
fi
echo "All checks passed"
