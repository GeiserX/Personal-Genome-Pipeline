#!/usr/bin/env bash
# setup.sh --masked-reference: the pins in versions.env and the masking
# program, without Docker or the network.
#   1. GRC_EXCLUSIONS_SHA256 is the sha256 of the GRC's list as GIAB
#      publishes it (the 12 rows below, read from the file on 2026-11-04);
#      the same check on a copy of versions.env with one changed digit fails.
#      GRC_MASKED_REF_SHA256 is a sha256.
#   2. MASK_FASTA_PY, run with python3 exactly as setup.sh holds it, on a
#      small FASTA: every base of a row is N, every other byte is the same
#      (line breaks, lower case, the header and a sequence with no row), and
#      it names the row on a sequence the FASTA lacks. The same program with
#      the BED shifted by one base gives another file, so the comparison can
#      fail. The e2e case fixitfelix-chr21-1 runs the whole option on chr21.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
FAILS=0
pass() { echo "[PASS] $*"; }
fail() { echo "[FAIL] $*"; FAILS=$((FAILS + 1)); }
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }

# --- 1. the pins ---------------------------------------------------------------------
printf '%s\t%s\t%s\n' \
  chr9 70701806 70719798  chr13 86252979 86269524  chr13 111753855 111789114 \
  chr16 34339329 34456556  chr21 5010000 5166246  chr21 5966593 6161371 \
  chr21 6427258 6580181  chr21 6789085 6934219  chr21 7743700 7865746 \
  chrX 37085011 37099262  chr11_KI270721v1_random 1 100316  chr22_KI270734v1_random 1 165050 \
  > "${WORK}/grc.bed"
GRC_SHA=$(sha "${WORK}/grc.bed")

# pin_ok ENV: GRC_EXCLUSIONS_SHA256 in ENV is the sha256 of the 12 rows.
pin_ok() { [ "$(sed -n 's/^GRC_EXCLUSIONS_SHA256="\([0-9a-f]*\)"$/\1/p' "$1")" = "$GRC_SHA" ]; }
if pin_ok "${REPO}/versions.env"; then pass "GRC_EXCLUSIONS_SHA256 is the sha256 of the GRC list (${GRC_SHA})"
else fail "GRC_EXCLUSIONS_SHA256 in versions.env is not ${GRC_SHA}"; fi
sed 's/^\(GRC_EXCLUSIONS_SHA256="\)./\1x/' "${REPO}/versions.env" > "${WORK}/changed.env"
if ! pin_ok "${WORK}/changed.env" && ! cmp -s "${WORK}/changed.env" "${REPO}/versions.env"; then
  pass "control: a changed GRC_EXCLUSIONS_SHA256 fails the same check"
else fail "control: a changed GRC_EXCLUSIONS_SHA256 passed"; fi
if grep -Eq '^GRC_MASKED_REF_SHA256="[0-9a-f]{64}"$' "${REPO}/versions.env"; then pass "GRC_MASKED_REF_SHA256 is a sha256"
else fail "GRC_MASKED_REF_SHA256 is not a sha256"; fi
if grep -q '^GRC_EXCLUSIONS_URL="https://.*/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed"$' "${REPO}/versions.env"; then
  pass "GRC_EXCLUSIONS_URL names GIAB's copy of the list"
else fail "GRC_EXCLUSIONS_URL is not GIAB's T2Tv2 list"; fi

# --- 2. the masking program ----------------------------------------------------------
awk '/^MASK_FASTA_PY=\$\(cat <<.PY.$/ {on = 1; next} on && /^PY$/ {exit} on' "${REPO}/scripts/setup.sh" > "${WORK}/mask.py"
if grep -q 'def\|rows = {}' "${WORK}/mask.py"; then pass "MASK_FASTA_PY read from setup.sh ($(wc -l < "${WORK}/mask.py" | tr -d ' ') lines)"
else fail "MASK_FASTA_PY not found in setup.sh"; fi

printf '>s1 first\nACGTacg\ntACGTAC\nGT\n>s2\nAAAAAAA\nCC\n' > "${WORK}/in.fa"
printf 'browser position s1:1-3\n#seq\tstart\tend\ns1\t5\t9\ns1\t15\t16\nmissing\t0\t10\n' > "${WORK}/rows.bed"
# s1 is ACGTacgtACGTACGT: 0-based [5, 9) and [15, 16) are N.
printf '>s1 first\nACGTaNN\nNNCGTAC\nGN\n>s2\nAAAAAAA\nCC\n' > "${WORK}/want.fa"
OUT=$(python3 "${WORK}/mask.py" "${WORK}/in.fa" "${WORK}/rows.bed" "${WORK}/out.fa" 2>&1)
echo "  ${OUT}"
if cmp -s "${WORK}/out.fa" "${WORK}/want.fa"; then pass "the rows are N and every other byte is the same"
else fail "masked FASTA differs from the expected one"; diff "${WORK}/want.fa" "${WORK}/out.fa"; fi
case "$OUT" in
  "masked 2 of 3 rows (5 bases); rows on sequences the FASTA lacks: missing") pass "it reports 2 of 3 rows and names the missing sequence" ;;
  *) fail "report: '${OUT}'" ;;
esac
printf 's1\t6\t10\ns1\t15\t16\n' > "${WORK}/shifted.bed"
python3 "${WORK}/mask.py" "${WORK}/in.fa" "${WORK}/shifted.bed" "${WORK}/shifted.fa" > /dev/null 2>&1
if [ -s "${WORK}/shifted.fa" ] && ! cmp -s "${WORK}/shifted.fa" "${WORK}/want.fa"; then pass "control: a BED one base off gives another file"
else fail "control: a BED one base off gave the same file"; fi

if [ "$FAILS" -gt 0 ]; then echo "test_fixitfelix: ${FAILS} check(s) failed"; exit 1; fi
echo "test_fixitfelix: all checks passed"
