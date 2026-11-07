#!/usr/bin/env bash
# assets/callability_genes.tsv is what scripts/ci/build-callability.sh writes
# from its pinned downloads, byte for byte: a hand edit, or a change to
# assets/paralog_loci.tsv, assets/pgx_genes.tsv or bin/cn_loci.tsv that moves
# an owner without a rebuild, fails here. Control: a copy of the committed
# file with one owner cell changed must differ from the rebuild. The
# downloads come from the e2e data cache when the runner filled it (entry
# catalogs, mirror first; scripts/ci/e2e-run.sh), else from their sources.
. "$(dirname "$0")/lib.sh"

OUT="${CASE_TMP}/callability_genes.tsv"
COMMITTED="${REPO}/assets/callability_genes.tsv"
rm -f "$OUT"
echo "+ scripts/ci/build-callability.sh"
CACHE="${CASE_TMP}/downloads"
[ -z "${E2E_DATA_CACHE:-}" ] || CACHE="${E2E_DATA_CACHE}/catalogs/callability"
CALLABILITY_CACHE="$CACHE" bash "${REPO}/scripts/ci/build-callability.sh" "$OUT" 2>&1 | tee "$STEP_LOG"
check_eq "build-callability.sh exits 0" "${PIPESTATUS[0]}" 0
echo "${CASE_NAME}: $(grep -c '^Downloading ' "$STEP_LOG") of 7 downloads not found in ${CACHE}" >> "$E2E_NOTES"

# same A B: the two files hold the same bytes
same() { cmp -s "$1" "$2"; }
awk -F'\t' -v OFS='\t' '!/^#/ && $1 == "HLA-DRB3" {$7 = "pypgx"} 1' "$COMMITTED" > "${CASE_TMP}/edited.tsv"
if same "$OUT" "${CASE_TMP}/edited.tsv"; then
  fail "control: the rebuild matches a hand-edited catalog, so this check cannot fail"
else
  pass "control: the rebuild differs from a hand-edited catalog"
fi
if same "$OUT" "$COMMITTED"; then
  pass "the rebuilt catalog is the committed assets/callability_genes.tsv, byte for byte"
else
  fail "the rebuilt catalog differs from assets/callability_genes.tsv; rebuild it with scripts/ci/build-callability.sh"
  diff "$COMMITTED" "$OUT" | head -20
fi
check_ge "the rebuild lists the 273 CMRG genes" "$(awk -F'\t' '!/^#/ && $6 ~ /CMRG/' "$OUT" 2>/dev/null | wc -l | tr -d ' ')" 273

finish
