#!/usr/bin/env bash
# The interval rows of bin/cn_loci.tsv are what scripts/ci/build-cn-intervals.py
# writes from its pinned downloads (ClinGen's recurrent CNV BED V2.1 and three
# GIAB v3.6 stratifications, each checked against its sha256): a hand edit of
# an interval row fails here. Control: a copy of the catalog with one target
# moved by 1 kb must fail the same check. The downloads come from the e2e
# data cache when the runner filled it (entry catalogs, mirror first;
# scripts/ci/e2e-run.sh), else from their sources.
. "$(dirname "$0")/lib.sh"

CACHE="${CASE_TMP}/downloads"
[ -z "${E2E_DATA_CACHE:-}" ] || CACHE="${E2E_DATA_CACHE}/catalogs/cn-intervals"
mkdir -p "$CACHE"
echo "+ scripts/ci/build-cn-intervals.py --check"
python3 "${REPO}/scripts/ci/build-cn-intervals.py" --cache "$CACHE" --check 2>&1 | tee "$STEP_LOG"
check_eq "build-cn-intervals.py --check exits 0 on the committed catalog" "${PIPESTATUS[0]}" 0
check "the log says all 90 interval rows were checked" grep -q "the 90 interval rows are what the pinned downloads give" "$STEP_LOG"

# Control: one interval target moved by 1 kb
EDITED="${CASE_TMP}/edited.tsv"
awk -F'\t' -v OFS='\t' '$1 == "22q11.2_DGS_proximal_A-B" && !done {
    split($4, a, ";"); split(a[1], b, ":"); split(b[2], c, "-")
    a[1] = b[1] ":" c[1] + 1000 "-" c[2] + 1000 (b[3] != "" ? ":" b[3] : "")
    t = a[1]; for (i = 2; i in a; i++) t = t ";" a[i]; $4 = t; done = 1 } 1' \
  "${REPO}/bin/cn_loci.tsv" > "$EDITED"
# shellcheck disable=SC2016  # the positional arguments of bash -c
check "control: the edited copy differs from the catalog" bash -c '! cmp -s "$1" "$2"' _ "$EDITED" "${REPO}/bin/cn_loci.tsv"
python3 "${REPO}/scripts/ci/build-cn-intervals.py" --cache "$CACHE" --catalog "$EDITED" --check > "${CASE_TMP}/edited.log" 2>&1
check_eq "control: --check fails (exit 1) on a catalog with one interval target moved" "$?" 1
cat "${CASE_TMP}/edited.log"

finish
