# shellcheck shell=bash
# lib.sh — helpers shared by the e2e cases in tests/e2e/. Every case starts with
#   . "$(dirname "$0")/lib.sh"
# and ends with `finish`. scripts/ci/e2e-run.sh sets REPO, GENOME_DIR, SAMPLE,
# FIXTURE_DIR, E2E_WORK and E2E_NOTES before it runs a case.
#
# A case never stops at the first failed check: it runs them all, prints one
# [PASS] or [FAIL] line each, and `finish` exits 1 if any failed. Every check is
# a count, a column or a content match; an exit code alone is never enough.

set -uo pipefail
: "${REPO:?run the cases through scripts/ci/e2e-run.sh}"
: "${GENOME_DIR:?}" "${SAMPLE:?}" "${FIXTURE_DIR:?}" "${E2E_WORK:?}" "${E2E_NOTES:?}"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"

CASE_NAME="$(basename "$0" .sh)"
CASE_TMP="${E2E_WORK}/tmp/${CASE_NAME}"
mkdir -p "$CASE_TMP"
FAILS=0
STEP_RC=""
STEP_LOG="${CASE_TMP}/step.log"

pass() { echo "[PASS] $*"; }
fail() { echo "[FAIL] $*"; FAILS=$((FAILS + 1)); }

# check <description> <command...>: passes when the command succeeds.
check() {
  local desc=$1; shift
  if "$@"; then pass "$desc"; else fail "$desc"; fi
}
check_eq() {
  if [ "$2" = "$3" ]; then pass "$1 (${2})"; else fail "$1 (got '${2}', want '${3}')"; fi
}
check_ge() {
  if [[ "$2" =~ ^[0-9]+$ ]] && [ "$2" -ge "$3" ]; then pass "$1 (${2} >= ${3})"
  else fail "$1 (got '${2}', want >= ${3})"; fi
}
has()   { grep -Eq -- "$1" <<< "$2"; }
lacks() { ! grep -Eq -- "$1" <<< "$2"; }

# run_step <script> [args...]: runs scripts/<script>, output to the case log and
# to STEP_LOG; the exit code lands in STEP_RC.
run_step() {
  echo "+ scripts/$*"
  "${REPO}/scripts/$1" "${@:2}" 2>&1 | tee "$STEP_LOG"
  STEP_RC=${PIPESTATUS[0]}
  echo "+ exit ${STEP_RC}"
}
check_step_exit() { check_eq "scripts/$1 exits 0" "$STEP_RC" 0; }

# Tools from the pinned images; paths are relative to GENOME_DIR.
in_genome() {
  local image=$1; shift
  docker run --rm -i -v "${GENOME_DIR}:/genome" -w /genome "$image" "$@"
}
bcf() { in_genome "$BCFTOOLS_IMAGE" bcftools "$@"; }
sam() { in_genome "$SAMTOOLS_IMAGE" samtools "$@"; }
vcf_ok()    { bcf view -h "$1" >/dev/null 2>&1; }
vcf_count() { bcf view -H "$@" 2>/dev/null | wc -l | tr -d ' '; }
nonempty()  { [ -s "${GENOME_DIR}/$1" ]; }

# planted <column>: a field of the planted ClinVar record (chrom pos ref alt gene clinvar_id).
planted() {
  awk -F'\t' -v k="$1" 'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i} NR == 2 {print $col[k]}' \
    "${FIXTURE_DIR}/planted.tsv"
}

# sample_side_hits <dir> <gene|"">: lines at the planted position in the files
# under <dir> that describe the sample, optionally only those naming <gene>.
# isec/0001.vcf and isec/0003.vcf are ClinVar's own records, and *_pass.vcf.gz
# is the sample's whole PASS call set, so none of them counts as a hit.
sample_side_hits() {
  local dir=$1 gene=$2
  [ -d "$dir" ] || { echo 0; return; }
  find "$dir" -type f \( -name '*.vcf' -o -name '*.vcf.gz' -o -name '*.tsv' -o -name '*.txt' \) \
      ! -path '*/isec/0001.vcf' ! -path '*/isec/0003.vcf' ! -name '*_pass.vcf.gz' -print0 \
    | xargs -0 -r -n1 gzip -cdf \
    | awk -F'\t' -v c="$(planted chrom)" -v p="$(planted pos)" '$1 == c && $2 == p' \
    | grep -c -- "$gene" || true
}

# pharmcat_called <report.json>: genes with a named diplotype (flat or nested
# `genes` map, sourceDiplotypes or recommendationDiplotypes), one per line.
pharmcat_called() {
  python3 - "$1" <<'PY'
import json, sys
try:
    data = json.load(open(sys.argv[1]))
except Exception as e:
    print(f"unreadable report: {e}", file=sys.stderr)
    sys.exit(0)
def called(g):
    for d in (g.get("sourceDiplotypes") or g.get("recommendationDiplotypes") or []):
        names = [(d.get(a) or {}).get("name", "") for a in ("allele1", "allele2")]
        if any(n and n.lower() not in ("unknown", "none", "?") for n in names):
            return True
    return False
genes = data.get("genes") or {}
out = set()
for key, val in genes.items():
    if not isinstance(val, dict):
        continue
    if "sourceDiplotypes" in val or "recommendationDiplotypes" in val:
        if called(val):
            out.add(key)
    else:
        out.update(name for name, g in val.items() if isinstance(g, dict) and called(g))
print("\n".join(sorted(out)))
PY
}

finish() {
  if [ "$FAILS" -gt 0 ]; then
    echo "${CASE_NAME}: ${FAILS} check(s) failed"
    exit 1
  fi
  echo "${CASE_NAME}: all checks passed"
  exit 0
}
