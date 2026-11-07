# shellcheck shell=bash
# lib.sh — helpers shared by the e2e cases in tests/e2e/. Every case starts with
#   . "$(dirname "$0")/lib.sh"
# and ends with `finish`. scripts/ci/e2e-run.sh sets REPO, GENOME_DIR, SAMPLE,
# FIXTURE_DIR, E2E_WORK and E2E_NOTES before it runs a case.
#
# A case never stops at the first failed check: it runs them all, prints one
# [PASS] or [FAIL] line each, and `finish` exits 1 if any failed. Every check is
# a count, a column or a content match; an exit code alone is never enough.
# A case that ran no check at all fails too, unless it set NO_CHECKS_REASON
# before `finish` (a case that has nothing to do on this run says why).

set -uo pipefail
: "${REPO:?run the cases through scripts/ci/e2e-run.sh}"
: "${GENOME_DIR:?}" "${SAMPLE:?}" "${FIXTURE_DIR:?}" "${E2E_WORK:?}" "${E2E_NOTES:?}"
# shellcheck source=../../versions.env
. "${REPO}/versions.env"

CASE_NAME="$(basename "$0" .sh)"
CASE_TMP="${E2E_WORK}/tmp/${CASE_NAME}"
mkdir -p "$CASE_TMP"
FAILS=0
CHECKS=0
NO_CHECKS_REASON=""
STEP_RC=""
STEP_LOG="${CASE_TMP}/step.log"

pass() { echo "[PASS] $*"; CHECKS=$((CHECKS + 1)); }
fail() { echo "[FAIL] $*"; CHECKS=$((CHECKS + 1)); FAILS=$((FAILS + 1)); }

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

# --- The ancestry panel, from the Actions cache or the mirror when it is there ---
# EBI's FTP server fails now and then (curl exit 7). The e2e workflow keeps the
# PGS Catalog's HAPNEST panel in the Actions cache (E2E_DATA_CACHE/hapnest,
# keyed on the pin below) for prs-2-ancestry and
# prs-scores-and-run-options-2-panel-chrx. The pin is the sha256 of the file
# EBI served on 2026-11-04 (281,185,251 bytes, md5
# 3630f2ec3343599f3c55827ddbe38f60 as EBI's md5s.txt lists it).
HAPNEST_SHA256=8ba05cffc4b816f8f8aafd78add52b90ab80315e930936277e09861c5d134aa3
# NET_BLOCKED: an env prefix that sends every HTTP(S) request to a closed port.
NET_BLOCKED=(env http_proxy=http://127.0.0.1:9 https_proxy=http://127.0.0.1:9 HTTP_PROXY=http://127.0.0.1:9
             HTTPS_PROXY=http://127.0.0.1:9 no_proxy= NO_PROXY= FETCH_TRIES=1 FETCH_WAIT=0)
_sha256() { sha256sum "$1" 2>/dev/null | cut -d' ' -f1; }
# setup_ancestry_panel NAME PANEL: setup.sh --ancestry-panel for panel NAME
# (installed at PANEL), its log in STEP_LOG, with checks. The HAPNEST panel
# comes from the Actions cache, an earlier install at PANEL, or else from
# E2E_DATA_MIRROR (the e2e-data-1 release, which holds a copy), each checked
# against HAPNEST_SHA256; then setup.sh runs with the network blocked. When
# none has it, setup.sh downloads it from EBI as a user's run does. With E2E_HAPNEST_KEEP=true a
# panel that matches the pin and did not come from the cache goes into it.
setup_ancestry_panel() {
  local name=$1 panel=$2 cached="${E2E_DATA_CACHE:-}/hapnest/${1}.tar.zst" from=""
  if [ "$name" = GRCh38_HAPNEST_reference ] && [ -n "${E2E_DATA_CACHE:-}" ] && [ -s "$cached" ]; then
    if [ "$(_sha256 "$cached")" = "$HAPNEST_SHA256" ]; then
      mkdir -p "${panel%/*}"
      ln -f "$cached" "$panel" 2>/dev/null || cp -f "$cached" "$panel"
      from="the Actions cache"
    else
      echo "Data cache: ${cached} does not match its pin; trying the mirror, then setup.sh downloads it"
      rm -f "$cached"
    fi
  fi
  # A panel an earlier case installed, when it matches the pin, is used as it is.
  if [ -z "$from" ] && [ "$name" = GRCh38_HAPNEST_reference ] && [ -s "$panel" ] \
     && [ "$(_sha256 "$panel")" = "$HAPNEST_SHA256" ]; then
    from="an earlier install"
  fi
  if [ -z "$from" ] && [ "$name" = GRCh38_HAPNEST_reference ] && [ -n "${E2E_DATA_MIRROR:-}" ]; then
    mkdir -p "${panel%/*}"
    # One try, cut when it stalls or takes over 15 minutes: EBI is next.
    if curl -fsSL --connect-timeout 30 --speed-limit 1024 --speed-time 120 --max-time 900 \
         -o "${panel}.mirror" "${E2E_DATA_MIRROR%/}/${name}.tar.zst" \
       && [ "$(_sha256 "${panel}.mirror")" = "$HAPNEST_SHA256" ]; then
      mv -f "${panel}.mirror" "$panel"
      from="the mirror"
    else
      rm -f "${panel}.mirror"
      echo "Data cache: ${name} not on the mirror or not matching its pin; setup.sh downloads it from EBI"
    fi
  fi
  echo "+ ANCESTRY_PANEL_NAME=${name} scripts/setup.sh --ancestry-panel (panel from: ${from:-EBI, by setup.sh})"
  if [ -n "$from" ]; then
    "${NET_BLOCKED[@]}" ANCESTRY_PANEL_NAME="$name" "${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
    check_eq "setup.sh --ancestry-panel exits 0 from the panel of ${from}, with EBI unreachable" "${PIPESTATUS[0]}" 0
    check "and downloads nothing" lacks 'Downloading' "$(cat "$STEP_LOG")"
  else
    ANCESTRY_PANEL_NAME="$name" "${REPO}/scripts/setup.sh" --ancestry-panel "$GENOME_DIR" 2>&1 | tee "$STEP_LOG"
    check_eq "setup.sh --ancestry-panel exits 0" "${PIPESTATUS[0]}" 0
    if [ "$name" = GRCh38_HAPNEST_reference ]; then
      check_eq "the downloaded panel matches its pin" "$(_sha256 "$panel")" "$HAPNEST_SHA256"
    fi
  fi
  if [ "$name" = GRCh38_HAPNEST_reference ] && [ "$from" != "the Actions cache" ] && [ -n "${E2E_DATA_CACHE:-}" ] \
     && [ "${E2E_HAPNEST_KEEP:-false}" = true ] && [ "$(_sha256 "$panel")" = "$HAPNEST_SHA256" ]; then
    mkdir -p "${cached%/*}"
    ln -f "$panel" "$cached" 2>/dev/null || cp -f "$panel" "$cached"
    echo "Data cache: the panel is in ${cached} for the Actions cache"
  fi
  echo "${CASE_NAME}: ancestry panel from ${from:-EBI}" >> "$E2E_NOTES"
}

# finish: exit 1 when a check failed, or when no check ran and the case gave
# no NO_CHECKS_REASON; scripts/ci/e2e-run.sh reads the last line it prints.
finish() {
  if [ "$FAILS" -gt 0 ]; then
    echo "${CASE_NAME}: ${FAILS} check(s) failed"
    exit 1
  fi
  if [ "$CHECKS" -eq 0 ]; then
    if [ -n "$NO_CHECKS_REASON" ]; then
      echo "${CASE_NAME}: no checks run: ${NO_CHECKS_REASON}"
      echo "${CASE_NAME}: no checks run: ${NO_CHECKS_REASON}" >> "$E2E_NOTES"
      exit 0
    fi
    echo "${CASE_NAME}: no check ran, so the case proves nothing (set NO_CHECKS_REASON when that is meant)"
    exit 1
  fi
  echo "${CASE_NAME}: all ${CHECKS} checks passed"
  exit 0
}
