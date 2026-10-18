#!/usr/bin/env bash
# check-coupled-versions.sh: fail when a data version in versions.env no longer
# matches the image it belongs to.
#
# Some data releases only work with one image release, so a bump of one
# without the other breaks a step at run time, far from the bump:
#   VEP_CACHE_RELEASE            must equal the major of VEP_IMAGE (release_116.0 -> 116)
#   PYPGX_BUNDLE_VERSION         must equal the version of PYPGX_IMAGE (0.26.0--pyh... -> 0.26.0)
#   ANNOTSV_ANNOTATIONS_VERSION  must equal major.minor of ANNOTSV_IMAGE (3.5.10--h... -> 3.5)
#   PCGR_IMAGE, PCGR_VEP_CACHE_RELEASE and PCGR_DATA_BUNDLE must all be set
#   YLEAF_DATA_VERSION           must equal the version of YLEAF_IMAGE (3.2.1--pyh... -> 3.2.1)
#   GENCODE_RELEASE              (scripts/lib/common.sh) must be the GENCODE release of the
#                                Ensembl release of VEP_IMAGE, its major minus 66 (116 -> 50)
#   VB2_PANEL_FILES              (scripts/setup.sh) must name the VerifyBamID tag of the
#                                version of VERIFYBAMID2_IMAGE (2.0.3--hc004090_0 -> v2.0.3)
#
# Usage:
#   scripts/ci/check-coupled-versions.sh [FILE]     check FILE (default: versions.env),
#                                                   with scripts/lib/common.sh and scripts/setup.sh
#   scripts/ci/check-coupled-versions.sh --self-test
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# check FILE [COMMON SETUP]: print one line per rule, exit 1 if one fails. FILE
# is sourced in a clean shell, so a variable exported by the caller cannot
# stand in for it. COMMON and SETUP (default: the repository's) are read as
# text: sourcing them would run them.
check() {
  local vals fail=0 common=${2:-${ROOT}/scripts/lib/common.sh} setup=${3:-${ROOT}/scripts/setup.sh}
  # shellcheck disable=SC2016  # $1 and ${!v} expand in the inner shell
  vals=$(env -i bash --noprofile --norc -c '. "$1"; for v in VEP_IMAGE VEP_CACHE_RELEASE PYPGX_IMAGE PYPGX_BUNDLE_VERSION \
      ANNOTSV_IMAGE ANNOTSV_ANNOTATIONS_VERSION PCGR_IMAGE PCGR_VEP_CACHE_RELEASE PCGR_DATA_BUNDLE \
      YLEAF_IMAGE YLEAF_DATA_VERSION VERIFYBAMID2_IMAGE; do
      printf "%s=%s\n" "$v" "${!v-}"; done' _ "$1") || { echo "FAIL: cannot read $1"; return 1; }
  get() { sed -n "s/^$1=//p" <<<"$vals"; }
  tag() { local i; i=$(get "$1"); i=${i##*/}; [[ "$i" == *:* ]] && echo "${i#*:}"; }
  # rule NAME HAVE WANT
  rule() {
    if [ -n "$3" ] && [ "$2" = "$3" ]; then
      echo "OK   $1: $2"
    else
      echo "FAIL $1: ${2:-(unset)}, want ${3:-(cannot be derived)}"
      fail=1
    fi
  }
  local t
  t=$(tag VEP_IMAGE || true); t=${t#release_}
  rule "VEP_CACHE_RELEASE matches the major of VEP_IMAGE" "$(get VEP_CACHE_RELEASE)" "${t%%.*}"
  t=$(tag PYPGX_IMAGE || true)
  rule "PYPGX_BUNDLE_VERSION matches the version of PYPGX_IMAGE" "$(get PYPGX_BUNDLE_VERSION)" "${t%%--*}"
  t=$(tag ANNOTSV_IMAGE || true); t=${t%%--*}
  rule "ANNOTSV_ANNOTATIONS_VERSION matches major.minor of ANNOTSV_IMAGE" "$(get ANNOTSV_ANNOTATIONS_VERSION)" \
    "$(sed -nE 's/^([0-9]+\.[0-9]+).*/\1/p' <<<"$t")"
  t=$(tag YLEAF_IMAGE || true)
  rule "YLEAF_DATA_VERSION matches the version of YLEAF_IMAGE" "$(get YLEAF_DATA_VERSION)" "${t%%--*}"
  t=$(tag VEP_IMAGE || true); t=${t#release_}; t=${t%%.*}
  rule "GENCODE_RELEASE (scripts/lib/common.sh) is the GENCODE release of VEP_IMAGE's Ensembl release" \
    "$(sed -nE 's/^GENCODE_RELEASE=(\$\{GENCODE_RELEASE:-)?([0-9]+)\}?.*/\2/p' "$common" 2>/dev/null | head -n1)" \
    "$( [[ "$t" =~ ^[0-9]+$ ]] && echo $((t - 66)))"
  # Every panel file must come from the same tag, so one stray URL fails too.
  local tags
  tags=$(sed -nE 's|.*raw\.githubusercontent\.com/Griffan/VerifyBamID/([^/]+)/resource/.*|\1|p' "$setup" 2>/dev/null | sort -u | paste -sd, -)
  t=$(tag VERIFYBAMID2_IMAGE || true); t=${t%%--*}
  rule "VB2_PANEL_FILES (scripts/setup.sh) come from the VerifyBamID tag of VERIFYBAMID2_IMAGE" "$tags" "${t:+v$t}"
  local v missing=""
  for v in PCGR_IMAGE PCGR_VEP_CACHE_RELEASE PCGR_DATA_BUNDLE; do
    [ -n "$(get "$v")" ] || missing="${missing} ${v}"
  done
  if [ -z "$missing" ]; then
    echo "OK   PCGR_IMAGE, PCGR_VEP_CACHE_RELEASE and PCGR_DATA_BUNDLE are set together"
  else
    echo "FAIL PCGR variables must be set together; unset:${missing}"
    fail=1
  fi
  return "$fail"
}

# self_test: each broken copy of versions.env must fail with its rule named,
# and the unchanged file must pass.
self_test() {
  local tmp fail=0 rc out
  tmp=$(mktemp -d)
  # shellcheck disable=SC2064
  trap "rm -rf '$tmp'" EXIT
  # expect NAME SED PATTERN [FILE]: FILE (versions.env by default, or common or
  # setup) is the copy the sed edits; the other two are the repository's.
  expect() {
    local f=${4:-env} env=${ROOT}/versions.env common=${ROOT}/scripts/lib/common.sh setup=${ROOT}/scripts/setup.sh
    case $f in
      env) sed -E "$2" "$env" > "${tmp}/$1.env"; env=${tmp}/$1.env ;;
      common) sed -E "$2" "$common" > "${tmp}/$1.common"; common=${tmp}/$1.common ;;
      setup) sed -E "$2" "$setup" > "${tmp}/$1.setup"; setup=${tmp}/$1.setup ;;
    esac
    rc=0; out=$(check "$env" "$common" "$setup" 2>&1) || rc=$?
    if [ "$rc" -ne 1 ] || ! grep -qE -- "$3" <<<"$out"; then
      echo "self-test: '$1' exited ${rc} and did not report /$3/:"; printf '%s\n' "$out"; fail=1
    else
      echo "self-test: '$1' caught: $(grep -E -- "$3" <<<"$out")"
    fi
  }
  rc=0; out=$(check "${ROOT}/versions.env" 2>&1) || rc=$?
  [ "$rc" -eq 0 ] || { echo "self-test: versions.env itself failed:"; printf '%s\n' "$out"; fail=1; }
  # The values the unchanged file holds, so a legitimate bump keeps this
  # test green; each planted value differs from them.
  local vep pypgx annotsv yleaf vb2
  vep=$(sed -nE 's/^VEP_IMAGE="[^"]*:release_([0-9]+).*/\1/p' "${ROOT}/versions.env")
  pypgx=$(sed -nE 's/^PYPGX_IMAGE="[^"]*:([^"]*)--.*/\1/p' "${ROOT}/versions.env")
  annotsv=$(sed -nE 's/^ANNOTSV_IMAGE="[^"]*:([0-9]+\.[0-9]+).*/\1/p' "${ROOT}/versions.env")
  yleaf=$(sed -nE 's/^YLEAF_IMAGE="[^"]*:([^"]*)--.*/\1/p' "${ROOT}/versions.env")
  vb2=$(sed -nE 's/^VERIFYBAMID2_IMAGE="[^"]*:([^"]*)--.*/\1/p' "${ROOT}/versions.env")
  if [ -z "$vep" ] || [ -z "$pypgx" ] || [ -z "$annotsv" ] || [ -z "$yleaf" ] || [ -z "$vb2" ]; then
    echo "self-test: cannot read VEP_IMAGE, PYPGX_IMAGE, ANNOTSV_IMAGE, YLEAF_IMAGE or VERIFYBAMID2_IMAGE from versions.env"; return 1
  fi
  expect vep-cache "s/^VEP_CACHE_RELEASE=.*/VEP_CACHE_RELEASE=\"$((vep - 1))\"/" "^FAIL VEP_CACHE_RELEASE .*: $((vep - 1)), want ${vep}\$"
  expect vep-image "s/^(VEP_IMAGE=\"[^\"]*:release_)[0-9]+/\\1$((vep + 1))/" "^FAIL VEP_CACHE_RELEASE .*: ${vep}, want $((vep + 1))\$"
  expect pypgx 's/^PYPGX_BUNDLE_VERSION=.*/PYPGX_BUNDLE_VERSION="9.9.9"/' "^FAIL PYPGX_BUNDLE_VERSION .*: 9\\.9\\.9, want ${pypgx//./\\.}\$"
  expect annotsv 's/^ANNOTSV_ANNOTATIONS_VERSION=.*/ANNOTSV_ANNOTATIONS_VERSION="0.1"/' "^FAIL ANNOTSV_ANNOTATIONS_VERSION .*: 0\\.1, want ${annotsv//./\\.}\$"
  expect pcgr '/^PCGR_DATA_BUNDLE=/d' '^FAIL PCGR variables must be set together; unset: PCGR_DATA_BUNDLE'
  expect yleaf 's/^YLEAF_DATA_VERSION=.*/YLEAF_DATA_VERSION="9.9.9"/' "^FAIL YLEAF_DATA_VERSION .*: 9\\.9\\.9, want ${yleaf//./\\.}\$"
  expect gencode "s/^(GENCODE_RELEASE=)[^ ]*/\\1$((vep - 67))/" "^FAIL GENCODE_RELEASE .*: $((vep - 67)), want $((vep - 66))\$" common
  expect gencode-vep "s/^(VEP_IMAGE=\"[^\"]*:release_)[0-9]+/\\1$((vep + 1))/" "^FAIL GENCODE_RELEASE .*: $((vep - 66)), want $((vep - 65))\$"
  # One of the three panel URLs moved to another tag.
  expect vb2 's|(VerifyBamID/)v[0-9.]+(/resource/[^ ]*\.dat\.UD )|\1v9.9.9\2|' \
    "^FAIL VB2_PANEL_FILES .*: (v${vb2//./\\.},v9\\.9\\.9|v9\\.9\\.9,v${vb2//./\\.}), want v${vb2//./\\.}\$" setup
  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -*) echo "usage: $0 [FILE | --self-test]" >&2; exit 2 ;;
  *) check "${1:-${ROOT}/versions.env}" ;;
esac
