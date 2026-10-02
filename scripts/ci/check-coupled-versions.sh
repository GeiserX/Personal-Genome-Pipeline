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
#
# Usage:
#   scripts/ci/check-coupled-versions.sh [FILE]     check FILE (default: versions.env)
#   scripts/ci/check-coupled-versions.sh --self-test
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# check FILE: print one line per rule, exit 1 if one fails. FILE is sourced in
# a clean shell, so a variable exported by the caller cannot stand in for it.
check() {
  local vals fail=0
  # shellcheck disable=SC2016  # $1 and ${!v} expand in the inner shell
  vals=$(env -i bash --noprofile --norc -c '. "$1"; for v in VEP_IMAGE VEP_CACHE_RELEASE PYPGX_IMAGE PYPGX_BUNDLE_VERSION \
      ANNOTSV_IMAGE ANNOTSV_ANNOTATIONS_VERSION PCGR_IMAGE PCGR_VEP_CACHE_RELEASE PCGR_DATA_BUNDLE; do
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
  # expect NAME SED PATTERN
  expect() {
    sed -E "$2" "${ROOT}/versions.env" > "${tmp}/$1.env"
    rc=0; out=$(check "${tmp}/$1.env" 2>&1) || rc=$?
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
  local vep pypgx annotsv
  vep=$(sed -nE 's/^VEP_IMAGE="[^"]*:release_([0-9]+).*/\1/p' "${ROOT}/versions.env")
  pypgx=$(sed -nE 's/^PYPGX_IMAGE="[^"]*:([^"]*)--.*/\1/p' "${ROOT}/versions.env")
  annotsv=$(sed -nE 's/^ANNOTSV_IMAGE="[^"]*:([0-9]+\.[0-9]+).*/\1/p' "${ROOT}/versions.env")
  if [ -z "$vep" ] || [ -z "$pypgx" ] || [ -z "$annotsv" ]; then
    echo "self-test: cannot read VEP_IMAGE, PYPGX_IMAGE or ANNOTSV_IMAGE from versions.env"; return 1
  fi
  expect vep-cache "s/^VEP_CACHE_RELEASE=.*/VEP_CACHE_RELEASE=\"$((vep - 1))\"/" "^FAIL VEP_CACHE_RELEASE .*: $((vep - 1)), want ${vep}\$"
  expect vep-image "s/^(VEP_IMAGE=\"[^\"]*:release_)[0-9]+/\\1$((vep + 1))/" "^FAIL VEP_CACHE_RELEASE .*: ${vep}, want $((vep + 1))\$"
  expect pypgx 's/^PYPGX_BUNDLE_VERSION=.*/PYPGX_BUNDLE_VERSION="9.9.9"/' "^FAIL PYPGX_BUNDLE_VERSION .*: 9\\.9\\.9, want ${pypgx//./\\.}\$"
  expect annotsv 's/^ANNOTSV_ANNOTATIONS_VERSION=.*/ANNOTSV_ANNOTATIONS_VERSION="0.1"/' "^FAIL ANNOTSV_ANNOTATIONS_VERSION .*: 0\\.1, want ${annotsv//./\\.}\$"
  expect pcgr '/^PCGR_DATA_BUNDLE=/d' '^FAIL PCGR variables must be set together; unset: PCGR_DATA_BUNDLE'
  [ "$fail" -eq 0 ] && echo "self-test: OK"
  return "$fail"
}

case "${1:-}" in
  --self-test) self_test ;;
  -*) echo "usage: $0 [FILE | --self-test]" >&2; exit 2 ;;
  *) check "${1:-${ROOT}/versions.env}" ;;
esac
