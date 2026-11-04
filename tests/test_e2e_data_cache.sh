#!/usr/bin/env bash
# The data cache entries hla and pgs of scripts/ci/e2e-run.sh: --data-check
# passes a whole entry and fails a changed, missing or partial one (the e2e
# workflow saves and uses only an entry that passes), and --data-key names
# the releases. No network, no docker.
set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
RUN="${REPO}/scripts/ci/e2e-run.sh"
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
failed=0
# expect DESC 0|1 NAME: --data-check NAME on $T exits as said.
expect() {
  local rc=0
  E2E_DATA_CACHE=$T "$RUN" --data-check "$3" > /dev/null 2>&1 || rc=$?
  if { [ "$2" = 0 ] && [ "$rc" -eq 0 ]; } || { [ "$2" = 1 ] && [ "$rc" -ne 0 ]; }; then echo "PASS: $1"
  else echo "FAIL: $1 (exit ${rc})"; failed=1; fi
}

# hla: four files and the SHA256SUMS data_fetch writes for them.
mkdir -p "${T}/hla/hla/IPD-IMGT-HLA_1/" "${T}/hla/hla/IPD-IMGT-HLA_2" "${T}/hla/reference"
for f in hla/IPD-IMGT-HLA_1/hla.dat hla/IPD-IMGT-HLA_2/hla.dat reference/genes.gtf reference/exons.bed; do
  echo "$f" > "${T}/hla/${f}"
done
(cd "${T}/hla" && find . -type f ! -name SHA256SUMS | LC_ALL=C sort | xargs sha256sum > SHA256SUMS)
expect "hla: a whole entry passes" 0 hla
echo changed >> "${T}/hla/reference/genes.gtf"
expect "hla: a changed file fails" 1 hla
sed -i.bak '$d' "${T}/hla/reference/genes.gtf" && rm -f "${T}/hla/reference/genes.gtf.bak"
expect "hla: the file put back passes again" 0 hla
sed -i.bak '$d' "${T}/hla/SHA256SUMS" && rm -f "${T}/hla/SHA256SUMS.bak"
expect "hla: a SHA256SUMS with three of the four files fails" 1 hla
rm -f "${T}/hla/SHA256SUMS"
expect "hla: no SHA256SUMS fails" 1 hla

# pgs: each score id of the case, with its .md5 beside it.
ids=$(sed -n 's/^NEW_IDS="\(.*\)"$/\1/p' "${REPO}/tests/e2e/prs-scores-and-run-options-1-strand-sex.sh")
[ -n "$ids" ] || { echo "FAIL: no NEW_IDS line in prs-scores-and-run-options-1-strand-sex.sh"; exit 1; }
mkdir -p "${T}/pgs"
for id in $ids; do
  echo "$id" | gzip -c > "${T}/pgs/${id}.txt.gz"
  echo "$(md5sum < "${T}/pgs/${id}.txt.gz" | cut -d' ' -f1)  ${id}_hmPOS_GRCh38.txt.gz" > "${T}/pgs/${id}.txt.gz.md5"
done
expect "pgs: every score with a matching .md5 passes" 0 pgs
last=${ids##* }
echo x | gzip -c > "${T}/pgs/${last}.txt.gz"
expect "pgs: a score that does not match its .md5 fails" 1 pgs
rm -f "${T}/pgs/${last}.txt.gz" "${T}/pgs/${last}.txt.gz.md5"
expect "pgs: a missing score fails" 1 pgs
expect "an unknown entry fails" 1 nosuch

key=$("$RUN" --data-key hla)
pin=$(sed -n 's/^HLA_DB_RELEASE=\${HLA_DB_RELEASE:-\(.*\)}$/\1/p' "${REPO}/scripts/lib/common.sh")
if [[ "$key" == e2e-data-hla-*-"${pin}"-gencode[0-9]*-* ]]; then echo "PASS: the hla key names the pin ${pin}: ${key}"
else echo "FAIL: the hla key does not name the pin ${pin}: ${key}"; failed=1; fi
if "$RUN" --data-key pgs | grep -qE '^e2e-data-pgs-[0-9a-f]{16}$'; then echo "PASS: the pgs key is a hash"
else echo "FAIL: the pgs key"; failed=1; fi
exit "$failed"
