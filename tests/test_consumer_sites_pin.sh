#!/usr/bin/env bash
# The pin of the step-40 site table: consumer_sites_digest in
# scripts/lib/common.sh and CONSUMER_SITES_TABLE_SHA256 in versions.env.
#   1. two gzip files of the same table written at different times differ as
#      files and get the same digest (the gzip header holds the write time);
#   2. the digest is the sha256 of the table's content;
#   3. a table with one row changed gets another digest;
#   4. the pin is a sha256, and setup.sh and validate-setup.sh both compare
#      the table with it.
# It first proves it can fail: a digest of the gzip file itself must fail
# case 1.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# The helpers alone, so the test needs neither GENOME_DIR nor docker.
for f in _digest consumer_sites_digest; do
  sed -n "/^${f}() {/,/^}/p" "${REPO}/scripts/lib/common.sh"
done > "${WORK}/helpers.sh"
for f in '_digest() {' 'consumer_sites_digest() {'; do
  grep -qF -- "$f" "${WORK}/helpers.sh" || { echo "FAIL: no ${f%% *} in scripts/lib/common.sh"; exit 1; }
done
# shellcheck source=/dev/null
. "${WORK}/helpers.sh"

BAD=0
ok()  { echo "ok   $*"; }
bad() { echo "FAIL $*"; BAD=$((BAD + 1)); }

# gz TEXT_FILE OUT MTIME: OUT holds TEXT_FILE gzipped, with MTIME in its header
gz() {
  python3 -c 'import gzip, sys; gzip.GzipFile(sys.argv[2], "wb", mtime=int(sys.argv[3])).write(open(sys.argv[1], "rb").read())' "$1" "$2" "$3"
}
sha() { _digest sha256 "$1"; }

printf '## consumer sites: test\n#rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\nrs1\t1\t100\tA\t200\tA\t+\nrs2\t2\t300\tC\t400\tG\t-\n' > "${WORK}/table.tsv"
gz "${WORK}/table.tsv" "${WORK}/a.tsv.gz" 1
gz "${WORK}/table.tsv" "${WORK}/b.tsv.gz" 2
awk -F'\t' -v OFS='\t' '$1 == "rs2" { $3 = 301 } 1' "${WORK}/table.tsv" > "${WORK}/changed.tsv"
gz "${WORK}/changed.tsv" "${WORK}/c.tsv.gz" 1

# same_table DIGEST_FUNCTION: case 1 against that function
same_table() {
  [ "$(sha "${WORK}/a.tsv.gz")" != "$(sha "${WORK}/b.tsv.gz")" ] || return 2
  [ "$("$1" "${WORK}/a.tsv.gz")" = "$("$1" "${WORK}/b.tsv.gz")" ]
}

# Red control: hashing the gzip file itself sees the write time.
file_digest() { sha "$1"; }
if same_table file_digest; then
  bad "control: a digest of the gzip file passed case 1, so case 1 cannot fail"
else
  ok "control: a digest of the gzip file fails case 1"
fi

rc=0; same_table consumer_sites_digest || rc=$?
case $rc in
  0) ok "1. the same table gzipped at two times gets one digest" ;;
  2) bad "1. the two gzip files are equal, so the case tests nothing" ;;
  *) bad "1. the same table gzipped at two times got two digests" ;;
esac

want=$(sha "${WORK}/table.tsv")
got=$(consumer_sites_digest "${WORK}/a.tsv.gz")
if [ "$got" = "$want" ]; then ok "2. the digest is the sha256 of the content"; else bad "2. digest ${got}, content sha256 ${want}"; fi

if [ "$(consumer_sites_digest "${WORK}/c.tsv.gz")" != "$got" ]; then
  ok "3. a changed row changes the digest"
else
  bad "3. a changed row kept the digest"
fi

pin=$(sed -n 's/^CONSUMER_SITES_TABLE_SHA256="\(.*\)"$/\1/p' "${REPO}/versions.env")
if [[ "$pin" =~ ^[0-9a-f]{64}$ ]]; then ok "4. CONSUMER_SITES_TABLE_SHA256 is a sha256"; else bad "4. CONSUMER_SITES_TABLE_SHA256 is '${pin}', not a sha256"; fi
for s in setup.sh validate-setup.sh; do
  if grep -q 'consumer_sites_digest' "${REPO}/scripts/${s}" && grep -q 'CONSUMER_SITES_TABLE_SHA256' "${REPO}/scripts/${s}"; then
    ok "4. scripts/${s} compares the table with the pin"
  else
    bad "4. scripts/${s} does not compare the table with the pin"
  fi
done

if [ "$BAD" -eq 0 ]; then echo "RESULT: ALL PASS"; else echo "RESULT: ${BAD} FAILED"; exit 1; fi
