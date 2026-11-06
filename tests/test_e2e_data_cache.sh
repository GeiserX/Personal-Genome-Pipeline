#!/usr/bin/env bash
# The data cache entries hla and pgs of scripts/ci/e2e-run.sh: --data-check
# passes a whole entry and fails a changed, missing or partial one (the e2e
# workflow saves and uses only an entry that passes), and --data-key names
# the releases; and .github/actions/free-disk names a cleanup that failed in
# its summary (fake sudo, docker and df). No network, no docker.
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
other=$(sed -n 's/^OTHER=//p' "${REPO}/tests/e2e/bash-step-08-hla.sh")
gencode=$(sed -n 's/^GENCODE_RELEASE=//p' "${REPO}/scripts/lib/common.sh")
if [ -n "$other" ] && [ -n "$pin" ] && [ -n "$gencode" ] \
   && [[ "$key" =~ ^e2e-data-hla-${other//./\\.}-${pin//./\\.}-gencode${gencode}-[0-9a-f]{16}$ ]]; then
  echo "PASS: the hla key names ${other}, ${pin}, GENCODE ${gencode} and a hash: ${key}"
else echo "FAIL: the hla key does not match its inputs (${other:-?} ${pin:-?} ${gencode:-?}): ${key}"; failed=1; fi
if "$RUN" --data-key pgs | grep -qE '^e2e-data-pgs-[0-9a-f]{16}$'; then echo "PASS: the pgs key is a hash"
else echo "FAIL: the pgs key"; failed=1; fi

# --- free-disk: a cleanup that fails is named in the summary -----------------
# The action's script, run as GitHub runs a composite bash step, with sudo,
# docker and df faked on PATH: sudo fails for the planted path, docker's prune
# fails in the second run. The control run, with nothing failing, must not
# name a failure.
FD="${T}/free-disk"
mkdir -p "${FD}/bin" "${FD}/ok" "${FD}/planted-fail"
awk '/^      run: \|/ {on = 1; next} on && /^        / {sub(/^        /, ""); print; next} on && NF {exit}' \
  "${REPO}/.github/actions/free-disk/action.yml" > "${FD}/run.sh"
grep -q 'failed' "${FD}/run.sh" || { echo "FAIL: could not read the free-disk script"; failed=1; }
cat > "${FD}/bin/sudo" <<'EOF'
#!/usr/bin/env bash
case "$*" in
  *planted-fail*) echo "rm: cannot remove: Permission denied" >&2; exit 1 ;;
  *"${FD_ROOT:?}"/*) exec "$@" ;;
  *) echo "fake sudo: refused, outside the test folder: $*" >&2; exit 9 ;;
esac
EOF
cat > "${FD}/bin/docker" <<'EOF'
#!/usr/bin/env bash
[ "${FAKE_PRUNE_FAILS:-false}" = true ] && exit 3
exit 0
EOF
cat > "${FD}/bin/df" <<'EOF'
#!/usr/bin/env bash
printf 'Avail\n10485760\n'
EOF
chmod +x "${FD}/bin/"*
freedisk() {  # freedisk PATHS [PRUNE_FAILS]: the summary the script wrote
  : > "${FD}/summary"
  FD_ROOT=$FD PATH="${FD}/bin:${PATH}" FREE_DISK_PATHS="$1" FAKE_PRUNE_FAILS="${2:-false}" GITHUB_STEP_SUMMARY="${FD}/summary" \
    bash --noprofile --norc -eo pipefail "${FD}/run.sh" > "${FD}/log" 2>&1
  echo "exit $?: $(cat "${FD}/summary")"
}
out=$(freedisk "${FD}/ok")
if [[ "$out" == "exit 0: Free disk:"* ]] && [[ "$out" != *failed* ]] && [ ! -e "${FD}/ok" ]; then
  echo "PASS: free-disk control: nothing failed, nothing named (${out})"
else echo "FAIL: free-disk control: ${out}"; failed=1; fi
mkdir -p "${FD}/ok"
out=$(freedisk "${FD}/ok ${FD}/planted-fail" true)
if [[ "$out" == "exit 0: "* ]] && [[ "$out" == *"cleanup failed"*"rm -rf ${FD}/planted-fail (exit 1)"* ]] \
   && [[ "$out" == *"docker image prune -a (exit 3)"* ]] && [[ "$out" != *"rm -rf ${FD}/ok "* ]]; then
  echo "PASS: free-disk names the planted failing rm and the failed prune, and still exits 0 (${out})"
else echo "FAIL: free-disk did not name the planted failures: ${out}"; failed=1; fi

exit "$failed"
