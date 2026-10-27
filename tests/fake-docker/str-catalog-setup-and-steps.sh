#!/usr/bin/env bash
# The STRchive catalog of steps 09 and 9b: setup.sh installs it, step 09 runs
# on it by default, step 9b grades with the catalog step 09 ran on, and
# EH_CATALOG=bundled keeps the call of before. The fake download serves the
# two real STRchive files of versions.env, which this case fetches itself; a
# hook runs bin/strchive_convert.py for real on the host and plays Stranger.
#   1. A copy of the ExpansionHunter catalog with one byte flipped: setup.sh
#      --strchive-data exits 1, keeps neither that file nor a catalog, and
#      never runs the converter (the control that the sha256 check fails).
#   2. The real files: setup.sh --strchive-data writes catalog.json with 82
#      loci through bin/strchive_convert.py in PYTHON_IMAGE, and a second run
#      downloads nothing and converts nothing.
#   3. Step 09 with no EH_CATALOG runs ExpansionHunter on catalog.json, names
#      the loci it does not grade (MUC1) and copies the catalog beside the
#      VCF; step 9b grades with that copy and counts the motif of the
#      flagged FAME locus with it; STRANGER_REPEATS=bundled goes back to
#      Stranger's GRCh38 catalog.
#   4. EH_CATALOG=bundled: the ExpansionHunter call is the one step 09 makes
#      without the STRchive file (the call of before this change), byte for
#      byte, and leaves no catalog copy, so step 9b grades with Stranger's
#      bundled GRCh38 catalog.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook   # writes host-path.sh, which the hook below reads
DIR="${GENOME_DIR}/reference/strchive-${STRCHIVE_VERSION}"
EH_FILE="STRchive-disease-loci-${STRCHIVE_VERSION}.hg38.expansionhunter.json"
LOCI_FILE="STRchive-loci-${STRCHIVE_VERSION}.json"
SERVED="${CASE_WORK}/served" BAD="${CASE_WORK}/served-bad"
mkdir -p "$SERVED" "$BAD"
python3 - "$SERVED" "$STRCHIVE_VERSION" "$EH_FILE" "$LOCI_FILE" <<'PY'
import os, sys, time, urllib.request
out, rel = sys.argv[1], sys.argv[2]
for name in sys.argv[3:]:
    url = f"https://github.com/dashnowlab/STRchive/releases/download/{rel}/{name}"
    for i in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                open(os.path.join(out, name), "wb").write(r.read())
            break
        except OSError as e:
            print(f"download {i + 1}/3 of {url} failed: {e}")
            time.sleep(5 * (i + 1))
PY
sha() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi; }
[ "$(sha "${SERVED}/${EH_FILE}")" = "$STRCHIVE_EH_SHA256" ] || fail "the served ${EH_FILE} is not the pinned file"
[ "$(sha "${SERVED}/${LOCI_FILE}")" = "$STRCHIVE_LOCI_SHA256" ] || fail "the served ${LOCI_FILE} is not the pinned file"
cp "${SERVED}/${LOCI_FILE}" "${BAD}/"
python3 -c 'import sys; d = bytearray(open(sys.argv[1], "rb").read()); d[len(d) // 2] ^= 1; open(sys.argv[2], "wb").write(d)' \
  "${SERVED}/${EH_FILE}" "${BAD}/${EH_FILE}"

cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
case " $* " in
  *" /pgp-bin/strchive_convert.py "*)
    # The real converter, with the container paths mapped to the host.
    args=()
    for a in "$@"; do
      case "$a" in
        /pgp-bin/*) args+=("${REPO_ROOT}/bin/${a#/pgp-bin/}") ;;
        /genome*) args+=("$(host_path "$a")") ;;
        *) args+=("$a") ;;
      esac
    done
    exec "${args[@]}" ;;
  *" /pgp-bin/rfc1_motif.py "*" --flagged "*)
    prev=""
    for a in "$@"; do [ "$prev" = --out ] && out=$(host_path "$a"); prev=$a; done
    printf 'locus\tverdict\nSAMD12\tpathogenic_motif\n' > "$out" ;;
  *" stranger "*)
    # One FAME record Stranger flags, so step 9b counts its motif.
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
    printf 'chr8\t118366913\tSAMD12\tA\t<STR40>\t.\tPASS\tEND=118366918;VARID=SAMD12;REPID=SAMD12;STR_STATUS=pre_mutation\tGT\t0/1\n' ;;
  *) exec "${CASE_WORK}/hook-outputs" image "$@" ;;
esac
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"
fetches() { grep -cE "^(curl|wget) [^ ]*/STRchive-" "$FAKE_DOCKER_LOG" || true; }
converts() { grep -c 'strchive_convert\.py' "$FAKE_DOCKER_LOG" || true; }

# --- 1. a flipped byte: no file, no catalog, no conversion --------------------------------
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$BAD" run_expect 1 setup-bad "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup-bad "sha256 checksum of ${EH_FILE//./\\.} is "
output_has setup-bad "\[WARN\] Could not install STRchive's ${EH_FILE//./\\.}"
[ ! -e "${DIR}/${EH_FILE}" ] || fail "setup.sh kept the ExpansionHunter catalog that failed its sha256 check"
[ ! -e "${DIR}/catalog.json" ] || fail "setup.sh wrote a catalog from a file that failed its sha256 check"
[ "$(converts)" -eq 0 ] || fail "setup.sh ran the converter on a file that failed its sha256 check"

# --- 2. the real files: converted once ------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup "\[OK\] STRchive ${STRCHIVE_VERSION//./\\.} catalog \(steps 09 and 9b\): "
[ "$(fetches)" -eq 2 ] || fail "setup.sh downloaded the STRchive files $(fetches) times, expected 2"
docker_log_has "image=${PYTHON_IMAGE//./\\.}.* python3 /pgp-bin/strchive_convert\\.py " "setup.sh did not run the converter in PYTHON_IMAGE"
N=$(python3 -c 'import json, sys; print(len(json.load(open(sys.argv[1]))))' "${DIR}/catalog.json" 2>/dev/null || echo 0)
[ "$N" -eq 82 ] || fail "catalog.json has ${N} loci, expected 82"
: > "$FAKE_DOCKER_LOG"
FAKE_DOWNLOAD_DIR="$SERVED" run_expect 0 setup-again "${SCRIPTS}/setup.sh" --strchive-data "$GENOME_DIR"
output_has setup-again "already present"
[ "$(fetches)" -eq 0 ] || fail "the second setup.sh downloaded the STRchive files again"
[ "$(converts)" -eq 0 ] || fail "the second setup.sh converted again"

# --- 3. steps 09 and 9b on the STRchive catalog -------------------------------------------------
EH_DIR="${GENOME_DIR}/sample1/expansion_hunter"
eh_call() { grep -e '--variant-catalog' "$FAKE_DOCKER_LOG" | tail -n 1; }  # the logged docker run of ExpansionHunter
: > "$FAKE_DOCKER_LOG"
run_expect 0 eh-strchive "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-strchive "Variant catalog: STRchive ${STRCHIVE_VERSION//./\\.}, "
output_has eh-strchive "Not graded \(no call from this count\): MUC1: the reference repeat is 3734 bases"
docker_log_has "--variant-catalog /genome/reference/strchive-${STRCHIVE_VERSION//./\\.}/catalog\\.json " \
  "step 09 did not run ExpansionHunter on the STRchive catalog"
cmp -s "${DIR}/catalog.json" "${EH_DIR}/sample1_eh_catalog.json" || fail "step 09 did not copy the catalog it ran on beside the VCF"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' > "${EH_DIR}/sample1_eh.vcf"
: > "$FAKE_DOCKER_LOG"
run_expect 0 stranger-strchive "${SCRIPTS}/09b-stranger.sh" sample1
output_has stranger-strchive "Repeat catalog: the one step 09 ran on, .*sample1_eh_catalog\.json \(STRchive-${STRCHIVE_VERSION//./\\.}\)"
docker_log_has "stranger --repeats-file /genome/sample1/expansion_hunter/sample1_eh_catalog\\.json " \
  "step 9b did not grade with the catalog step 09 ran on"
docker_log_has "rfc1_motif\\.py .*--catalog /genome/sample1/expansion_hunter/sample1_eh_catalog\\.json --flagged" \
  "step 9b did not count the motif of the flagged FAME locus with the catalog"
grep -q $'^SAMD12\tpathogenic_motif' "${EH_DIR}/sample1_str_motifs.tsv" 2>/dev/null || fail "step 9b wrote no sample1_str_motifs.tsv"
rm -f "${EH_DIR}/sample1_eh_stranger.vcf"
: > "$FAKE_DOCKER_LOG"
STRANGER_REPEATS=bundled run_expect 0 stranger-bundled "${SCRIPTS}/09b-stranger.sh" sample1
docker_log_has 'variant_catalog_grch38\.json' "STRANGER_REPEATS=bundled did not pass Stranger's GRCh38 catalog"
[ ! -e "${EH_DIR}/sample1_str_motifs.tsv" ] || fail "with Stranger's bundled catalog step 9b kept the STRchive motif table"

# --- 4. EH_CATALOG=bundled: the call of before ---------------------------------------------------
mv "$DIR" "${DIR}.away"
: > "$FAKE_DOCKER_LOG"
run_expect 0 eh-before "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-before "Variant catalog: the 31 loci inside the image \(reference/strchive-"
BEFORE=$(eh_call)
mv "${DIR}.away" "$DIR"
: > "$FAKE_DOCKER_LOG"
EH_CATALOG=bundled run_expect 0 eh-bundled "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
output_has eh-bundled "Variant catalog: the 31 loci inside the image \(EH_CATALOG=bundled\)"
AFTER=$(eh_call)
echo "without STRchive:    ${BEFORE}"
echo "EH_CATALOG=bundled:  ${AFTER}"
[ -n "$BEFORE" ] && [ "$BEFORE" = "$AFTER" ] || fail "EH_CATALOG=bundled does not make the ExpansionHunter call of before"
case "$AFTER" in
  *"--variant-catalog /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json "*) ;;
  *) fail "EH_CATALOG=bundled does not pass the catalog inside the image: ${AFTER}" ;;
esac
[ ! -e "${EH_DIR}/sample1_eh_catalog.json" ] || fail "EH_CATALOG=bundled left the STRchive catalog copy beside the VCF"
rm -f "${EH_DIR}/sample1_eh_stranger.vcf"
: > "$FAKE_DOCKER_LOG"
run_expect 0 stranger-after-bundled "${SCRIPTS}/09b-stranger.sh" sample1
docker_log_has 'variant_catalog_grch38\.json' "after EH_CATALOG=bundled step 9b did not grade with Stranger's GRCh38 catalog"
