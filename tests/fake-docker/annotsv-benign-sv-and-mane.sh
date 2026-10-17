#!/usr/bin/env bash
# setup.sh --annotsv-data on a GENOME_DIR that already holds AnnotSV's bundle:
# it downloads gnomAD-SV's sites and NCBI's MANE summary from the URLs of the
# versions.env pins and checks each against its sha256. The fake downloads are
# placeholders, so both checks fail: setup must say so for each file, exit
# non-zero and leave neither the files AnnotSV reads nor their stamps. Without
# the bundle it downloads the bundle first and stops when that fails.
# (tests/test_annotsv_benign_mane.sh runs the same function on real files.)
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

G="${CASE_WORK}/genome"
A="${G}/annotsv_annotations"
B="${A}/Annotations_Human/SVincludedInFt/BenignSV/GRCh38"
mkdir -p "${A}/Annotations_Human/Genes/GRCh38"

: > "$FAKE_DOCKER_LOG"
run_expect 1 annotsv-data "${SCRIPTS}/setup.sh" --annotsv-data "$G"
output_has annotsv-data '\[OK\] AnnotSV annotation data already present'
v=${GNOMAD_SV_VERSION//./\\.}
docker_log_has "^curl https://storage\\.googleapis\\.com/gcp-public-data--gnomad/release/${v}/genome_sv/gnomad\\.v${v}\\.sv\\.sites\\.vcf\\.gz " \
  "setup.sh did not download gnomAD-SV ${GNOMAD_SV_VERSION}'s sites file from Google's gnomAD bucket"
m=${MANE_RELEASE//./\\.}
docker_log_has "^curl https://ftp\\.ncbi\\.nlm\\.nih\\.gov/refseq/MANE/MANE_human/release_${m}/MANE\\.GRCh38\\.v${m}\\.summary\\.txt\\.gz " \
  "setup.sh did not download the MANE ${MANE_RELEASE} summary from NCBI"
output_has annotsv-data "sha256 checksum of gnomad\\.v${v}\\.sv\\.sites\\.vcf\\.gz is [0-9a-f]+, expected ${GNOMAD_SV_SHA256}"
output_has annotsv-data "sha256 checksum of MANE\\.GRCh38\\.v${m}\\.summary\\.txt\\.gz is [0-9a-f]+, expected ${MANE_SUMMARY_SHA256}"
output_has annotsv-data "\\[WARN\\] Could not download gnomAD-SV ${v}"
output_has annotsv-data "\\[WARN\\] Could not download the MANE ${m} summary"
for f in "${B}/gnomad-sv.installed" "${A}/mane_select_refseq.txt" "${A}/mane_select_refseq.release"; do
  [ ! -e "$f" ] || fail "a failed checksum left ${f}"
done
left=$(find "$A" -name '*.part' -o -name 'gnomad.v4.0.sv.chr*' | head -n 1)
[ -z "$left" ] || fail "a failed checksum left ${left}"

# No bundle: setup tries the bundle and stops there.
rm -rf "$A"
: > "$FAKE_DOCKER_LOG"
run_expect 1 annotsv-data-no-bundle "${SCRIPTS}/setup.sh" --annotsv-data "$G"
docker_log_has '^curl https://www\.lbgi\.fr/~geoffroy/Annotations/Annotations_Human_' "setup.sh did not download AnnotSV's bundle"
output_has annotsv-data-no-bundle 'Step 5 \(AnnotSV\) will be skipped'
if grep -qE 'gnomad\.v|MANE\.GRCh38' "$FAKE_DOCKER_LOG"; then
  fail "setup.sh fetched gnomAD-SV or MANE without a bundle to add them to"
fi
