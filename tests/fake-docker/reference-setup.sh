#!/usr/bin/env bash
# setup.sh with its default reference: NCBI's GRCh38 no-ALT analysis set.
#   - it downloads the .fna.gz and the .fna.fai NCBI publishes, checks both
#     against the md5 NCBI lists for them in md5checksums.txt, unpacks the
#     FASTA to reference/GRCh38_no_alt_analysis_set.fasta and counts its
#     sequences;
#   - a download whose md5 is not the listed one is removed and setup fails,
#     leaving no reference;
#   - a checksum file that does not list the download fails setup too.
# The fake curl serves the files of FAKE_DOWNLOAD_DIR by base name.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# lib.sh points every other case at a made-up reference URL; this one tests
# the default.
unset REF_FASTA_URL REF_FASTA_MD5 REF_FAI_MD5

NAME=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
md5_of() { fake_md5 < "$1"; }

# serve DIR: a three-sequence FASTA as NCBI publishes it (gzip, plus its .fai)
# and an md5checksums.txt listing both, with other files around them.
serve() {
  mkdir -p "$1"
  printf '>chr1\nACGTACGTAC\nGT\n>chr2\nTTTTGGGG\n>chrM\nGATC\n' > "${CASE_WORK}/ref.fa"
  gzip -nc "${CASE_WORK}/ref.fa" > "${1}/${NAME}.gz"
  printf 'chr1\t12\t6\t10\t11\nchr2\t8\t25\t8\t9\nchrM\t4\t40\t4\t5\n' > "${1}/${NAME}.fai"
  {
    printf '%s  ./%s\n' "$(md5_of "${1}/${NAME}.gz")" "${NAME%.fna}_plus_hs38d1.fna.gz"
    printf '%s  ./%s\n' "$(md5_of "${1}/${NAME}.fai")" "${NAME}.fai"
    printf '%s  ./%s\n' "$(md5_of "${1}/${NAME}.gz")" "${NAME}.gz"
  } > "${1}/md5checksums.txt"
}

# --- the default reference ------------------------------------------------------
G="${CASE_WORK}/genome"
mkdir -p "$G"
serve "${CASE_WORK}/served"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/served" run_expect 0 setup "${SCRIPTS}/setup.sh" "$G"
docker_log_has "^curl ${NCBI//./\\.}/${NAME//./\\.}\\.gz -> " "setup.sh did not download the no-ALT analysis set from NCBI"
docker_log_has "^curl ${NCBI//./\\.}/md5checksums\\.txt -> " "setup.sh did not read NCBI's md5checksums.txt"
docker_log_has "^curl ${NCBI//./\\.}/${NAME//./\\.}\\.fai -> " "setup.sh did not download NCBI's .fai"
REF="${G}/reference/GRCh38_no_alt_analysis_set.fasta"
cmp -s "$REF" "${CASE_WORK}/ref.fa" || fail "${REF} is not the unpacked download"
cmp -s "${REF}.fai" "${CASE_WORK}/served/${NAME}.fai" || fail "${REF}.fai is not NCBI's .fai"
[ ! -e "${REF}.download.gz" ] || fail "setup.sh kept the compressed download"
output_has setup "\\[OK\\] Reference: ${REF//./\\.} \\(3 sequences\\)"

# --- a download whose md5 is not the listed one ---------------------------------------
G2="${CASE_WORK}/genome2"
mkdir -p "$G2"
serve "${CASE_WORK}/served-bad"
sed -i.orig "s|^[0-9a-f]*  ./${NAME}.gz\$|00000000000000000000000000000000  ./${NAME}.gz|" "${CASE_WORK}/served-bad/md5checksums.txt"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/served-bad" run_rc setup-bad-md5 "${SCRIPTS}/setup.sh" "$G2"
[ "$RC" -ne 0 ] || fail "setup.sh exited 0 although the reference did not match NCBI's md5"
output_has setup-bad-md5 "md5 checksum of ${NAME//./\\.}\\.gz is [0-9a-f]+, expected 0{32}"
[ -z "$(find "${G2}/reference" -name 'GRCh38_no_alt_analysis_set*' 2>/dev/null)" ] \
  || fail "a reference download with a wrong md5 was kept: $(find "${G2}/reference" -name 'GRCh38_no_alt_analysis_set*')"

# --- a checksum file that does not list the download -------------------------------------
G3="${CASE_WORK}/genome3"
mkdir -p "$G3"
serve "${CASE_WORK}/served-unlisted"
grep -v "  ./${NAME}.gz\$" "${CASE_WORK}/served-unlisted/md5checksums.txt" > "${CASE_WORK}/unlisted.txt"
mv "${CASE_WORK}/unlisted.txt" "${CASE_WORK}/served-unlisted/md5checksums.txt"
FAKE_DOWNLOAD_DIR="${CASE_WORK}/served-unlisted" run_rc setup-unlisted "${SCRIPTS}/setup.sh" "$G3"
[ "$RC" -ne 0 ] || fail "setup.sh exited 0 although md5checksums.txt does not list the reference"
output_has setup-unlisted "could not read the checksum of ${NAME//./\\.}\\.gz"
[ ! -e "${G3}/reference/GRCh38_no_alt_analysis_set.fasta" ] || fail "an unchecked reference was stored"

echo "reference-setup: the default reference is NCBI's no-ALT analysis set, unpacked and checked against md5checksums.txt"
