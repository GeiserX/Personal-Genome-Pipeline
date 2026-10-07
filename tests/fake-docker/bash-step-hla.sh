#!/usr/bin/env bash
# Step 08 (T1K HLA typing) on a pinned IPD-IMGT/HLA release:
#   - without the release installed by setup.sh it is skipped, with the reason;
#   - it builds its index from the installed hla.dat and GENCODE gene lines
#     into t1k_idx/t1k-<version>_imgt-<release>_gencode-<release>/, and a
#     second run reuses it;
#   - it writes database_release.txt with the release line of hla.dat;
#   - another HLA_DB_RELEASE builds another index;
#   - an index whose coordinate file has "-1 -1" for a typed gene is refused
#     and not kept (T1K would extract no reads for that gene);
#   - THREADS and ALIGN_DIR reach the container (also for steps 09 and 10),
#     and step 10 passes the GRCh38 bands once they are installed.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/sample1/aligned_bwamem2"
cp "${GENOME_DIR}/sample1/aligned/"* "${GENOME_DIR}/sample1/aligned_bwamem2/"
use_output_hook
cat > "${CASE_WORK}/t1k-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
case "$args" in
  *" t1k-build.pl "*)
    echo "t1k-build.pl" >> "${CASE_WORK}/t1k-builds"
    o=$(host_path "$(sed -n 's/.* -o \([^ ]*\) .*/\1/p' <<<"$args")")
    mkdir -p "$o"
    printf '>HLA-A*01:01:01:01\nACGT\n' > "${o}/hla_dna_seq.fa"
    : > "${o}/hla_dna_coord.fa"
    for g in HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1; do
      if [ "$g" = "${FAKE_NO_COORD:-}" ]; then
        printf '>%s*01:01 chr6 -1 -1 +\nACGT\n' "$g" >> "${o}/hla_dna_coord.fa"
      else
        printf '>%s*01:01 chr6 29942532 29945870 +\nACGT\n' "$g" >> "${o}/hla_dna_coord.fa"
      fi
    done
    exit 0 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/t1k-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/t1k-hook"

# install_hla RELEASE: what setup.sh installs for that release.
install_hla() {
  mkdir -p "${GENOME_DIR}/hla/IPD-IMGT-HLA_$1"
  printf 'ID   HLA00001; SV 4; standard; DNA; HUM; 3503 BP.\nCC   IPD-IMGT/HLA Release Version %s\n' "$1" \
    > "${GENOME_DIR}/hla/IPD-IMGT-HLA_$1/hla.dat"
}
GENES="${GENOME_DIR}/reference/gencode.v50.basic.genes.gtf"
grep -q '^GENCODE_RELEASE=50$' "${REPO_ROOT}/scripts/lib/common.sh" \
  || fail "common.sh no longer pins GENCODE 50; update this case"

# --- not installed: skipped, not failed -------------------------------------------------
run_expect 0 hla-not-installed "${SCRIPTS}/08-hla-typing.sh" sample1
output_has hla-not-installed 'SKIPPED: IPD-IMGT/HLA [0-9.]+ is not installed'
[ ! -e "${CASE_WORK}/t1k-builds" ] || fail "step 08 built an index without an installed database"

# --- installed: built once, reused, release recorded ---------------------------------------
install_hla 3.65.0
printf 'chr6\tHAVANA\tgene\t29941260\t29949572\t.\t+\t.\tgene_name "HLA-A";\n' > "$GENES"
: > "$FAKE_DOCKER_LOG"
THREADS=2 ALIGN_DIR=aligned_bwamem2 run_expect 0 hla "${SCRIPTS}/08-hla-typing.sh" sample1
IDX=$(find "${GENOME_DIR}/t1k_idx" -mindepth 1 -maxdepth 1 -type d -name 't1k-*_imgt-3.65.0_gencode-50')
[ -n "$IDX" ] || fail "step 08 built no index directory keyed on T1K version and release 3.65.0"
docker_log_has '^run image=[^ ]*t1k.* t1k-build\.pl -d /genome/hla/IPD-IMGT-HLA_3\.65\.0/hla\.dat -g /genome/reference/gencode\.v50\.basic\.genes\.gtf ' \
  "step 08 did not build from the installed hla.dat and GENCODE genes"
docker_log_has '^run image=[^ ]*t1k.* --cpus 2 .* run-t1k -b /genome/sample1/aligned_bwamem2/sample1_sorted\.bam .* -t 2 ' \
  "step 08 did not use THREADS=2 and ALIGN_DIR=aligned_bwamem2"
REL="${GENOME_DIR}/sample1/hla_t1k/database_release.txt"
grep -q '^database: IPD-IMGT/HLA Release Version 3\.65\.0$' "$REL" 2>/dev/null \
  || fail "database_release.txt does not name release 3.65.0: $(cat "$REL" 2>/dev/null)"
THREADS=2 ALIGN_DIR=aligned_bwamem2 run_expect 0 hla-again "${SCRIPTS}/08-hla-typing.sh" sample1
[ "$(grep -c . "${CASE_WORK}/t1k-builds")" -eq 1 ] || fail "step 08 rebuilt an index that already existed"

# --- another release: another index -------------------------------------------------------
install_hla 3.64.0
HLA_DB_RELEASE=3.64.0 run_expect 0 hla-other-release "${SCRIPTS}/08-hla-typing.sh" sample1
[ "$(grep -c . "${CASE_WORK}/t1k-builds")" -eq 2 ] || fail "a new HLA_DB_RELEASE did not build a new index"
[ -n "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -type d -name 't1k-*_imgt-3.64.0_gencode-50')" ] \
  || fail "no index directory for release 3.64.0"
grep -q '3\.64\.0' "$REL" || fail "database_release.txt still names the old release"

# --- a typed gene without coordinates -----------------------------------------------------
install_hla 3.63.0
HLA_DB_RELEASE=3.63.0 FAKE_NO_COORD=HLA-B run_rc hla-no-coord "${SCRIPTS}/08-hla-typing.sh" sample1
[ "$RC" -ne 0 ] || fail "step 08 accepted an index where HLA-B has no coordinates"
output_has hla-no-coord 'genes without GRCh38 coordinates in .*: HLA-B'
[ -z "$(find "${GENOME_DIR}/t1k_idx" -maxdepth 1 -name 't1k-*_imgt-3.63.0_gencode-50')" ] \
  || fail "step 08 kept the index with a gene without coordinates"

# --- steps 09 and 10: THREADS and ALIGN_DIR; step 10: the GRCh38 bands --------------------
: > "$FAKE_DOCKER_LOG"
THREADS=2 ALIGN_DIR=aligned_bwamem2 run_rc eh "${SCRIPTS}/09-expansion-hunter.sh" sample1 male
docker_log_has '^run image=[^ ]*expansionhunter.* --cpus 2 .* --reads /genome/sample1/aligned_bwamem2/sample1_sorted\.bam ' \
  "step 09 did not use THREADS=2 and ALIGN_DIR=aligned_bwamem2"
THREADS=2 ALIGN_DIR=aligned_bwamem2 run_expect 0 telomere-no-bands "${SCRIPTS}/10-telomere-hunter.sh" sample1
output_has telomere-no-bands 'GRCh38 chromosome bands are not installed'
docker_log_has '^run image=[^ ]*telomerehunter.* --cpus 2 .* -ibt /genome/sample1/aligned_bwamem2/sample1_sorted\.bam ' \
  "step 10 did not use THREADS=2 and ALIGN_DIR=aligned_bwamem2"
if awk '/^run / && /telomerehunter/ && / -b / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 10 passed -b without an installed band file"
fi
printf 'chr1\t0\t2300000\tp36.33\tgneg\n' > "${GENOME_DIR}/reference/cytoBand.hg38.txt"
: > "$FAKE_DOCKER_LOG"
run_expect 0 telomere "${SCRIPTS}/10-telomere-hunter.sh" sample1
docker_log_has '^run image=[^ ]*telomerehunter.* -b /genome/reference/cytoBand\.hg38\.txt' "step 10 did not pass the GRCh38 bands"
