#!/usr/bin/env bash
# A step whose tool fails half way must not leave a file the next run takes
# for a finished result, and must not report success on an older file.
#   Stranger (09b): the tool writes part of its VCF and exits 1. No VCF is
#     left behind, and the rerun annotates again instead of printing "already
#     exists". An empty VCF left by an older version is redone too.
#   Sniffles2 (04c): the tool writes part of its VCF and exits 1; the rerun
#     must replace it (Sniffles2 refuses an existing output without
#     --allow-overwrite).
#   Clair3 (03e): the tool exits 0 without merge_output.vcf.gz while an older
#     sample1.vcf.gz is present. The step fails instead of printing "complete".
#   FreeBayes (03b): bcftools sort writes part of its VCF and fails. The step
#     fails, keeps the raw VCF and leaves no sorted VCF, whole or truncated (a
#     pipe without pipefail hid the failure before, and the sort then wrote
#     straight to the final name).
# Each scenario failed against the scripts before this change.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/sample1/aligned_longread" "${GENOME_DIR}/sample1/expansion_hunter"
cp "${GENOME_DIR}/sample1/aligned/"* "${GENOME_DIR}/sample1/aligned_longread/"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' \
  > "${GENOME_DIR}/sample1/expansion_hunter/sample1_eh.vcf"

# The tool behaviour for this run: FAKE_MODE=fail or ok.
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=" $* "
header='##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n'
record='chr1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n'
# bgzf PATH: a gzip VCF ending with the BGZF end-of-file block, as htslib writes.
bgzf() {
  mkdir -p "$(dirname "$1")"
  { printf "${header}${record}" | gzip -c
    printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "$1"
}
# word_after NAME: the word after NAME in the command.
word_after() {
  local -a w
  read -r -a w <<<"$args"
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]}"; return 0; fi
  done
  return 1
}
case "$args" in
  *" stranger "*)
    if [ "${FAKE_MODE:-ok}" = fail ]; then printf '##fileformat=VCFv4.2\n'; exit 1; fi
    printf "${header}${record}" ;;
  *" sniffles "*)
    out=$(host_path "$(word_after -v)")
    if [ -e "$out" ] && [[ "$args" != *" --allow-overwrite "* ]]; then
      echo "fake Sniffles2: output ${out} exists; refusing to overwrite it" >&2
      exit 1
    fi
    mkdir -p "$(dirname "$out")"
    if [ "${FAKE_MODE:-ok}" = fail ]; then printf "${header}" | gzip -c | head -c 20 > "$out"; exit 1; fi
    bgzf "$out"; : > "${out}.tbi" ;;
  *"/opt/bin/run_clair3.sh "*)
    out=$(host_path "$(sed -n 's/.* --output=\([^ ]*\) .*/\1/p' <<<"$args")")
    mkdir -p "$out"
    if [ "${FAKE_MODE:-ok}" = ok ]; then bgzf "${out}/merge_output.vcf.gz"; : > "${out}/merge_output.vcf.gz.tbi"; fi ;;
  *" freebayes "*)
    printf "${header}${record}" ;;
  *"bcftools sort "*"| bcftools view"*)
    # The pipe of the script before this change: without pipefail the shell
    # reports the exit code of bcftools view, which succeeds.
    exit 0 ;;
  *" bcftools sort "*)
    out=$(host_path "$(word_after -o)")
    if [ "${FAKE_MODE:-ok}" = fail ]; then
      # As a killed or failing sort does: part of the output, then an error.
      mkdir -p "$(dirname "$out")"; printf "${header}" | gzip -c | head -c 20 > "$out"
      echo "fake bcftools sort: failed" >&2; exit 1
    fi
    bgzf "$out" ;;
  *" bcftools stats "*)
    printf 'SN\t0\tnumber of records:\t1\n' ;;
  *" bcftools index "*)
    p=$(host_path "${*: -1}") && : > "${p}.tbi" ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/tools-hook"
use_output_hook   # writes host-path.sh
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

# --- Stranger -------------------------------------------------------------------------
OUT="${GENOME_DIR}/sample1/expansion_hunter/sample1_eh_stranger.vcf"
FAKE_MODE=fail run_rc stranger-fail "${SCRIPTS}/09b-stranger.sh" sample1
[ "$RC" -ne 0 ] || fail "step 09b exited 0 although Stranger failed"
[ ! -e "$OUT" ] || fail "step 09b left ${OUT} behind after Stranger failed"
[ ! -e "${OUT}.tmp" ] || fail "step 09b left ${OUT}.tmp behind"
run_expect 0 stranger-rerun "${SCRIPTS}/09b-stranger.sh" sample1
output_lacks stranger-rerun 'already exists'
grep -q '^#CHROM' "$OUT" || fail "the rerun of step 09b wrote no VCF header"
docker_log_has 'variant_catalog_grch38\.json' "step 09b did not pass Stranger's GRCh38 catalog"
run_expect 0 stranger-again "${SCRIPTS}/09b-stranger.sh" sample1
output_has stranger-again 'already exists'
: > "$OUT"   # an empty file, as `stranger ... > OUT` left after a failed pull
run_expect 0 stranger-empty "${SCRIPTS}/09b-stranger.sh" sample1
output_lacks stranger-empty 'already exists'
[ -s "$OUT" ] || fail "step 09b kept an empty ${OUT} as its result"

# --- Sniffles2 ------------------------------------------------------------------------
SV="${GENOME_DIR}/sample1/sv_sniffles/sample1_sv.vcf.gz"
FAKE_MODE=fail run_rc sniffles-fail "${SCRIPTS}/04c-sniffles2.sh" sample1
[ "$RC" -ne 0 ] || fail "step 04c exited 0 although Sniffles2 failed"
run_expect 0 sniffles-rerun "${SCRIPTS}/04c-sniffles2.sh" sample1
have_bgzf_eof() { [ "$(tail -c 28 "$1" | od -An -v -tx1 | tr -d ' \n')" = 1f8b08040000000000ff0600424302001b0003000000000000000000 ]; }
have_bgzf_eof "$SV" || fail "the rerun of step 04c did not replace the half-written ${SV}"

# --- Clair3 ---------------------------------------------------------------------------
C3="${GENOME_DIR}/sample1/vcf_clair3"
mkdir -p "$C3"
printf 'an older result\n' | gzip -c > "${C3}/sample1.vcf.gz"
PLATFORM=hifi FAKE_MODE=fail run_rc clair3-nothing "${SCRIPTS}/03e-clair3.sh" sample1
[ "$RC" -ne 0 ] || fail "step 03e exited 0 although Clair3 wrote no merge_output.vcf.gz"
output_lacks clair3-nothing 'Clair3 complete'
PLATFORM=hifi run_expect 0 clair3-ok "${SCRIPTS}/03e-clair3.sh" sample1
have_bgzf_eof "${C3}/sample1.vcf.gz" || fail "step 03e did not put merge_output.vcf.gz in place of the older VCF"
[ ! -e "${C3}/merge_output.vcf.gz" ] || fail "step 03e left merge_output.vcf.gz beside its renamed copy"

# --- FreeBayes ------------------------------------------------------------------------
FB="${GENOME_DIR}/sample1/vcf_freebayes"
FAKE_MODE=fail run_rc freebayes-sort-fails "${SCRIPTS}/03b-freebayes.sh" sample1
[ "$RC" -ne 0 ] || fail "step 03b exited 0 although bcftools sort failed"
[ -s "${FB}/sample1_raw.vcf" ] || fail "step 03b deleted the raw VCF although sorting failed"
[ ! -e "${FB}/sample1.vcf.gz" ] || fail "step 03b left a sorted VCF although sorting failed"
[ ! -e "${FB}/sample1.vcf.gz.tmp" ] || fail "step 03b left the half-written sort output behind"
run_expect 0 freebayes-ok "${SCRIPTS}/03b-freebayes.sh" sample1
have_bgzf_eof "${FB}/sample1.vcf.gz" || fail "step 03b wrote no sorted VCF"
[ ! -e "${FB}/sample1_raw.vcf" ] || fail "step 03b kept the raw VCF after a good run"
