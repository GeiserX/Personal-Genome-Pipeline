#!/usr/bin/env bash
# Step 13 (VEP) asks for the same annotation as the Nextflow module and leaves
# no stale file for the steps after it:
#   - vep gets --fasta (without it VEP turns HGVS off), --compress_output
#     bgzip, --cache_version of versions.env and --fork THREADS;
#   - the output is <sample>_vep.vcf.gz with its .tbi;
#   - a finished run removes vcfanno's <sample>_annotated.vcf.gz (built from
#     the previous annotation, and preferred by steps 30, 23 and 31) and the
#     uncompressed <sample>_vep.vcf of earlier versions;
#   - a failed VEP run keeps the previous result and the derived files;
#   - an input VCF without its index is refused before VEP starts;
#   - vep reads the PASS and '.' records only: a copy bcftools view -f PASS,.
#     writes, removed afterwards. The Nextflow VEP module's awk filter keeps
#     the same records of a VCF with PASS, '.', RefCall and LowQual ones.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38"
echo "species homo_sapiens" > "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt"
V="${GENOME_DIR}/sample1/vep"
mkdir -p "$V"
for f in sample1_annotated.vcf.gz sample1_annotated.vcf.gz.tbi sample1_vep.vcf; do
  printf 'from the previous VEP run\n' > "${V}/${f}"
done
use_output_hook
# FAKE_VEP_FAIL=1 makes vep exit 1; everything else goes to the generic hook.
cat > "${CASE_WORK}/vep-hook" <<'HOOK'
#!/usr/bin/env bash
if [ "${FAKE_VEP_FAIL:-}" = 1 ] && [ "${2:-}" = vep ]; then echo "fake vep: failed" >&2; exit 1; fi
# Record what vep reads, and the bcftools view -f calls (FILTER, output, input).
args=("$@")
for i in "${!args[@]}"; do
  [ "${args[$i]}" = --input_file ] && [ "${2:-}" = vep ] && printf '%s\n' "${args[$((i + 1))]}" > "${CASE_WORK}/vep-input"
done
if [ "${2:-} ${3:-} ${4:-}" = "bcftools view -f" ]; then
  printf '%s %s %s\n' "$5" "${8:-}" "${9:-}" > "${CASE_WORK}/pass-filter"
fi
exec "${CASE_WORK:?}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/vep-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/vep-hook"

THREADS=3 run_expect 0 vep "${SCRIPTS}/13-vep-annotation.sh" sample1
docker_log_has '^run image=[^ ]*ensembl-vep.* vep .*--fasta /genome/reference/GRCh38_no_alt_analysis_set\.fasta ' \
  "step 13 did not pass the reference FASTA to vep"
docker_log_has '^run image=[^ ]*ensembl-vep.* --compress_output bgzip ' "step 13 did not ask vep for bgzip output"
docker_log_has "^run image=[^ ]*ensembl-vep.* --cache_version ${VEP_CACHE_RELEASE} " "step 13 did not pass --cache_version ${VEP_CACHE_RELEASE}"
docker_log_has '^run image=[^ ]*ensembl-vep.* --cpus 3 .* --fork 3( |$)' "step 13 did not use THREADS=3 for --cpus and --fork"
[ -s "${V}/sample1_vep.vcf.gz" ] || fail "step 13 wrote no sample1_vep.vcf.gz"
[ -f "${V}/sample1_vep.vcf.gz.tbi" ] || fail "step 13 did not index sample1_vep.vcf.gz"
for f in sample1_annotated.vcf.gz sample1_annotated.vcf.gz.tbi sample1_vep.vcf; do
  [ ! -e "${V}/${f}" ] || fail "step 13 left ${f} from the previous annotation in place"
done
# VEP reads the PASS and '.' records, through a temporary copy.
[ "$(cat "${CASE_WORK}/pass-filter" 2>/dev/null)" = "PASS,. /genome/sample1/vep/sample1.pass.tmp.vcf.gz /genome/sample1/vcf/sample1.vcf.gz" ] \
  || fail "step 13 did not write the PASS,. records of the sample VCF to vep/sample1.pass.tmp.vcf.gz: $(cat "${CASE_WORK}/pass-filter" 2>/dev/null || echo 'no bcftools view -f call')"
[ "$(cat "${CASE_WORK}/vep-input" 2>/dev/null)" = /genome/sample1/vep/sample1.pass.tmp.vcf.gz ] \
  || fail "vep did not read the PASS records only: --input_file $(cat "${CASE_WORK}/vep-input" 2>/dev/null)"
[ ! -e "${V}/sample1.pass.tmp.vcf.gz" ] || fail "step 13 left its temporary PASS copy behind"

# A failed rerun keeps the previous result and the files built from it.
printf 'built from the result above\n' > "${V}/sample1_annotated.vcf.gz"
BEFORE=$(cksum < "${V}/sample1_vep.vcf.gz")
FAKE_VEP_FAIL=1 run_rc vep-fails "${SCRIPTS}/13-vep-annotation.sh" sample1
[ "$RC" -ne 0 ] || fail "step 13 exited 0 although vep failed"
[ "$(cksum < "${V}/sample1_vep.vcf.gz")" = "$BEFORE" ] || fail "a failed VEP run changed the previous sample1_vep.vcf.gz"
[ -e "${V}/sample1_annotated.vcf.gz" ] || fail "a failed VEP run removed the files built from the previous result"

# No index, no annotation: the VCF may be half written.
rm "${GENOME_DIR}/sample1/vcf/sample1.vcf.gz.tbi"
: > "$FAKE_DOCKER_LOG"
run_rc vep-no-index "${SCRIPTS}/13-vep-annotation.sh" sample1
[ "$RC" -ne 0 ] || fail "step 13 annotated a VCF without its index"
output_has vep-no-index 'sample1\.vcf\.gz\.tbi'
if awk '/^run / && /ensembl-vep/ { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "step 13 started VEP on a VCF without its index"
fi

# --- the Nextflow VEP module: the same records, selected with awk -------------------------
# The module's filter line, with the Groovy escapes undone and its own files
# replaced: bgzip -dc <vcf> | awk ... | bgzip -c > <out>
MOD="${REPO_ROOT}/modules/local/vep/main.nf"
LINE=$(awk '/^process VEP \{/ { on = 1 } on && /bgzip -dc .* \| awk / { print; exit }' "$MOD")
[ -n "$LINE" ] || fail "the VEP module selects no records before vep (no 'bgzip -dc ... | awk' line)"
FILTER=${LINE#*| }
FILTER=${FILTER% | bgzip -c*}
FILTER=$(sed -e 's/\\\$/$/g' -e 's/\\\\/\\/g' <<<"$FILTER")
echo "VEP module filter: ${FILTER}"
T="${CASE_WORK}/vep-module"
mkdir -p "$T"
printf '##fileformat=VCFv4.2\n##FILTER=<ID=RefCall,Description="x">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n' > "${T}/in.vcf"
printf 'chr1\t%s\t.\tA\tG\t50\t%s\t.\tGT\t0/1\n' 1 PASS 2 RefCall 3 . 4 LowQual 5 'PASS;x' >> "${T}/in.vcf"
gzip -c "${T}/in.vcf" > "${T}/in.vcf.gz"
gzip -dc "${T}/in.vcf.gz" | bash -c "$FILTER" > "${T}/out.vcf"
cat "${T}/out.vcf"
[ "$(grep -c '^#' "${T}/out.vcf")" -eq 3 ] || fail "the VEP module's filter dropped header lines"
# bcftools view -f PASS,. keeps a record whose FILTER list holds PASS (5).
[ "$(grep -v '^#' "${T}/out.vcf" | cut -f2 | tr '\n' ' ')" = "1 3 5 " ] \
  || fail "the VEP module's filter did not keep exactly the PASS and '.' records (kept positions: $(grep -v '^#' "${T}/out.vcf" | cut -f2 | tr '\n' ' '))"
