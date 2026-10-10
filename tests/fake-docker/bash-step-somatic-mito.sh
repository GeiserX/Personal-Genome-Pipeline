#!/usr/bin/env bash
# Step 29 (Mutect2 tumor-only) calls the CHIP driver genes by default, with
# the read orientation model and, when the common-sites VCF is installed, the
# contamination estimate. INTERVALS=genome calls the whole genome. A finished
# result is skipped only for the INTERVALS value and the optional resources it
# was called with.
# Step 20 (chrM) marks possible NuMTs with NuMTFilterTool at the median
# autosomal depth it reads from step 16b's mosdepth output. The Nextflow
# MITO_VARIANTS runs the same filter on the same mosdepth files: its script
# block, run here with a fake gatk, passes the same depth (30), and 0 when the
# workflow gives it no mosdepth files. The workflow joins MOSDEPTH's output
# in, and the module stops when mosdepth ran and the join gave it nothing
# (the stub runs in nextflow.yml go through that check).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n' > "${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.dict"
use_output_hook

# --- step 29 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 somatic "${SCRIPTS}/29-mutect2-somatic.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* Mutect2 .*--f1r2-tar-gz /genome/sample1/somatic/sample1_f1r2\.tar\.gz .*--intervals /genome/sample1/somatic/chip_genes_grch38\.bed' \
  "step 29 did not call the CHIP genes with --f1r2-tar-gz by default"
cmp -s "${REPO_ROOT}/assets/chip_genes_grch38.bed" "${GENOME_DIR}/sample1/somatic/chip_genes_grch38.bed" \
  || fail "step 29 did not use assets/chip_genes_grch38.bed"
docker_log_has '^run image=[^ ]*gatk.* LearnReadOrientationModel -I /genome/sample1/somatic/sample1_f1r2\.tar\.gz ' "step 29 did not learn the read orientation model"
docker_log_has '^run image=[^ ]*gatk.* FilterMutectCalls .*--ob-priors /genome/sample1/somatic/sample1_read-orientation-model\.tar\.gz ' \
  "FilterMutectCalls did not get the orientation priors"
output_has somatic 'Skipped: no common-sites VCF'
if grep -q 'GetPileupSummaries' "$FAKE_DOCKER_LOG"; then fail "step 29 ran GetPileupSummaries without a common-sites VCF"; fi

# Whole genome on request, with the contamination estimate once its VCF is there.
mkdir -p "${GENOME_DIR}/somatic"
printf 'placeholder\n' > "${GENOME_DIR}/somatic/small_exac_common_3.hg38.vcf.gz"
printf 'placeholder\n' > "${GENOME_DIR}/somatic/small_exac_common_3.hg38.vcf.gz.tbi"
rm -rf "${GENOME_DIR}/sample1/somatic"
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome run_expect 0 somatic-genome "${SCRIPTS}/29-mutect2-somatic.sh" sample1
if awk '/ Mutect2 / && /--intervals/ { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "INTERVALS=genome still passed --intervals"
fi
docker_log_has '^run image=[^ ]*gatk.* GetPileupSummaries .*-V /genome/somatic/small_exac_common_3\.hg38\.vcf\.gz ' "step 29 did not summarise pileups at the common sites"
docker_log_has '^run image=[^ ]*gatk.* CalculateContamination ' "step 29 did not estimate contamination"
docker_log_has '^run image=[^ ]*gatk.* FilterMutectCalls .*--contamination-table /genome/sample1/somatic/sample1_contamination\.table ' \
  "FilterMutectCalls did not get the contamination table"

# CHIP genes with the common sites: pileups only inside the genes.
rm -rf "${GENOME_DIR}/sample1/somatic"
: > "$FAKE_DOCKER_LOG"
run_expect 0 somatic-chip-contamination "${SCRIPTS}/29-mutect2-somatic.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* GetPileupSummaries .*-L /genome/sample1/somatic/chip_genes_grch38\.bed --interval-set-rule INTERSECTION ' \
  "GetPileupSummaries was not limited to the CHIP genes"

# A finished result is reused only for the intervals it was called on.
# bgzf_vcf PATH: a finished VCF, ending with the BGZF end-of-file block.
bgzf_vcf() {
  { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' | gzip -c
    printf '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'; } > "$1"
}
SO="${GENOME_DIR}/sample1/somatic"
grep -q '^INTERVALS=chip .* common_sites=yes$' "${SO}/sample1_somatic_filtered.run" 2>/dev/null \
  || fail "step 29 did not record INTERVALS=chip and the common sites beside its result"
bgzf_vcf "${SO}/sample1_somatic_filtered.vcf.gz"   # the CHIP run, finished
: > "$FAKE_DOCKER_LOG"
run_expect 0 somatic-chip-again "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_has somatic-chip-again 'Output already exists'
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome run_expect 0 somatic-genome-after-chip "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_lacks somatic-genome-after-chip 'Output already exists'
docker_log_has '^run image=[^ ]*gatk.* Mutect2 ' "INTERVALS=genome returned the CHIP result instead of calling the genome"
grep -q '^INTERVALS=genome ' "${SO}/sample1_somatic_filtered.run" 2>/dev/null \
  || fail "step 29 did not record INTERVALS=genome beside its result"
# A result filtered without the contamination estimate is redone once the
# common-sites VCF is there. Here the other way round: the VCF is removed.
bgzf_vcf "${SO}/sample1_somatic_filtered.vcf.gz"
mv "${GENOME_DIR}/somatic/small_exac_common_3.hg38.vcf.gz.tbi" "${CASE_WORK}/common.tbi"
: > "$FAKE_DOCKER_LOG"
INTERVALS=genome run_expect 0 somatic-resources-changed "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_lacks somatic-resources-changed 'Output already exists'
docker_log_has '^run image=[^ ]*gatk.* Mutect2 ' "step 29 reused a result filtered with the common sites after they were removed"
mv "${CASE_WORK}/common.tbi" "${GENOME_DIR}/somatic/small_exac_common_3.hg38.vcf.gz.tbi"
# A finished result with no record (an older version of this step) is called again.
bgzf_vcf "${SO}/sample1_somatic_filtered.vcf.gz"
rm -f "${SO}/sample1_somatic_filtered.run"
: > "$FAKE_DOCKER_LOG"
run_expect 0 somatic-no-record "${SCRIPTS}/29-mutect2-somatic.sh" sample1
output_lacks somatic-no-record 'Output already exists'
docker_log_has '^run image=[^ ]*gatk.* Mutect2 ' "step 29 reused a result without knowing its intervals"

# --- step 20 ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 mito-no-depth "${SCRIPTS}/20-mtoolbox.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* NuMTFilterTool .*--autosomal-coverage 0 ' "step 20 did not run NuMTFilterTool"
output_has mito-no-depth 'no mosdepth output'
M="${GENOME_DIR}/sample1/mosdepth"
mkdir -p "$M"
printf 'chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t1000\t30000\t30.00\t0\t60\nchr2\t1000\t30000\t30.00\t0\t60\nchrM\t16569\t1\t1000.00\t0\t2000\ntotal\t2000\t60000\t30.00\t0\t60\n' \
  > "${M}/sample1.mosdepth.summary.txt"
: > "${M}/sample1.mosdepth.global.dist.txt"
for c in chr1 chr2; do
  for d in 31 30 20 10 0; do
    case $d in 31) p=0.40 ;; 30) p=0.60 ;; 20) p=0.90 ;; 10) p=0.99 ;; 0) p=1.00 ;; esac
    printf '%s\t%s\t%s\n' "$c" "$d" "$p" >> "${M}/sample1.mosdepth.global.dist.txt"
  done
done
printf 'chrM\t2000\t1.00\ntotal\t31\t0.40\n' >> "${M}/sample1.mosdepth.global.dist.txt"
: > "$FAKE_DOCKER_LOG"
run_expect 0 mito "${SCRIPTS}/20-mtoolbox.sh" sample1
docker_log_has '^run image=[^ ]*gatk.* NuMTFilterTool .*--autosomal-coverage 30 ' "step 20 did not pass the median autosomal depth (30)"
docker_log_has '^run image=[^ ]*gatk.* NuMTFilterTool .*-O /genome/sample1/mito/sample1_chrM_filtered\.vcf\.gz' \
  "NuMTFilterTool does not write the file the reports read"

# --- MITO_VARIANTS (Nextflow) ------------------------------------------------------------
# The shell of the module's script block, with its Groovy values filled in:
# prefix, bam, reference and the mosdepth files ('' without mosdepth).
MOD="${REPO_ROOT}/modules/local/mito_variants/main.nf"
# shellcheck disable=SC2016  # the dollar signs are Groovy placeholders, matched literally
render_mito() {
  awk '/^    script:/ { on = 1; next } on && /^    """/ { n++; if (n == 2) exit; next } on && n == 1' "$MOD" \
    | sed -e 's/\${prefix}/sample1/g' -e 's/\${bam}/in.bam/g' -e 's/\${reference}/ref.fa/g' \
          -e "s|\\\${depth_files}|$1|g" -e '/^ *cat <<-END_VERSIONS/,$d' \
          -e 's/\\\$/$/g' -e 's/\\\\/\\/g'
}
NFW="${CASE_WORK}/mito-nf"
mkdir -p "${NFW}/bin"
cat > "${NFW}/bin/gatk" <<'FAKE'
#!/usr/bin/env bash
printf 'gatk %s\n' "$*" >> "${NFW_LOG:?}"
FAKE
chmod +x "${NFW}/bin/gatk"
cp "${M}/sample1.mosdepth.summary.txt" "${M}/sample1.mosdepth.global.dist.txt" "$NFW/"
# run_mito_module DEPTH_FILES: run the rendered block in $NFW; sets NUMT to the NuMTFilterTool call
run_mito_module() {
  render_mito "$1" > "${NFW}/script.sh"
  : > "${NFW}/gatk.log"
  (cd "$NFW" && PATH="${NFW}/bin:${PATH}" NFW_LOG="${NFW}/gatk.log" bash -euo pipefail script.sh) > "${NFW}/out.txt" 2>&1 \
    || { cat "${NFW}/script.sh" "${NFW}/out.txt"; fail "MITO_VARIANTS' script block failed with DEPTH_FILES='$1'"; }
  cat "${NFW}/out.txt"
  NUMT=$(grep ' NuMTFilterTool ' "${NFW}/gatk.log" || true)
}
run_mito_module "sample1.mosdepth.summary.txt sample1.mosdepth.global.dist.txt"
[ -n "$NUMT" ] || fail "MITO_VARIANTS does not run NuMTFilterTool"
grep -q -- '--autosomal-coverage 30 ' <<<"${NUMT} " || fail "MITO_VARIANTS did not pass the median autosomal depth (30), as step 20 does: ${NUMT}"
grep -q -- '-V sample1_chrM_mutect2_filtered\.vcf\.gz ' <<<"${NUMT} " || fail "NuMTFilterTool does not read FilterMutectCalls' output: ${NUMT}"
grep -q -- '-O sample1_chrM_filtered\.vcf\.gz$' <<<"$NUMT" || fail "NuMTFilterTool does not write the file MITO_HAPLOGROUP and the report read: ${NUMT}"
grep -q ' FilterMutectCalls .*-O sample1_chrM_mutect2_filtered\.vcf\.gz' "${NFW}/gatk.log" \
  || fail "FilterMutectCalls does not write the intermediate file NuMTFilterTool reads"
run_mito_module ""
grep -q -- '--autosomal-coverage 0 ' <<<"${NUMT} " || fail "MITO_VARIANTS without mosdepth files did not run the filter at depth 0: ${NUMT:-no call}"
grep -q 'mosdepth is not in --tools' "${NFW}/out.txt" || fail "MITO_VARIANTS without mosdepth files does not say why the depth is 0"
# The workflow joins MOSDEPTH's summary and distribution in, and the module
# refuses to go on without them when mosdepth ran.
grep -q 'ch_bam.join(ch_mosdepth_depth)' "${REPO_ROOT}/workflows/bam_analysis.nf" \
  || fail "workflows/bam_analysis.nf does not join MOSDEPTH's output into MITO_VARIANTS"
[ "$(grep -c "contains('mosdepth') && !mosdepth_summary" "$MOD")" -eq 2 ] \
  || fail "MITO_VARIANTS' script and stub do not both stop when mosdepth ran and its output is missing"
echo "MITO_VARIANTS marks possible NuMTs at the depth step 20 uses."
