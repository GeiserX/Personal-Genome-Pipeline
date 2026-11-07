#!/usr/bin/env bash
# Two steps no e2e case runs, against the fake docker:
#   Strelka2 (03c): configures the workflow on the sample's BAM (or ALIGN_DIR's)
#     into vcf_strelka2/, runs it, counts the records; a finished run is not
#     repeated, an unfinished one is resumed without configuring again, a
#     run that writes no VCF fails, and a missing BAM index stops it first.
#   CNVpytor (18): asks the image for its data folder and mounts the pinned
#     resource files there read-only, imports the canonical chromosomes only,
#     runs -his, -partition and -call at 1 kb, feeds -view its commands on
#     stdin (docker -i), and reheaders the VCF from the reference's .fai. It
#     stops without the hg38 files, without one of the other five, and when
#     the image's data folder cannot be found.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook
cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
case "$args" in
  *"/runWorkflow.py "*)
    # Strelka2's workflow writes its VCF and index under the runDir.
    if [ "${FAKE_STRELKA_NO_VCF:-}" != 1 ]; then
      d=$(host_path /genome/sample1/vcf_strelka2/results/variants)
      mkdir -p "$d"
      printf 'placeholder\n' > "${d}/variants.vcf.gz"
      printf 'placeholder\n' > "${d}/variants.vcf.gz.tbi"
    fi ;;
  *" configureStrelkaGermlineWorkflow.py "*)
    d=$(host_path /genome/sample1/vcf_strelka2)
    mkdir -p "$d" && printf '#!/usr/bin/env python\n' > "${d}/runWorkflow.py" ;;
  *" bcftools stats "*)
    printf 'SN\t0\tnumber of samples:\t1\nSN\t0\tnumber of records:\t4242\n' ;;
  *" import cnvpytor, os"*)
    [ -n "${FAKE_CNVPYTOR_NO_DATA:-}" ] || echo /usr/local/lib/python3.12/site-packages/cnvpytor/data ;;
  *" -view 1000 "*)
    cat > "${CASE_WORK}/view-stdin.txt" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook"

# --- step 03c ---------------------------------------------------------------------------
SD="${GENOME_DIR}/sample1/vcf_strelka2"
mv "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai" "${CASE_WORK}/bai"
: > "$FAKE_DOCKER_LOG"
run_expect 1 strelka-no-index "${SCRIPTS}/03c-strelka2-germline.sh" sample1
output_has strelka-no-index 'ERROR: File not found: .*/aligned/sample1_sorted\.bam\.bai'
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "step 03c started a container without the BAM index"; fi
mv "${CASE_WORK}/bai" "${GENOME_DIR}/sample1/aligned/sample1_sorted.bam.bai"

: > "$FAKE_DOCKER_LOG"
run_expect 0 strelka "${SCRIPTS}/03c-strelka2-germline.sh" sample1
docker_log_has '^run image=[^ ]*strelka.* configureStrelkaGermlineWorkflow\.py --bam /genome/sample1/aligned/sample1_sorted\.bam --referenceFasta /genome/reference/GRCh38_no_alt_analysis_set\.fasta --runDir /genome/sample1/vcf_strelka2( |$)' \
  "step 03c did not configure Strelka2 on the sample's BAM into vcf_strelka2/"
docker_log_has '^run image=[^ ]*strelka.* /genome/sample1/vcf_strelka2/runWorkflow\.py -m local -j [0-9]+' \
  "step 03c did not run the configured workflow locally"
docker_log_has '^run image=[^ ]*bcftools.* bcftools stats /genome/sample1/vcf_strelka2/results/variants/variants\.vcf\.gz' \
  "step 03c did not count the records of Strelka2's VCF"
output_has strelka 'Total variants called: 4242$'

: > "$FAKE_DOCKER_LOG"
run_expect 0 strelka-again "${SCRIPTS}/03c-strelka2-germline.sh" sample1
output_has strelka-again 'Strelka2 already done for sample1'
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "step 03c ran a container again over a finished result"; fi

# An interrupted run (the workflow, no results) is resumed, not configured again.
rm -rf "${SD}/results"
: > "$FAKE_DOCKER_LOG"
run_expect 0 strelka-resume "${SCRIPTS}/03c-strelka2-germline.sh" sample1
output_has strelka-resume 'Found an unfinished Strelka2 run'
if grep -q 'configureStrelkaGermlineWorkflow' "$FAKE_DOCKER_LOG"; then fail "step 03c configured again over an unfinished run"; fi
docker_log_has '^run image=[^ ]*strelka.* /genome/sample1/vcf_strelka2/runWorkflow\.py ' "step 03c did not resume the unfinished run"

# A workflow that ends without a VCF fails the step.
rm -rf "${SD}/results"
: > "$FAKE_DOCKER_LOG"
FAKE_STRELKA_NO_VCF=1 run_expect 1 strelka-no-vcf "${SCRIPTS}/03c-strelka2-germline.sh" sample1
output_has strelka-no-vcf 'ERROR: Strelka2 finished without .*/vcf_strelka2/results/variants/variants\.vcf\.gz'

# ALIGN_DIR picks the BWA-MEM2 alignment; a leftover folder with no workflow
# is removed first.
rm -rf "$SD"
mkdir -p "${SD}/junk" "${GENOME_DIR}/sample1/aligned_bwamem2"
printf 'placeholder\n' > "${GENOME_DIR}/sample1/aligned_bwamem2/sample1_sorted.bam"
printf 'placeholder\n' > "${GENOME_DIR}/sample1/aligned_bwamem2/sample1_sorted.bam.bai"
: > "$FAKE_DOCKER_LOG"
ALIGN_DIR=aligned_bwamem2 run_expect 0 strelka-bwamem2 "${SCRIPTS}/03c-strelka2-germline.sh" sample1
docker_log_has '^run image=[^ ]*strelka.* rm -rf /genome/sample1/vcf_strelka2' "step 03c did not remove the leftover folder"
docker_log_has '^run image=[^ ]*strelka.* configureStrelkaGermlineWorkflow\.py --bam /genome/sample1/aligned_bwamem2/sample1_sorted\.bam ' \
  "ALIGN_DIR=aligned_bwamem2 did not reach Strelka2"

# --- step 18 ----------------------------------------------------------------------------
CD="${GENOME_DIR}/reference/cnvpytor"
mkdir -p "$CD"
: > "$FAKE_DOCKER_LOG"
run_expect 1 cnvpytor-no-data "${SCRIPTS}/18-cnvpytor.sh" sample1
output_has cnvpytor-no-data "CNVpytor resource 'gc_hg38.pytor' missing or empty"
printf 'x\n' > "${CD}/gc_hg38.pytor"
printf 'x\n' > "${CD}/mask_hg38.pytor"
for f in gc_hg19.pytor mask_hg19.pytor gc_chm13v2.0.pytor gc_chm13v1.1.pytor; do : > "${CD}/${f}"; done
run_expect 1 cnvpytor-one-missing "${SCRIPTS}/18-cnvpytor.sh" sample1
output_has cnvpytor-one-missing "CNVpytor resource 'gc_kn99.pytor' missing .*all 7 resource files must exist"
if grep -q '^run ' "$FAKE_DOCKER_LOG"; then fail "step 18 started a container without its resource files"; fi
: > "${CD}/gc_kn99.pytor"

FAKE_CNVPYTOR_NO_DATA=1 run_expect 1 cnvpytor-no-image-data "${SCRIPTS}/18-cnvpytor.sh" sample1
output_has cnvpytor-no-image-data 'could not find the cnvpytor data directory'

: > "$FAKE_DOCKER_LOG"
run_expect 0 cnvpytor "${SCRIPTS}/18-cnvpytor.sh" sample1
MNT="${CD}:/usr/local/lib/python3\\.12/site-packages/cnvpytor/data:ro"
PYTOR=/genome/sample1/cnvpytor/sample1\\.pytor
docker_log_has "^run image=[^ ]*cnvpytor.*-v ${MNT} .* cnvpytor -root ${PYTOR} -rd /genome/sample1/aligned/sample1_sorted\\.bam -chrom chr1 chr2 .* chr22 chrX chrY -j 4( |$)" \
  "step 18 did not import the canonical chromosomes, with the resource files on the image's data folder"
if grep -E ' -rd .* -chrom .*(chrM|chrUn|_random|chrEBV)' "$FAKE_DOCKER_LOG"; then fail "step 18 imported a contig the GC data does not cover"; fi
for a in '-his 1000' '-partition 1000' '-call 1000'; do
  docker_log_has "^run image=[^ ]*cnvpytor.*-v ${MNT} .* cnvpytor -root ${PYTOR} ${a}( |$)" "step 18 did not run cnvpytor ${a}"
done
docker_log_has "^run image=[^ ]*cnvpytor.* -i .*-v ${MNT} .* cnvpytor -root ${PYTOR} -view 1000( |$)" "step 18 did not run -view with docker -i"
grep -qx 'set print_filename /genome/sample1/cnvpytor/sample1_cnvs.raw.vcf' "${CASE_WORK}/view-stdin.txt" 2>/dev/null \
  || fail "cnvpytor -view did not get the file name to print on its stdin"
grep -qx 'print calls' "${CASE_WORK}/view-stdin.txt" || fail "cnvpytor -view was not told to print the calls"
docker_log_has '^run image=[^ ]*bcftools.*bcftools reheader --fai /genome/reference/GRCh38_no_alt_analysis_set\.fasta\.fai ' \
  "step 18 did not reheader the VCF from the reference's .fai"
test -f "${GENOME_DIR}/sample1/cnvpytor/sample1_cnvs.txt" || fail "step 18 wrote no sample1_cnvs.txt"
output_has cnvpytor 'Calls \(VCF\): .*/sample1/cnvpytor/sample1_cnvs\.vcf\.gz'
echo "PASS steps 03c and 18 against the fake docker"
