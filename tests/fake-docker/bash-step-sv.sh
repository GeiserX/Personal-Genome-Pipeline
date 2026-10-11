#!/usr/bin/env bash
# The SV steps:
#   TIDDIT (04a): with only a BWA-MEM2 index next to the reference, assembly
#     is skipped and the log says so (TIDDIT's assembly calls classic bwa);
#     with the classic index it runs with assembly. When TIDDIT exits 0
#     without a VCF, as it does when one of its own checks fails, the step
#     fails and prints the end of TIDDIT's log.
#   GRIDSS (04b): runs with --workingdir sv_gridss/work, removed once the VCF
#     is written; it is skipped, with the reason, when Docker has less memory
#     than GRIDSS_MIN_MEM_GB (default 32; the fake Docker reports 16 GiB).
#   Delly (19): gets -x with the exclude map setup.sh installs; without it,
#     runs as before and warns.
#   duphold (15): writes duphold/<sample>_sv_filtered.vcf.gz with the filter
#     expression of the Nextflow DUPHOLD_FILTER, and its index.
#   AnnotSV (05): annotates that file, the input of the Nextflow ANNOTSV;
#     without it, or when Manta's calls are newer, Manta's calls with a note.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
REF="${GENOME_DIR}/reference/GRCh38_no_alt_analysis_set.fasta"
use_output_hook
cat > "${CASE_WORK}/sv-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args=" ${*:2} "
case "$args" in
  *" tiddit "*)
    if [ "${FAKE_TIDDIT_NO_VCF:-}" = 1 ]; then
      echo "fake TIDDIT: the BAM index does not match the BAM; stopping"
      exit 0
    fi
    o=$(sed -n 's/.* -o \([^ ]*\) .*/\1/p' <<<"$args")
    h=$(host_path "$o")
    mkdir -p "$(dirname "$h")"
    printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n' > "${h}.vcf" ;;
  *" gridss "*|*"gridss.sh "*|*"/gridss "*)
    w=$(sed -n 's/.* --workingdir \([^ ]*\) .*/\1/p' <<<"$args")
    if [ -n "$w" ]; then mkdir -p "$(host_path "$w")"; : > "$(host_path "$w")/intermediate.bam"; fi ;;
  *" bcftools stats "*)
    printf 'SN\t0\tnumber of records:\t1\n' ;;
  *" bcftools view -e "*)
    # The expression and the files of a filter call, one per line
    printf '%s\n' "${5}" "${*: -1}" > "${CASE_WORK}/bcftools-filter" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/sv-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/sv-hook"

# --- TIDDIT ---------------------------------------------------------------------------
: > "${REF}.bwt.2bit.64"
: > "$FAKE_DOCKER_LOG"
run_expect 0 tiddit-mem2 "${SCRIPTS}/04a-tiddit.sh" sample1
output_has tiddit-mem2 'assembly skipped'
docker_log_has '^run image=[^ ]*tiddit.* --skip_assembly ' "TIDDIT ran without --skip_assembly although only a BWA-MEM2 index exists"
for ext in amb ann bwt pac sa; do : > "${REF}.${ext}"; done
: > "$FAKE_DOCKER_LOG"
run_expect 0 tiddit-bwa "${SCRIPTS}/04a-tiddit.sh" sample1
output_has tiddit-bwa 'local assembly enabled'
if grep -q -- '--skip_assembly' "$FAKE_DOCKER_LOG"; then fail "TIDDIT skipped assembly although the classic BWA index exists"; fi
FAKE_TIDDIT_NO_VCF=1 run_rc tiddit-no-vcf "${SCRIPTS}/04a-tiddit.sh" sample1
[ "$RC" -ne 0 ] || fail "step 04a exited 0 although TIDDIT wrote no VCF"
output_has tiddit-no-vcf 'the BAM index does not match'

# --- GRIDSS ---------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 gridss-low-mem "${SCRIPTS}/04b-gridss.sh" sample1
output_has gridss-low-mem 'SKIPPED: Docker has [0-9]+ GB of memory; GRIDSS needs 32 GB'
if awk '/^run / && /gridss/ { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "GRIDSS was started although Docker has less memory than it needs"
fi
GRIDSS_MIN_MEM_GB=8 run_expect 0 gridss "${SCRIPTS}/04b-gridss.sh" sample1
docker_log_has '^run image=[^ ]*gridss.* --workingdir /genome/sample1/sv_gridss/work ' "GRIDSS ran without --workingdir sv_gridss/work"
docker_log_has "^run image=[^ ]*gridss.* -t ${THREADS:-8} " "GRIDSS did not get THREADS"
[ ! -e "${GENOME_DIR}/sample1/sv_gridss/work" ] || fail "step 04b kept sv_gridss/work after GRIDSS wrote its VCF"

# --- Delly ----------------------------------------------------------------------------
: > "$FAKE_DOCKER_LOG"
run_expect 0 delly-no-map "${SCRIPTS}/19-delly.sh" sample1
output_has delly-no-map "WARNING: Delly's exclude map is not installed"
if awk '/^run / && /delly sr/ && / -x / { found = 1 } END { exit !found }' "$FAKE_DOCKER_LOG"; then
  fail "Delly got -x although no exclude map is installed"
fi
printf 'chr1\t0\t10000\ttelomere\n' > "${GENOME_DIR}/reference/delly_human.hg38.excl.tsv"
: > "$FAKE_DOCKER_LOG"
run_expect 0 delly "${SCRIPTS}/19-delly.sh" sample1
docker_log_has '^run image=[^ ]*delly.* delly sr .*-x /genome/reference/delly_human\.hg38\.excl\.tsv ' "Delly ran without its exclude map"

# --- duphold --------------------------------------------------------------------------
mkdir -p "${GENOME_DIR}/sample1/manta/results/variants"
printf 'placeholder\n' > "${GENOME_DIR}/sample1/manta/results/variants/diploidSV.vcf.gz"
run_expect 0 duphold "${SCRIPTS}/15-duphold.sh" sample1
D="${GENOME_DIR}/sample1/duphold"
[ -f "${D}/sample1_sv_filtered.vcf.gz" ] && [ -f "${D}/sample1_sv_filtered.vcf.gz.tbi" ] \
  || fail "step 15 did not write duphold/sample1_sv_filtered.vcf.gz with its index"
# The same expression as DUPHOLD_FILTER, on the annotated VCF
WANT=$(awk -F"'" '/^ *-e .*SVTYPE/ {print $2; exit}' "${REPO_ROOT}/modules/local/duphold/main.nf")
[ -n "$WANT" ] || fail "no filter expression found in modules/local/duphold/main.nf"
[ "$(sed -n 1p "${CASE_WORK}/bcftools-filter" 2>/dev/null)" = "$WANT" ] \
  || fail "step 15 did not filter with DUPHOLD_FILTER's expression (${WANT}): $(sed -n 1p "${CASE_WORK}/bcftools-filter" 2>/dev/null)"
[ "$(sed -n 2p "${CASE_WORK}/bcftools-filter" 2>/dev/null)" = /genome/sample1/duphold/sample1_sv_duphold.vcf ] \
  || fail "step 15 did not filter duphold's annotated VCF"
output_has duphold 'Filtered: .*sample1_sv_filtered\.vcf\.gz'

# --- AnnotSV --------------------------------------------------------------------------
mkdir -p "${GENOME_DIR}/annotsv_annotations/Annotations_Human/Genes/GRCh38"
: > "$FAKE_DOCKER_LOG"
run_expect 0 annotsv "${SCRIPTS}/05-annotsv.sh" sample1
docker_log_has '^run image=[^ ]*annotsv.* -SVinputFile /genome/sample1/duphold/sample1_sv_filtered\.vcf\.gz ' \
  "step 05 did not annotate step 15's filtered calls"
output_lacks annotsv 'NOTE:'
# Manta's calls newer than the filtered ones: Manta's, with a note.
touch -t 209901010000 "${GENOME_DIR}/sample1/manta/results/variants/diploidSV.vcf.gz"
: > "$FAKE_DOCKER_LOG"
run_expect 0 annotsv-stale "${SCRIPTS}/05-annotsv.sh" sample1
docker_log_has '^run image=[^ ]*annotsv.* -SVinputFile /genome/sample1/manta/results/variants/diploidSV\.vcf\.gz ' \
  "step 05 annotated filtered calls older than Manta's"
output_has annotsv-stale 'NOTE: .*older than Manta'
# No filtered calls at all: Manta's, with a note.
rm -f "${D}/sample1_sv_filtered.vcf.gz" "${D}/sample1_sv_filtered.vcf.gz.tbi"
: > "$FAKE_DOCKER_LOG"
run_expect 0 annotsv-no-duphold "${SCRIPTS}/05-annotsv.sh" sample1
docker_log_has '^run image=[^ ]*annotsv.* -SVinputFile /genome/sample1/manta/results/variants/diploidSV\.vcf\.gz ' \
  "step 05 did not fall back to Manta's calls without step 15"
output_has annotsv-no-duphold 'NOTE: no duphold-filtered calls'
# SV_VCF still wins.
mkdir -p "${GENOME_DIR}/sample1/sniffles"
printf 'placeholder\n' > "${GENOME_DIR}/sample1/sniffles/sv.vcf.gz"
: > "$FAKE_DOCKER_LOG"
SV_VCF="${GENOME_DIR}/sample1/sniffles/sv.vcf.gz" run_expect 0 annotsv-sv-vcf "${SCRIPTS}/05-annotsv.sh" sample1
docker_log_has '^run image=[^ ]*annotsv.* -SVinputFile /genome/sample1/sniffles/sv\.vcf\.gz ' "SV_VCF did not pick step 05's input"
