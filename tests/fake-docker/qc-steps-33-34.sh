#!/usr/bin/env bash
# Steps 33 (sample identity and contamination) and 34 (CRAM archive), with a
# hook that plays somalier, VerifyBamID2 and samtools, and runs the verdict
# (bin/collect_summary.py sample-qc) for real on the host.
#
# Step 33: a sex somalier infers that differs from the declared one stops the
# step; SEX_CHECK=warn goes on; a sex somalier cannot tell, or reads as X and Y
# disagreeing (-2), is not checked; somalier gets the declared sex as a
# one-line pedigree (--ped), and the pedigree's sex it keeps when chrX cannot
# tell is not checked either;
# FREEMIX above 0.03 warns and exits 0; "Insufficient Available markers"
# reruns VerifyBamID2 with --DisableSanityCheck and records it; any other
# VerifyBamID2 failure, or missing data, stops the step.
# Step 34: a CRAM whose flagstat differs from the BAM's, or that fails
# quickcheck, is removed and the BAM kept, also with --delete-bam; the BAM is
# deleted only after the check passed; --restore refuses to write over a BAM
# and leaves no partial BAM when it fails; archive and restore both stop
# while another run holds the sample's lock, and the lock goes with its run.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
use_output_hook   # writes host-path.sh, which the hook below reads
mkdir -p "${GENOME_DIR}/reference/somalier" "${GENOME_DIR}/reference/verifybamid2"
printf 'placeholder\n' > "${GENOME_DIR}/reference/somalier/sites.hg38.vcf.gz"
for e in UD mu bed; do
  printf 'placeholder\n' > "${GENOME_DIR}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.${e}"
done

cat > "${CASE_WORK}/hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
image=$1; shift
arg() {  # arg FLAG: the value after FLAG in the call
  local prev=""
  for a in "$@"; do [ "$prev" = "$FLAG" ] && { printf '%s' "$a"; return; }; prev=$a; done
}
case "$*" in
  "somalier extract"*)
    FLAG=-d; d=$(arg "$@"); h=$(host_path "$d"); mkdir -p "$h"
    printf 'fake' > "${h}/${FAKE_SM:-sample1}.somalier" ;;
  "somalier relate"*)
    FLAG=-o; o=$(host_path "$(arg "$@")")
    # As somalier: the pedigree's sex is original_pedigree_sex, and the sex
    # column starts from it (FAKE_SOMALIER_SEX plays what the reads change it to).
    FLAG=--ped; ped=$(arg "$@"); ped_code=-9; ped_sex=-9
    if [ -n "$ped" ]; then
      ped_code=$(awk -F'\t' '{print $5; exit}' "$(host_path "$ped")")
      case "$ped_code" in 1) ped_sex=male ;; 2) ped_sex=female ;; *) ped_sex=unknown ;; esac
    fi
    printf '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\toriginal_pedigree_sex\tgt_depth_mean\tn_hom_ref\tn_het\tn_hom_alt\tX_depth_mean\tX_n\tX_hom_ref\tX_het\tX_hom_alt\tY_depth_mean\tY_n\n%s\t%s\t-9\t-9\t%s\t-9\t%s\t30.0\t5000\t6000\t4000\t15.0\t%s\t40\t%s\t260\t14.0\t17\n' \
      "${FAKE_SM:-sample1}" "${FAKE_SM:-sample1}" "${FAKE_SOMALIER_SEX:-$ped_code}" "$ped_sex" \
      "${FAKE_X_N:-300}" "${FAKE_X_HET:-0}" > "${o}.samples.tsv"
    printf '#sample_a\tsample_b\trelatedness\n' > "${o}.pairs.tsv" ;;
  "verifybamid2 "*)
    FLAG=--Output; o=$(host_path "$(arg "$@")")
    case " $* " in *" --DisableSanityCheck "*) few=false ;; *) few=${FAKE_VB2_FEW:-false} ;; esac
    if [ -n "${FAKE_VB2_FAIL:-}" ]; then echo "Segmentation fault"; exit 139; fi
    if $few; then echo "WARNING - Insufficient Available markers, check input bam depth"; exit 1; fi
    printf '#SEQ_ID\tRG\tCHIP_ID\t#SNPS\t#READS\tAVG_DP\tFREEMIX\n%s\tNA\tNA\t%s\t1\t30\t%s\n' \
      "${FAKE_SM:-sample1}" "$([ "${FAKE_VB2_FEW:-false}" = true ] && echo 320 || echo 99000)" "${FAKE_FREEMIX:-0.004}" > "${o}.selfSM" ;;
  "python3 /pgp-bin/"*)
    # The real verdict, with the container paths mapped to the host.
    args=()
    for a in "$@"; do
      case "$a" in
        /pgp-bin/*) args+=("${REPO_ROOT}/bin/${a#/pgp-bin/}") ;;
        /genome*) args+=("$(host_path "$a")") ;;
        *) args+=("$a") ;;
      esac
    done
    exec "${args[@]}" ;;
  "samtools view"*" -C "*)
    FLAG=-o; h=$(host_path "$(arg "$@")"); printf 'fake cram' > "$h" ;;
  "samtools view"*" -b "*)
    FLAG=-o; h=$(host_path "$(arg "$@")"); printf 'fake bam' > "$h"
    if [ -n "${FAKE_VIEW_FAIL:-}" ]; then echo "samtools view: disk full" >&2; exit 1; fi ;;
  "samtools index"*)
    for a in "$@"; do last=$a; done
    h=$(host_path "$last"); case "$h" in *.cram) : > "${h}.crai" ;; *) : > "${h}.bai" ;; esac ;;
  "samtools quickcheck"*)
    for a in "$@"; do last=$a; done
    if [ "${FAKE_QUICKCHECK_FAIL:-}" = "${last##*.}" ]; then echo "${last} had no EOF block" >&2; exit 1; fi ;;
  "samtools flagstat"*)
    for a in "$@"; do last=$a; done
    n=1000
    case "$last" in *.cram) n=${FAKE_CRAM_READS:-1000} ;; esac
    printf '%s + 0 in total (QC-passed reads + QC-failed reads)\n%s + 0 primary\n' "$n" "$n" ;;
esac
HOOK
chmod +x "${CASE_WORK}/hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook"
T="${GENOME_DIR}/sample1/qc/sample1_sample_qc.tsv"
tval() { awk -F'\t' -v k="$1" '$1 == k { print $2; exit }' "$T"; }

# --- step 33 -----------------------------------------------------------------------------
FAKE_SOMALIER_SEX=1 run_expect 1 qc-mismatch "${SCRIPTS}/33-sample-qc.sh" sample1 female
output_has qc-mismatch 'SEX CHECK MISMATCH: declared female, somalier infers male'
output_has qc-mismatch 'chrX sites 300: 0 heterozygous, 260 homozygous ALT'
[ "$(tval sex_check)" = mismatch ] || fail "the table does not record the mismatch"

FAKE_SOMALIER_SEX=1 SEX_CHECK=warn run_expect 0 qc-mismatch-warn "${SCRIPTS}/33-sample-qc.sh" sample1 female
output_has qc-mismatch-warn 'SEX_CHECK=warn: continuing'

FAKE_SOMALIER_SEX=1 run_expect 0 qc-match "${SCRIPTS}/33-sample-qc.sh" sample1 male
output_has qc-match 'Sex check: OK'
docker_log_has 'somalier relate --infer --ped /genome/sample1/qc/somalier/sample1\.declared\.ped ' "somalier relate got no pedigree"
[ "$(cat "${GENOME_DIR}/sample1/qc/somalier/sample1.declared.ped")" = "$(printf 'sample1\tsample1\t-9\t-9\t1\t-9')" ] \
  || fail "the pedigree does not hold the declared male: $(cat "${GENOME_DIR}/sample1/qc/somalier/sample1.declared.ped")"
[ "$(tval contamination)" = ok ] || fail "FREEMIX 0.004 was not reported ok"
[ "$(tval verifybamid2_marker_check)" = passed ] || fail "the marker check was not recorded as passed"

FAKE_SOMALIER_SEX=-9 run_expect 0 qc-unknown "${SCRIPTS}/33-sample-qc.sh" sample1 female
output_has qc-unknown 'Sex check: not done \(somalier could not tell the sex'

# 2 chrX sites: somalier keeps the pedigree's male, which is no call from the reads
FAKE_X_N=2 run_expect 0 qc-ped-kept "${SCRIPTS}/33-sample-qc.sh" sample1 male
output_has qc-ped-kept 'Sex check: not done \(somalier could not tell the sex from 2 chrX sites'
[ "$(tval inferred_sex)" = unknown ] || fail "the pedigree's sex somalier kept was read as its call: $(tval inferred_sex)"

# No declared sex: the pedigree says -9
run_expect 0 qc-ped-none "${SCRIPTS}/33-sample-qc.sh" sample1
[ "$(awk -F'\t' '{print $5}' "${GENOME_DIR}/sample1/qc/somalier/sample1.declared.ped")" = -9 ] \
  || fail "the pedigree of an undeclared sample is not -9"

# somalier's -2: chrX heterozygous like a female, yet chrY has reads
FAKE_SOMALIER_SEX=-2 FAKE_X_HET=150 run_expect 0 qc-x-and-y "${SCRIPTS}/33-sample-qc.sh" sample1 female
output_has qc-x-and-y 'Sex check: not done \(chrX is heterozygous like a female sample but chrY has reads'
FAKE_SOMALIER_SEX=-2 FAKE_X_HET=150 run_expect 0 qc-x-and-y-undeclared "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-x-and-y-undeclared 'Sex check: not done \(chrX is heterozygous like a female sample but chrY has reads'

FAKE_SOMALIER_SEX=2 FAKE_FREEMIX=0.12 run_expect 0 qc-contaminated "${SCRIPTS}/33-sample-qc.sh" sample1 female
output_has qc-contaminated 'FREEMIX 0\.1200 is above 0\.03: about 12% of the reads may come'
[ "$(tval contamination)" = warn ] || fail "FREEMIX 0.12 was not reported as a warning"

FAKE_FREEMIX=0.05 FREEMIX_WARN=0.08 run_expect 0 qc-threshold "${SCRIPTS}/33-sample-qc.sh" sample1
[ "$(tval contamination)" = ok ] || fail "FREEMIX_WARN=0.08 did not move the threshold"
FREEMIX_WARN=3 run_expect 1 qc-bad-threshold "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-bad-threshold 'FREEMIX_WARN must be a fraction between 0 and 1'

: > "$FAKE_DOCKER_LOG"
FAKE_VB2_FEW=true run_expect 0 qc-few-markers "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-few-markers 'Fewer than 1,000 panel markers have reads'
[ "$(grep -c ' verifybamid2 --SVDPrefix ' "$FAKE_DOCKER_LOG")" -eq 2 ] || fail "VerifyBamID2 did not run exactly twice"
docker_log_has 'verifybamid2 .*--DisableSanityCheck' "the second VerifyBamID2 call lacks --DisableSanityCheck"
[ "$(tval verifybamid2_marker_check)" = skipped ] || fail "the skipped marker check was not recorded"
[ "$(tval panel_markers)" = 320 ] || fail "the panel size was not recorded"

FAKE_VB2_FAIL=1 run_expect 1 qc-vb2-crash "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-vb2-crash 'ERROR: VerifyBamID2 failed'

FAKE_SM=SM_from_header run_expect 0 qc-sm-name "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-sm-name "somalier reads the sample as 'SM_from_header'"

mv "${GENOME_DIR}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu" "${CASE_WORK}/mu"
run_expect 1 qc-no-panel "${SCRIPTS}/33-sample-qc.sh" sample1
output_has qc-no-panel 'Install it with: ./scripts/setup.sh --sample-qc-data '
mv "${CASE_WORK}/mu" "${GENOME_DIR}/reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.mu"

# --- step 34 -----------------------------------------------------------------------------
A="${GENOME_DIR}/sample1/aligned"
FAKE_CRAM_READS=999 run_expect 1 cram-fewer "${SCRIPTS}/34-cram-archive.sh" sample1 --delete-bam
output_has cram-fewer 'The CRAM does not hold the same reads as the BAM'
[ -s "${A}/sample1_sorted.bam" ] || fail "the BAM was deleted although the CRAM held fewer reads"
[ ! -e "${A}/sample1_sorted.cram" ] && [ ! -e "${A}/sample1_sorted.part.cram" ] || fail "a CRAM that failed the check was left"

FAKE_QUICKCHECK_FAIL=cram run_expect 1 cram-truncated "${SCRIPTS}/34-cram-archive.sh" sample1 --delete-bam
output_has cram-truncated 'fails samtools quickcheck'
[ -s "${A}/sample1_sorted.bam" ] || fail "the BAM was deleted although the CRAM failed quickcheck"

run_expect 0 cram-keep "${SCRIPTS}/34-cram-archive.sh" sample1
[ -s "${A}/sample1_sorted.cram" ] && [ -e "${A}/sample1_sorted.cram.crai" ] || fail "no checked CRAM"
[ -s "${A}/sample1_sorted.bam" ] || fail "the BAM was deleted without --delete-bam"

# A second archive or restore while another run holds the sample's lock
# (aligned/sample1_sorted.lock) stops before it touches a file: no samtools
# call, the BAM and the CRAM as they were. This shell holds the lock on fd 8.
LOCK="${A}/sample1_sorted.lock"
command -v flock >/dev/null || fail "flock is not on PATH; the lock cases need it"
: >> "$LOCK"
exec 8<"$LOCK"
flock -n 8 || fail "the test could not take ${LOCK}"
CRAM_SUM=$(cksum < "${A}/sample1_sorted.cram")
: > "$FAKE_DOCKER_LOG"
run_expect 1 cram-locked "${SCRIPTS}/34-cram-archive.sh" sample1 --delete-bam
output_has cram-locked 'another archive or restore of sample1 is running'
# --restore checks the lock before it looks at the BAM, so a held lock is
# what it reports here, not the BAM that exists.
run_expect 1 cram-restore-locked "${SCRIPTS}/34-cram-archive.sh" sample1 --restore
output_has cram-restore-locked 'another archive or restore of sample1 is running'
output_lacks cram-restore-locked 'exists already'
if grep -q 'samtools' "$FAKE_DOCKER_LOG"; then fail "a run that found the lock held called samtools"; fi
[ -s "${A}/sample1_sorted.bam" ] && [ -e "${A}/sample1_sorted.bam.bai" ] || fail "a run that found the lock held deleted the BAM"
[ "$(cksum < "${A}/sample1_sorted.cram")" = "$CRAM_SUM" ] || fail "a run that found the lock held changed the CRAM"
exec 8<&-
# The lock goes with the run that holds it: once this shell lets go, the
# step runs. (The failed runs above, cram-fewer and cram-truncated, let go
# too: cram-keep after them took the lock.)

run_expect 0 cram-delete "${SCRIPTS}/34-cram-archive.sh" sample1 --delete-bam
[ ! -e "${A}/sample1_sorted.bam" ] && [ ! -e "${A}/sample1_sorted.bam.bai" ] || fail "--delete-bam left the BAM"

run_expect 0 cram-restore "${SCRIPTS}/34-cram-archive.sh" sample1 --restore
[ -s "${A}/sample1_sorted.bam" ] && [ -e "${A}/sample1_sorted.bam.bai" ] || fail "--restore wrote no BAM"
run_expect 1 cram-restore-again "${SCRIPTS}/34-cram-archive.sh" sample1 --restore
output_has cram-restore-again 'exists already; there is nothing to restore'

rm -f "${A}/sample1_sorted.bam" "${A}/sample1_sorted.bam.bai"
FAKE_CRAM_READS=999 run_expect 1 cram-restore-differs "${SCRIPTS}/34-cram-archive.sh" sample1 --restore
[ ! -e "${A}/sample1_sorted.bam" ] && [ ! -e "${A}/sample1_sorted.part.bam" ] || fail "--restore kept a BAM that does not match the CRAM"

FAKE_VIEW_FAIL=1 run_expect 1 cram-restore-write-fails "${SCRIPTS}/34-cram-archive.sh" sample1 --restore
output_has cram-restore-write-fails 'samtools could not write the BAM from the CRAM'
[ ! -e "${A}/sample1_sorted.part.bam" ] || fail "a failed --restore left its partial BAM"

run_expect 2 cram-bad-option "${SCRIPTS}/34-cram-archive.sh" sample1 --delete
echo "PASS: steps 33 and 34"
