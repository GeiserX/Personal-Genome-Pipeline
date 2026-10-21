#!/usr/bin/env bash
# Step 33 (sample identity, sex from the reads, nuclear contamination) in a
# default run-all.sh, and validate-setup.sh's warning about an index older
# than its data file.
#   1. With somalier's sites and VerifyBamID2's three panel files where
#      setup.sh installs them, step 33 runs: sample_qc is in --tools and the
#      pipeline gets --somalier_sites and --verifybamid2_panel.
#   2. Without the panel's .mu file it is skipped as 'data not installed',
#      and neither argument is passed (the control for 1).
#   3. A failed SOMALIER task (a process not named after the tool) marks
#      step 33 failed, not ok.
#   4. validate-setup.sh warns about a .tbi older than its VCF (then about
#      the gVCF's too, then about extra.bai beside extra.bam), and says
#      nothing about a pair in the right order or about an index in the
#      Nextflow work directory.
# A fake nextflow stands in; a wrapper in front of it writes the two files
# main.nf's completion handler would (as in run-all-kir-step.sh).
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
ST="${G}/sample1/logs/run_status.tsv"
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
SQ="${G}/reference/somalier/sites.hg38.vcf.gz" VB="${G}/reference/verifybamid2" VP=1000g.phase3.100k.b38.vcf.gz.dat
mkdir -p "${G}/reference/somalier" "$VB"
printf 'placeholder\n' > "$SQ"
for e in UD mu bed; do printf 'placeholder\n' > "${VB}/${VP}.${e}"; done

mkdir -p "${CASE_WORK}/nfbin"
cat > "${CASE_WORK}/nfbin/nextflow" <<'NF'
#!/usr/bin/env bash
[ -z "${PROCESS_STATUS:-}" ] || printf 'process\tsucceeded\tcached\tfailed\tignored\tunfinished\n%b' "$PROCESS_STATUS" > process_status.tsv
[ -z "${FAILED_TASKS:-}" ] || printf 'process\tstatus\n%b' "$FAILED_TASKS" > failed_tasks.tsv
exec "${REPO_ROOT}/scripts/ci/fake-docker/nextflow" "$@"
NF
chmod +x "${CASE_WORK}/nfbin/nextflow"
NFPATH="${CASE_WORK}/nfbin:${PATH}"
last_nf() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -n 1; }

# --- 1. both files installed: step 33 runs -------------------------------------------
run_expect 0 qc-on "${SCRIPTS}/run-all.sh" sample1 male
output_has qc-on '^  33 Sample identity and QC +runs$'
NFLINE=$(last_nf)
grep -qE -- "--tools [^ ]*sample_qc" <<< "$NFLINE" || fail "sample_qc is not in --tools: ${NFLINE}"
grep -qF -- "--somalier_sites $(printf '%q' "$SQ") --verifybamid2_panel $(printf '%q' "$VB") " <<< "$NFLINE" \
  || fail "nextflow was not given the sites and the panel folder: ${NFLINE}"
grep -q $'^step\t33\tok$' "$ST" || fail "run_status.tsv lacks step 33 ok: $(cat "$ST")"

# --- 2. the panel incomplete: skipped, no argument (the control) ----------------------
rm "${VB}/${VP}.mu"
run_expect 0 qc-off "${SCRIPTS}/run-all.sh" sample1 male
output_has qc-off "^  33 Sample identity and QC +skipped +\\(data not installed: reference/verifybamid2/${VP//./\\.}\\.mu\\)$"
NFLINE=$(last_nf)
! grep -qE -- "sample_qc|--somalier_sites|--verifybamid2_panel" <<< "$NFLINE" \
  || fail "nextflow got step 33 without its panel: ${NFLINE}"
grep -q $'^step\t33\tskipped (data not installed: ' "$ST" || fail "run_status.tsv lacks step 33 skipped: $(cat "$ST")"
printf 'placeholder\n' > "${VB}/${VP}.mu"

# --- 3. a failed SOMALIER task fails step 33 ------------------------------------------
row() { printf 'BAM_ANALYSIS:%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "$@"; }
run_expect 3 qc-failed env PATH="$NFPATH" FAKE_NEXTFLOW_RC=3 TOOLS=sample_qc FAILED_TASKS='BAM_ANALYSIS:SOMALIER\tfailed\n' \
  PROCESS_STATUS="$(row SOMALIER 0 0 1 0 0)$(row VERIFYBAMID2 1 0 0 0 0)$(row SOMALIER_RELATE 0 0 0 0 0)$(row SAMPLE_QC 0 0 0 0 0)" \
  "${SCRIPTS}/run-all.sh" sample1 male
got=$(awk -F'\t' '$2 == "33" {print $3}' "$ST")
[ "$got" = failed ] || fail "step 33 after a failed SOMALIER task is '${got}', want 'failed'"
output_has qc-failed 'sample_qc \(SOMALIER\) failed'

# --- 4. validate-setup.sh: an index older than its data file ---------------------------
# (The BAM pair keeps its order: with no @HD SO: tag an older .bai is a FAIL of its own.)
V="${G}/sample1/vcf/sample1.vcf.gz" GV="${G}/sample1/vcf/sample1.g.vcf.gz"
printf 'placeholder\n' > "$GV"; printf 'placeholder\n' > "${GV}.tbi"
touch -d 2026-01-01T00:00:00 "${V}.tbi" "$GV" "${GV}.tbi"
touch -d 2026-01-02T00:00:00 "$V"   # the VCF is newer than its .tbi; the gVCF pair has one time
mkdir -p "${G}/sample1/nextflow/work/ab"
printf 'placeholder\n' > "${G}/sample1/nextflow/work/ab/x.bam.bai"
printf 'placeholder\n' > "${G}/sample1/nextflow/work/ab/x.bam"
touch -d 2026-01-01T00:00:00 "${G}/sample1/nextflow/work/ab/x.bam.bai"
run_expect 0 stale "${SCRIPTS}/validate-setup.sh" sample1
output_has stale '\[WARN\].*Index older than its data file: sample1/vcf/sample1\.vcf\.gz\.tbi$'
output_has stale 'index it again \(samtools index for a BAM or CRAM, bcftools index for a VCF\); if it was only copied, touch the index'
[ "$(grep -c 'Index older than its data file: ' "${CASE_WORK}/stale.out")" -eq 1 ] \
  || fail "validate-setup.sh named an index in the right order, or one in the work directory"
touch -d 2026-01-03T00:00:00 "$GV"   # now the gVCF's .tbi is older too
run_expect 0 stale-gvcf "${SCRIPTS}/validate-setup.sh" sample1
output_has stale-gvcf 'Index older than its data file: sample1/vcf/sample1\.g\.vcf\.gz\.tbi$'
[ "$(grep -c 'Index older than its data file: ' "${CASE_WORK}/stale-gvcf.out")" -eq 2 ] || fail "validate-setup.sh did not name both older indexes"
# An index named sample.bai belongs to sample.bam (samtools reads that name too).
printf 'placeholder\n' > "${G}/sample1/aligned/extra.bai"; printf 'placeholder\n' > "${G}/sample1/aligned/extra.bam"
touch -d 2026-01-01T00:00:00 "${G}/sample1/aligned/extra.bai"
run_expect 0 stale-bai "${SCRIPTS}/validate-setup.sh" sample1
output_has stale-bai 'Index older than its data file: sample1/aligned/extra\.bai$'
touch "${V}.tbi" "${GV}.tbi" "${G}/sample1/aligned/extra.bai"
run_expect 0 fresh "${SCRIPTS}/validate-setup.sh" sample1
output_lacks fresh 'Index older than its data file'
echo "Step 33 runs in a default run when its data is installed, is skipped without it, and validate-setup.sh names an index older than its data."
