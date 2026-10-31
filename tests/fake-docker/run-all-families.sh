#!/usr/bin/env bash
# run-all.sh with FAMILIES=<families file>: the file says who is related to
# whom, and run-all.sh runs the family steps when their members are ready.
#   1. partner_x of family_b (tests/fixtures/families/two_families.tsv) runs
#      first: bin/family_report.py plan writes GENOME_DIR/family_b/family_plan.tsv,
#      and steps 38 and 39 are recorded as skipped, waiting for partner_y's
#      gVCF; neither script runs;
#   2. partner_y runs next: both members have a gVCF and ClinVar hits, so
#      step 38 (38-family.sh family_b <file>) and step 39 run, and both are
#      recorded ok;
#   3. partner_x again with the sex the file contradicts: exit 2 before
#      nextflow starts;
#   4. the real 39-couple.sh given the families file writes family_b's own
#      rows to GENOME_DIR/family_b/pedigree.tsv and hands that file, not the
#      two-family one, to bin/couple_report.py.
# run-all.sh runs from a copy of the repository whose 38-family.sh and
# 39-couple.sh only log their arguments; family_report.py runs for real, with
# the host's python3, through a hook on the fake docker.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
FAMS="${CASE_WORK}/families.tsv"  # outside GENOME_DIR, as a user may keep it
cp "${REPO_ROOT}/tests/fixtures/families/two_families.tsv" "$FAMS"

# The copy: tracked and new files, with the two family scripts logging their arguments.
R="${CASE_WORK}/repo"
mkdir -p "$R"
(cd "$REPO_ROOT" && git ls-files -z -co --exclude-standard | tar --null -T - -cf -) | tar -xf - -C "$R"
for n in 38-family 39-couple; do
  printf '#!/usr/bin/env bash\necho "%s $*" >> "${CASE_WORK}/family-steps.log"\n' "${n%%-*}" > "${R}/scripts/${n}.sh"
done
STEPS_LOG="${CASE_WORK}/family-steps.log"

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
case "${*:2}" in
  "python3 /pgp-bin/family_report.py"*)
    # The real family_report.py, on the host paths of every /... argument.
    shift
    mapped=()
    for a in "$@"; do
      if [[ "$a" == /* ]] && h=$(host_path "$a"); then mapped+=("$h"); else mapped+=("$a"); fi
    done
    exec "${mapped[@]}" ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

# member S: what the pipeline leaves for S: the BAM and VCF, the gVCF with its index, step 06's hits
member() {
  seed_sample "$G" "$1"
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$1" | gzip -c > "${G}/$1/vcf/$1.g.vcf.gz"
  : > "${G}/$1/vcf/$1.g.vcf.gz.tbi"
  mkdir -p "${G}/$1/clinvar"
  printf '##fileformat=VCFv4.2\n' > "${G}/$1/clinvar/$1_clinvar_hits.vcf"
}
status() { cat "${G}/$1/logs/run_status.tsv"; }
nf_runs() { grep -c '^nextflow :: .* run ' "$FAKE_DOCKER_LOG" || true; }
PLAN="${G}/family_b/family_plan.tsv"

# --- 1. the first member: the family steps wait -----------------------------------------
member partner_x
run_rc first env FAMILIES="$FAMS" "${R}/scripts/run-all.sh" partner_x male
[ "$RC" -le 1 ] || fail "the first run exited ${RC}, expected 0 or 1"
output_has first '^\[Family\] partner_x is in family family_b'
[ -s "$PLAN" ] || fail "no ${PLAN} after the first run"
grep -q $'^family_b\tjoint_call\tpartner_x,partner_y\twaiting\tpartner_y has no gVCF yet' "$PLAN" \
  || fail "the plan does not wait for partner_y's gVCF: $(cat "$PLAN" 2>/dev/null)"
status partner_x | grep -q $'^step\t38\tskipped (waiting: partner_y has no gVCF yet' \
  || fail "step 38 is not recorded as waiting for partner_y: $(status partner_x)"
status partner_x | grep -q $'^step\t39\tskipped (waiting: partner_y has no gVCF yet' \
  || fail "step 39 is not recorded as waiting for partner_y: $(status partner_x)"
[ ! -e "$STEPS_LOG" ] || fail "a family script ran while partner_y has no gVCF: $(cat "$STEPS_LOG")"
# run-all.sh never passes the other members to the pipeline
grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | grep -q partner_y && fail "partner_y reached the pipeline of partner_x's run"
grep -q '^partner_x,' "${G}/partner_x/nextflow/samplesheet.csv" || fail "the samplesheet is not partner_x's"
[ "$(wc -l < "${G}/partner_x/nextflow/samplesheet.csv")" -eq 2 ] || fail "the samplesheet has more than partner_x's row"

# --- 2. the second member: both steps run --------------------------------------------------
member partner_y
run_rc second env FAMILIES="$FAMS" "${R}/scripts/run-all.sh" partner_y female
[ "$RC" -le 1 ] || fail "the second run exited ${RC}, expected 0 or 1"
grep -q $'^family_b\tjoint_call\tpartner_x,partner_y\truns\t' "$PLAN" || fail "the plan does not run the joint call: $(cat "$PLAN")"
grep -q $'^family_b\tcouple_carrier\tpartner_x,partner_y\truns\t' "$PLAN" || fail "the plan does not run the couple screen: $(cat "$PLAN")"
grep -qx "38 family_b ${FAMS}" "$STEPS_LOG" || fail "38-family.sh did not run on family_b with the file: $(cat "$STEPS_LOG" 2>/dev/null)"
grep -qx "39 family_b ${FAMS}" "$STEPS_LOG" || fail "39-couple.sh did not run on family_b with the file: $(cat "$STEPS_LOG" 2>/dev/null)"
status partner_y | grep -q $'^step\t38\tok$' || fail "step 38 is not recorded ok: $(status partner_y)"
status partner_y | grep -q $'^step\t39\tok$' || fail "step 39 is not recorded ok: $(status partner_y)"
if status partner_y | grep -q $'\tplanned$'; then fail "a planned row is left: $(status partner_y)"; fi
[ "$(grep -c '' "$STEPS_LOG")" -eq 2 ] || fail "the family scripts ran $(grep -c '' "$STEPS_LOG") times, want 2"

# --- 3. a sex the file contradicts stops before nextflow -------------------------------------
before=$(nf_runs)
run_expect 2 badsex env FAMILIES="$FAMS" "${R}/scripts/run-all.sh" partner_x female
output_has badsex "the pedigree says 'partner_x' is male, the sex given to run-all.sh says female"
[ "$(nf_runs)" -eq "$before" ] || fail "nextflow ran although the sex contradicts the families file"
run_expect 2 nofile env FAMILIES="${CASE_WORK}/none.tsv" "${R}/scripts/run-all.sh" partner_x male
output_has nofile 'FAMILIES=.*no such file'

# --- 4. the real step 39 reads only its family's rows -----------------------------------------
rm -f "${G}/family_b/pedigree.tsv"
run_expect 0 couple "${REPO_ROOT}/scripts/39-couple.sh" family_b "$FAMS"
[ "$(cut -f1 "${G}/family_b/pedigree.tsv" 2>/dev/null | tr '\n' ' ')" = "sample partner_x partner_y " ] \
  || fail "family_b/pedigree.tsv is not family_b's rows: $(cat "${G}/family_b/pedigree.tsv" 2>/dev/null)"
cmp -s "$FAMS" "${G}/family_b/families.tsv" || fail "the families file was not kept beside it as families.tsv"
docker_log_has '/pgp-bin/couple_report\.py screen --pedigree /genome/family_b/pedigree\.tsv ' \
  "couple_report.py was not given family_b's own pedigree"
docker_log_has '--clinvar partner_y=/genome/partner_y/clinvar/partner_y_clinvar_hits\.vcf' \
  "step 39 did not hand partner_y's ClinVar hits over"
if grep 'couple_report' "$FAKE_DOCKER_LOG" | grep -q child_a; then fail "step 39 handed family_a's child_a over"; fi
run_expect 2 nofam "${REPO_ROOT}/scripts/39-couple.sh" family_c "$FAMS"
output_has nofam "has no family 'family_c'; its families are: family_a, family_b"
echo "FAMILIES: the family steps waited for the second member, then ran on family_b's own rows."
