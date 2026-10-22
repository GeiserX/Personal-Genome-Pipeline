#!/usr/bin/env bash
# Step 25's run options, against the fake docker and nextflow:
#   - a rerun keeps the previous run: its summary under its time
#     (<sample>_prs_summary.<YYYYMMDD-HHMM>.tsv, with the earlier content)
#     and its pgsc_calc folder as pgsc_calc.<the same time>, with the earlier
#     content; the new summary is written beside them;
#   - PRS_KEEP_AMBIGUOUS=true reaches prs-format (--keep-ambiguous with the
#     reference in the container) and pgsc_calc (--keep_ambiguous); without
#     it neither gets the flag;
#   - PRS_KEEP_AMBIGUOUS=true without the gVCF stops before any container
#     (WGS only), and a value other than true or false is refused;
#   - the declared sex reaches prs-table (--sex female).
# Control: a copy of step 25 with the "previous run" block cut out must lose
# the earlier summary and folder, so the checks above can fail.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

G="${CASE_WORK}/genome"
mkdir -p "$G"
export GENOME_DIR="$G" NXF_HOME="${CASE_WORK}/nxf"
seed_reference "$G"
seed_sample "$G" sample1
D="${G}/sample1/prs"
mkdir -p "${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}" "${NXF_HOME}/plugins/nf-schema-${PGSC_CALC_NF_SCHEMA}" "${G}/prs_scores"
: > "${G}/tools/pgsc_calc-${PGSC_CALC_VERSION}/main.nf"
for id in $(awk -F'\t' '$1 ~ /^PGS[0-9]+$/ {print $1}' "${REPO_ROOT}/assets/pgs_scores.tsv"); do
  printf '#HmPOS_build=GRCh38\nhm_chr\thm_pos\teffect_allele\teffect_weight\n20\t100\tA\t1\n' | gzip -c > "${G}/prs_scores/${id}.txt.gz"
done
printf 'placeholder\n' > "${G}/sample1/vcf/sample1.g.vcf.gz"
printf 'placeholder\n' > "${G}/sample1/vcf/sample1.g.vcf.gz.tbi"

# The tools' outputs: prs-format's score folder and allele list, one genotyped
# position from bcftools (on stdout, as step 25 reads it), and prs-table's
# summary, whose second row carries RUN_TAG so each run's table is told apart.
use_output_hook
cat > "${CASE_WORK}/hook-prs" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
shift   # the image
args=("$@")
val() { local i; for ((i = 0; i < ${#args[@]} - 1; i++)); do [ "${args[i]}" = "$1" ] && { host_path "${args[i + 1]}"; return; }; done; return 1; }
case " $* " in
  *" prs-format "*)
    out=$(val --out); alleles=$(val --alleles)
    mkdir -p "$out"; : > "${out}/PGS000018.txt.gz"
    printf 'chr20\t100\tA\nchr20\t100\tG\n' > "$alleles" ;;
  *" prs-table "*)
    out=$(val --out)
    printf 'Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tMatched_Pct\tPercentile\tAncestry_Group\tChrX\tAmbiguous\tApplies_To\tSex_Check\tInput\n' > "$out"
    printf 'Coronary artery disease\tPGS000018\t%s\t1\t1\t100.0\tNA\tNA\t0\t0\t\t\tgvcf\n' "${RUN_TAG:?}" >> "$out" ;;
  *" bcftools "*)
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\nchr20\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n' ;;
esac
exit 0
HOOK
chmod +x "${CASE_WORK}/hook-prs"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/hook-prs"

summary() { awk -F'\t' 'NR == 2 {print $3}' "$1" 2>/dev/null; }

# --- first run, female, default options ------------------------------------------------
: > "$FAKE_DOCKER_LOG"
RUN_TAG=first run_expect 0 first "${SCRIPTS}/25-prs.sh" sample1 female
[ "$(summary "${D}/sample1_prs_summary.tsv")" = first ] || fail "the first run wrote no summary"
FIRST_DATED=$(find "$D" -maxdepth 1 -name 'sample1_prs_summary.*.tsv' | head -n 1)
[ -n "$FIRST_DATED" ] || fail "the first run kept no dated copy of its summary"
[[ "$FIRST_DATED" =~ sample1_prs_summary\.[0-9]{8}-[0-9]{4}\.tsv$ ]] || fail "dated name is not <sample>_prs_summary.YYYYMMDD-HHMM.tsv: ${FIRST_DATED}"
docker_log_has 'prs-table .*--sex female' "the declared sex did not reach prs-table"
awk '/prs-format/ && /keep-ambiguous/ {f = 1} END {exit f}' "$FAKE_DOCKER_LOG" || fail "prs-format got --keep-ambiguous without PRS_KEEP_AMBIGUOUS"
awk '/^nextflow / && /keep_ambiguous/ {f = 1} END {exit f}' "$FAKE_DOCKER_LOG" || fail "pgsc_calc got --keep_ambiguous without PRS_KEEP_AMBIGUOUS"
docker_log_has '^nextflow .*run .*pgsc_calc' "step 25 did not start pgsc_calc"
echo first > "${D}/pgsc_calc/marker"
# Back-date the first run so the second gets another minute.
touch -t 202601010101 "${D}/sample1_prs_summary.tsv" "$FIRST_DATED"
mv "$FIRST_DATED" "${D}/sample1_prs_summary.202601010101.tsv.first"

# --- second run, PRS_KEEP_AMBIGUOUS=true ---------------------------------------------------
: > "$FAKE_DOCKER_LOG"
RUN_TAG=second PRS_KEEP_AMBIGUOUS=true run_expect 0 second "${SCRIPTS}/25-prs.sh" sample1
[ "$(summary "${D}/sample1_prs_summary.tsv")" = second ] || fail "the second run did not write its summary"
[ "$(summary "${D}/sample1_prs_summary.20260101-0101.tsv")" = first ] \
  || fail "the rerun did not keep the first summary as sample1_prs_summary.20260101-0101.tsv: $(ls "$D")"
[ "$(cat "${D}/pgsc_calc.20260101-0101/marker" 2>/dev/null)" = first ] \
  || fail "the rerun did not keep the first pgsc_calc folder as pgsc_calc.20260101-0101: $(ls "$D")"
[ ! -e "${D}/pgsc_calc/marker" ] || fail "the new pgsc_calc folder holds the old run's files"
output_has second 'previous run.s pgsc_calc folder is kept as .*/pgsc_calc\.20260101-0101'
docker_log_has 'prs-format .*--keep-ambiguous --reference /genome/reference/GRCh38_no_alt_analysis_set\.fasta' \
  "PRS_KEEP_AMBIGUOUS=true did not reach prs-format with the reference"
docker_log_has '^nextflow .*run .*--keep_ambiguous' "PRS_KEEP_AMBIGUOUS=true did not reach pgsc_calc"

# --- a third run in the same minute as an earlier summary ------------------------------------
# Dated summaries of this minute and the next, unless one is there already:
# the third run's own dated copy must take a new name, never overwrite them.
for s in "$(date +%Y%m%d-%H%M)" "$(python3 -c 'import time; print(time.strftime("%Y%m%d-%H%M", time.localtime(time.time() + 60)))')"; do
  [ -e "${D}/sample1_prs_summary.${s}.tsv" ] || printf 'h\nx\tx\tother\n' > "${D}/sample1_prs_summary.${s}.tsv"
done
mapfile -t EARLIER < <(cd "$D" && ls sample1_prs_summary.*.tsv)
BEFORE=$(cd "$D" && cksum "${EARLIER[@]}")
RUN_TAG=third run_expect 0 third "${SCRIPTS}/25-prs.sh" sample1
AFTER=$(cd "$D" && cksum "${EARLIER[@]}")
[ "$BEFORE" = "$AFTER" ] || fail "the third run overwrote an earlier dated summary: ${BEFORE} / ${AFTER}"
THIRD=$(for f in "$D"/sample1_prs_summary.*.tsv; do if [ "$(summary "$f")" = third ]; then echo "$f"; fi; done)
[[ "$THIRD" =~ sample1_prs_summary\.[0-9]{8}-[0-9]{4}-[0-9]+\.tsv$ ]] \
  || fail "the third run's dated copy did not take a -N name beside the earlier one of its minute: ${THIRD:-none}"

# --- refusals --------------------------------------------------------------------------------
mv "${G}/sample1/vcf/sample1.g.vcf.gz" "${CASE_WORK}/gvcf.hidden"
: > "$FAKE_DOCKER_LOG"
RUN_TAG=x PRS_KEEP_AMBIGUOUS=true run_expect 1 no-gvcf "${SCRIPTS}/25-prs.sh" sample1
output_has no-gvcf 'PRS_KEEP_AMBIGUOUS=true needs step 3.s gVCF'
if grep -qE '^run ' "$FAKE_DOCKER_LOG"; then fail "the refused run started a container"; fi
mv "${CASE_WORK}/gvcf.hidden" "${G}/sample1/vcf/sample1.g.vcf.gz"
RUN_TAG=x PRS_KEEP_AMBIGUOUS=yes run_expect 1 bad-value "${SCRIPTS}/25-prs.sh" sample1
output_has bad-value "PRS_KEEP_AMBIGUOUS must be 'true' or 'false'"

# --- control: step 25 without the "previous run" block -------------------------------------
# A copy of the checkout's scripts, bin, assets and versions.env (common.sh
# finds them beside itself), with the block cut out and the old rm -rf back.
B="${CASE_WORK}/broken"
mkdir -p "$B"
cp -R "${REPO_ROOT}/scripts" "${REPO_ROOT}/bin" "${REPO_ROOT}/assets" "${REPO_ROOT}/versions.env" "$B/"
awk '/^# >>> previous run/ {skip = 1} !skip {print} /^# <<< previous run/ {skip = 0}' "${SCRIPTS}/25-prs.sh" \
  | sed 's|^mkdir -p "\${WORK}/scores"$|rm -rf "$WORK"; mkdir -p "${WORK}/scores"|' > "${B}/scripts/25-prs.sh"
if grep -q '>>> previous run' "${B}/scripts/25-prs.sh" || ! grep -q '^rm -rf "\$WORK"; mkdir' "${B}/scripts/25-prs.sh"; then
  fail "control: the block was not cut out, or the old rm -rf was not put back"
fi
rm -rf "$D"
RUN_TAG=first run_expect 0 ctl-first bash "${B}/scripts/25-prs.sh" sample1
echo first > "${D}/pgsc_calc/marker"
find "$D" -maxdepth 1 -name 'sample1_prs_summary.*.tsv' -delete
touch -t 202601010101 "${D}/sample1_prs_summary.tsv"
RUN_TAG=second run_expect 0 ctl-second bash "${B}/scripts/25-prs.sh" sample1
if [ -n "$(find "$D" -maxdepth 1 -name 'pgsc_calc.*')" ] || [ -e "${D}/sample1_prs_summary.20260101-0101.tsv" ] \
   || [ -e "${D}/pgsc_calc/marker" ]; then
  fail "control: step 25 without the block still kept the first run: $(ls "$D")"
fi
echo "control: without the block the first run is lost, so the checks above can fail"
