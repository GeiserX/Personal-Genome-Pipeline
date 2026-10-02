#!/usr/bin/env bash
# run-all.sh on a sample that kept step 31's summary from an earlier run:
# this run skips step 31 (VEP is not installed), and both reports must mark the
# old slivar result stale, with its date, instead of showing it as current.
#
# The tools are faked as in base-run-all-default.sh, except the report
# renderer: the hook runs bin/render_report.py with the host's python3 on the
# host side of the call's mounts, so the reports are the real ones.
# Also checks the run manifest and the step status file run-all.sh writes.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1
mkdir -p "${GENOME_DIR}/prs_scores"
for id in PGS000018 PGS000014 PGS000004 PGS000662 PGS000016 PGS000334 PGS000027 PGS000017 PGS000055; do
  printf '#pgs_id=%s\n#HmPOS_build=GRCh38\neffect_allele\teffect_weight\thm_chr\thm_pos\nA\t0.1\t1\t1000\n' "$id" \
    | gzip -c > "${GENOME_DIR}/prs_scores/${id}.txt.gz"
done

# Step 31's output from an earlier run, 30 days old.
OLD="${GENOME_DIR}/sample1/slivar/sample1_slivar_summary.tsv"
mkdir -p "$(dirname "$OLD")"
printf 'CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\nchr1\t100\tA\tG\tHIGH\tOLDGENE\n' > "$OLD"
touch -d '30 days ago' "$OLD"

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args="${*:2}"
opt() {
  local -a w
  read -r -d '' -a w <<<"$args" || true
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]//[\'\"]/}"; return 0; fi
  done
  return 1
}
put() {
  local h
  h=$(host_path "$1")
  mkdir -p "$(dirname "$h")"
  printf '%b' "$2" > "$h"
}
case "$args" in
  "python3 /pgp-bin/render_report.py"*)
    # The real renderer, on the host paths of every /... argument.
    shift
    mapped=()
    for a in "$@"; do
      if [[ "$a" == /* ]] && h=$(host_path "$a"); then mapped+=("$h"); else mapped+=("$a"); fi
    done
    exec "${mapped[@]}" ;;
  configManta.py*)
    put "$(opt --runDir)/runWorkflow.py" '#!/usr/bin/env python\n' ;;
  */runWorkflow.py*)
    d="$(dirname "$2")/results/variants"
    for f in diploidSV candidateSV candidateSmallIndels; do
      h=$(host_path "${d}/${f}.vcf.gz")
      mkdir -p "$(dirname "$h")"
      printf '##fileformat=VCFv4.1\n' | gzip -c > "$h"
      : > "${h}.tbi"
    done ;;
  "goleft indexcov"*)
    d=$(opt --directory)
    put "${d}/${d##*/}-indexcov.ped" \
      '#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphenotype\tCNchrX\tCNchrY\nsample1\tsample1\t-9\t-9\t1\t-9\t1.02\t0.98\n' ;;
  *"cyrius "*)
    put "$(opt --outDir)/$(opt --prefix).tsv" 'Sample\tGenotype\tFilter\nsample1\t*1/*1\tPASS\n' ;;
  *pharmcat.jar*)
    put "$(opt -o)/$(opt -bf).report.json" '{}\n' ;;
  "plink2 "*--score*)
    put "$(opt --out).sscore" '#IID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM\nsample1\t2\t1\t0.05\t0.1\n' ;;
  "bcftools index -s "*)
    printf 'chr1\t248956422\t1\n' ;;
  "bcftools index -n "*|*"| wc -l"*)
    echo 1 ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

run_expect 0 run-all "${SCRIPTS}/run-all.sh" sample1 male
output_has run-all '^  31 slivar +skipped +needs VEP'

S="${GENOME_DIR}/sample1"
STATUS="${S}/logs/run_status.tsv"
[ -f "$STATUS" ] || fail "run-all.sh wrote no ${STATUS}"
grep -q $'^step\t31\tskipped' "$STATUS" || fail "run_status.tsv does not record step 31 as skipped: $(cat "$STATUS")"
grep -q $'^step\t06\tok' "$STATUS" || fail "run_status.tsv does not record step 06 as ok"
grep -q $'^meta\tdeclared_sex\tmale' "$STATUS" || fail "run_status.tsv does not record the declared sex"

MANIFEST="${S}/run_manifest.tsv"
[ -f "$MANIFEST" ] || fail "run-all.sh wrote no ${MANIFEST}"
grep -q $'^run\twritten_by\trun-all.sh' "$MANIFEST" || fail "the manifest does not say run-all.sh wrote it"
want_images=$(grep -cE '^[A-Z0-9_]+_IMAGE=' "${REPO_ROOT}/versions.env")
got_images=$(grep -c '^image' "$MANIFEST" || true)
[ "$got_images" -eq "$want_images" ] || fail "the manifest lists ${got_images} images, versions.env has ${want_images}"
grep -q $'^data\tpgs:PGS000018\tpgs_id=PGS000018;HmPOS_build=GRCh38' "$MANIFEST" \
  || fail "the manifest has no header line of PGS000018: $(grep pgs: "$MANIFEST" | head -2)"

TXT="${S}/sample1_report.txt"
HTML="${S}/sample1_report.html"
[ -s "$TXT" ] || fail "no text report"
[ -s "$HTML" ] || fail "no HTML report"
python3 "${REPO_ROOT}/tests/schema/validate.py" "${REPO_ROOT}/tests/schema/summary.schema.json" "${S}/summary.json" \
  || fail "summary.json does not validate"
grep -A2 '^## Variant Prioritization (slivar)' "$TXT" | grep -q 'STALE: from an earlier run.*step 31 was skipped' \
  || fail "the text report does not mark the old slivar result stale: $(grep -A3 'slivar' "$TXT" | tr '\n' '|')"
grep -q '<div class="stale">STALE: from an earlier run.*step 31 was skipped' "$HTML" \
  || fail "the HTML report does not mark the old slivar result stale"
# Steps this run did run are current, not stale.
if grep -A2 '^## Runs of Homozygosity' "$TXT" | grep -q STALE; then
  fail "ROH (step 11, ok in this run) is marked stale"
fi
echo "Both reports mark step 31's old result stale; the manifest and the status file are complete."
