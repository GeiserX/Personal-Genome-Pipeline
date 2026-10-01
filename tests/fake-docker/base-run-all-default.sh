#!/usr/bin/env bash
# run-all.sh with the default settings (SKIP_VALIDATION unset, no optional
# data) on a sample that already has a BAM and a VCF: the pre-flight
# validation passes, the steps whose data is not installed are skipped, and
# the run exits 0.
#
# The hook stands in for the tools: it writes the files each step checks
# after its container (Manta's VCFs, goleft's .ped, Cyrius's .tsv, PharmCAT's
# JSON, plink2's .sscore) and answers the bcftools calls whose output a step
# reads (contig list, record count). Everything else goes to the generic
# output hook from lib.sh.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
seed_reference "$GENOME_DIR"
seed_clinvar "$GENOME_DIR"
seed_sample "$GENOME_DIR" sample1

# Step 25 downloads PGS Catalog scoring files and refuses one without a
# GRCh38 #HmPOS_build header. Seed harmonised files so it scores, as it does
# on a second run.
mkdir -p "${GENOME_DIR}/prs_scores"
for id in PGS000018 PGS000014 PGS000004 PGS000662 PGS000016 PGS000334 PGS000027 PGS000017 PGS000055; do
  printf '#HmPOS_build=GRCh38\neffect_allele\teffect_weight\thm_chr\thm_pos\nA\t0.1\t1\t1000\n' \
    | gzip -c > "${GENOME_DIR}/prs_scores/${id}.txt.gz"
done

use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
. "${CASE_WORK:?}/host-path.sh"
args="${*:2}"
# opt NAME: the word after NAME in the command (a bash -c body spans lines),
# quotes stripped.
opt() {
  local -a w
  read -r -d '' -a w <<<"$args" || true
  local i
  for ((i = 0; i < ${#w[@]} - 1; i++)); do
    if [ "${w[i]}" = "$1" ]; then printf '%s' "${w[i + 1]//[\'\"]/}"; return 0; fi
  done
  return 1
}
put() {   # put CONTAINER_PATH CONTENT
  local h
  h=$(host_path "$1")
  mkdir -p "$(dirname "$h")"
  printf '%b' "$2" > "$h"
}
case "$args" in
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
output_lacks run-all 'Setup validation failed'
output_lacks run-all 'unbound variable'
output_has run-all ' skipped, 0 failed$'
