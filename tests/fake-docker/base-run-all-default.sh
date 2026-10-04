#!/usr/bin/env bash
# run-all.sh with the default settings on a sample that has a BAM and a VCF,
# ClinVar installed and no optional data. run-all.sh is a launcher: it runs
# validate-setup.sh, writes a one-row samplesheet and starts the Nextflow
# pipeline (a fake nextflow here, which logs its arguments) with -resume, the
# --tools list of a default run minus the steps whose data is missing, and
# the database parameters it found. Then the two reports run.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1

# The ExpansionHunter catalog comes out of its image (`cat` in the container).
use_output_hook
cat > "${CASE_WORK}/tools-hook" <<'HOOK'
#!/usr/bin/env bash
case "${*:2}" in
  "cat /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json") echo '[{"LocusId": "HTT"}]' ;;
esac
exec "${CASE_WORK}/hook-outputs" "$@"
HOOK
chmod +x "${CASE_WORK}/tools-hook"
export FAKE_DOCKER_RUN_HOOK="${CASE_WORK}/tools-hook"

run_expect 0 run-all "${SCRIPTS}/run-all.sh" sample1 male
output_has run-all '=== Summary ==='   # validate-setup.sh ran first
output_lacks run-all 'unbound variable'

NFLOG=$(grep '^nextflow :: ' "$FAKE_DOCKER_LOG" || true)
[ "$(grep -c . <<<"$NFLOG")" -eq 1 ] || fail "nextflow was called $(grep -c . <<<"$NFLOG") times, expected once: ${NFLOG}"
SHEET="${G}/sample1/nextflow/samplesheet.csv"
TOOLS='pharmcat,cpic,roh,mito_haplogroup,mosdepth,telomere_hunter,mito_variants,cyrius,manta,delly,duphold,survivor_merge,multiqc,clinvar,expansion_hunter,stranger'
CV="${G}/clinvar/clinvar_pathogenic_chr.vcf.gz"
want="nextflow :: cwd=$(cd "${G}/sample1/nextflow" && pwd) :: NXF_VER=25.10.8 :: $(printf '%q ' run "${REPO_ROOT}/main.nf" -profile docker -resume \
  --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --outdir "$G" --tools "$TOOLS" \
  --clinvar "$CV" --clinvar_index "${CV}.tbi" --expansion_catalog "${G}/reference/expansionhunter_variant_catalog.json")"
[ "$NFLOG" = "$want" ] || fail "nextflow arguments differ:
  got:  ${NFLOG}
  want: ${want}"
grep -q '^NEXTFLOW_VERSION="25.10.8"' "${REPO_ROOT}/versions.env" || fail "versions.env no longer pins Nextflow 25.10.8: update this case's NXF_VER"

B="${G}/sample1/aligned/sample1_sorted.bam" V="${G}/sample1/vcf/sample1.vcf.gz"
[ "$(cat "$SHEET")" = "sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex
sample1,,,${B},${B}.bai,${V},${V}.tbi,male" ] || fail "samplesheet: $(cat "$SHEET")"
output_has run-all 'NOTE: starting from the existing VCF'

# Exact counts: a step that turns from run into skipped, or back, fails here.
# Skipped: VEP, CPSR, CNVpytor, AnnotSV, pypgx, HLA and PRS (data not
# installed), and vcfanno, clinical filter and slivar (need VEP).
[ "$(grep -cE '^  [0-9]+b? .* runs$' "${CASE_WORK}/run-all.out")" -eq 16 ] || fail "not 16 steps run: $(grep -E ' runs$' "${CASE_WORK}/run-all.out" | tr '\n' '|')"
[ "$(grep -cE '^  [0-9]+b? .* skipped ' "${CASE_WORK}/run-all.out")" -eq 10 ] || fail "not 10 steps skipped"
output_has run-all '^  31 slivar +skipped +\(needs VEP, data not installed: vep_cache/'
output_has run-all '^  25 PRS +skipped +\(data not installed: prs_scores\)'

STATUS="${G}/sample1/logs/run_status.tsv"
grep -q $'^meta\tdeclared_sex\tmale$' "$STATUS" || fail "run_status.tsv lacks the declared sex"
grep -q $'^step\t13\tskipped (data not installed: ' "$STATUS" || fail "run_status.tsv lacks step 13 skipped"
for s in 06 07 16 21 27 28; do grep -q $'^step\t'"${s}"$'\tok$' "$STATUS" || fail "run_status.tsv lacks step ${s} ok: $(cat "$STATUS")"; done
grep -q $'^run\twritten_by\trun-all.sh' "${G}/sample1/run_manifest.tsv" || fail "no run manifest from run-all.sh"
for l in 24-html-report generate-report; do [ -s "${G}/sample1/logs/${l}.log" ] || fail "the report ${l} did not run"; done
echo "run-all.sh validated the setup, wrote the samplesheet and started nextflow once with the default flags."
