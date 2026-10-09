#!/usr/bin/env bash
# run-all.sh with the default settings on a sample that has a BAM and a VCF,
# ClinVar installed and no optional data. run-all.sh is a launcher: it runs
# validate-setup.sh, writes a one-row samplesheet and starts the Nextflow
# pipeline (a fake nextflow here, which logs its arguments) with -resume, the
# --tools list of a default run minus the steps whose data is missing, and
# the database parameters it found, and the host's CPU count and RAM as the
# caps, since Nextflow refuses a task that asks for more. Then the two reports
# run. prs_scores/ exists but holds no scoring file, as after a failed first
# download: PRS is skipped, not started on nothing.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"

export GENOME_DIR="${CASE_WORK}/genome"
G=$GENOME_DIR
seed_reference "$G"
seed_clinvar "$G"
seed_sample "$G" sample1
mkdir -p "${G}/prs_scores"

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
TOOLS='pharmcat,cpic,roh,mito_haplogroup,mosdepth,telomere_hunter,mito_variants,manta,delly,duphold,survivor_merge,multiqc,clinvar,expansion_hunter,stranger'
CV="${G}/clinvar/clinvar_pathogenic_chr.vcf.gz"
want="nextflow :: cwd=$(cd "${G}/sample1/nextflow" && pwd) :: NXF_VER=26.04.7 :: $(printf '%q ' run "${REPO_ROOT}/main.nf" -profile docker -resume \
  --input "$SHEET" --reference "${G}/reference/GRCh38_no_alt_analysis_set.fasta" --outdir "$G" --tools "$TOOLS" \
  --clinvar "$CV" --clinvar_index "${CV}.tbi" --expansion_catalog "${G}/reference/expansionhunter_variant_catalog.json" \
  --max_cpus "$(host_cpus)" --max_memory "$(host_mem_gb).GB")"
[ "$NFLOG" = "$want" ] || fail "nextflow arguments differ:
  got:  ${NFLOG}
  want: ${want}"
grep -q '^NEXTFLOW_VERSION="26.04.7"' "${REPO_ROOT}/versions.env" || fail "versions.env no longer pins Nextflow 26.04.7: update this case's NXF_VER"

B="${G}/sample1/aligned/sample1_sorted.bam" V="${G}/sample1/vcf/sample1.vcf.gz"
[ "$(cat "$SHEET")" = "sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex
sample1,,,${B},${B}.bai,${V},${V}.tbi,male" ] || fail "samplesheet: $(cat "$SHEET")"
output_has run-all 'NOTE: starting from the existing VCF, with no gVCF with its index beside it'

# Exact counts: a step that turns from run into skipped, or back, fails here.
# Skipped: VEP, CPSR, CNVpytor, AnnotSV, pypgx, HLA and PRS (data not
# installed), vcfanno, clinical filter and slivar (need VEP), ancestry (needs
# PRS), and Cyrius, Parascopy and the Y haplogroup (opt-in: only with TOOLS
# naming them).
[ "$(grep -cE '^  [0-9]+b? .* runs$' "${CASE_WORK}/run-all.out")" -eq 15 ] || fail "not 15 steps run: $(grep -E ' runs$' "${CASE_WORK}/run-all.out" | tr '\n' '|')"
[ "$(grep -cE '^  [0-9]+b? .* skipped ' "${CASE_WORK}/run-all.out")" -eq 14 ] || fail "not 14 steps skipped"
output_has run-all '^  26 Ancestry \(pgsc_calc\) +skipped +\(needs PRS\)$'
output_has run-all '^  37 Y haplogroup \(Yleaf\) +skipped +\(opt-in: add y_haplogroup to TOOLS\)$'
output_has run-all '^  21 Cyrius CYP2D6 +skipped +\(opt-in: add cyrius to TOOLS\)$'
output_has run-all '^  35 Parascopy SMN1/SMN2 +skipped +\(opt-in: add parascopy to TOOLS\)$'
output_has run-all '^  31 slivar +skipped +\(needs VEP, data not installed: vep_cache/'
output_has run-all '^  25 PRS +skipped +\(data not installed: prs_scores/<PGS id>\.txt\.gz\)$'

STATUS="${G}/sample1/logs/run_status.tsv"
grep -q $'^meta\tdeclared_sex\tmale$' "$STATUS" || fail "run_status.tsv lacks the declared sex"
grep -q $'^step\t13\tskipped (data not installed: ' "$STATUS" || fail "run_status.tsv lacks step 13 skipped"
grep -q $'^step\t21\tskipped (opt-in: add cyrius to TOOLS)$' "$STATUS" || fail "run_status.tsv lacks step 21 skipped (opt-in)"
for s in 06 07 16 27 28; do grep -q $'^step\t'"${s}"$'\tok$' "$STATUS" || fail "run_status.tsv lacks step ${s} ok: $(cat "$STATUS")"; done
grep -q $'^run\twritten_by\trun-all.sh' "${G}/sample1/run_manifest.tsv" || fail "no run manifest from run-all.sh"
for l in 24-html-report generate-report; do [ -s "${G}/sample1/logs/${l}.log" ] || fail "the report ${l} did not run"; done
echo "run-all.sh validated the setup, wrote the samplesheet and started nextflow once with the default flags."
