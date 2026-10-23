#!/usr/bin/env bash
# run-all.sh passes step 20's two annotation files to the pipeline, so
# MITO_VARIANTS_ANNOTATE annotates the chrM calls as step 20 does:
#   both installed   --mito_clinvar clinvar/clinvar_chr.vcf.gz (with its .tbi)
#                    and --gnomad_mito gnomad_mito/<GNOMAD_MITO_URL's file>
#   neither          no flag: the module then names each file as not given
#   ClinVar without its index   no --mito_clinvar (main.nf would stop on the
#                    missing .tbi), --gnomad_mito still passed
# A fake nextflow logs its arguments.
# shellcheck source=../../scripts/ci/fake-docker/lib.sh
. "${REPO_ROOT:?}/scripts/ci/fake-docker/lib.sh"
# shellcheck source=../../versions.env
. "${REPO_ROOT}/versions.env"

export GENOME_DIR="${CASE_WORK}/genome" SKIP_VALIDATION=true
G=$GENOME_DIR
seed_reference "$G"
seed_sample "$G" sample1
use_output_hook   # the reports at the end read the files their containers write
last() { grep '^nextflow :: ' "$FAKE_DOCKER_LOG" | tail -1; }
has_args() { grep -qF -- "$(printf '%q ' "$@")" <<<"$(last)" || fail "nextflow lacks '$*': $(last)"; }
lacks_arg() { if grep -qF -- " $1 " <<<"$(last)"; then fail "$1 passed $2: $(last)"; fi; }
MC="${G}/clinvar/clinvar_chr.vcf.gz" GM="${G}/gnomad_mito/${GNOMAD_MITO_URL##*/}"

# Neither file: no flag
run_rc none "${SCRIPTS}/run-all.sh" sample1 male
grep -q '^nextflow :: ' "$FAKE_DOCKER_LOG" || fail "run-all.sh never started nextflow"
lacks_arg --mito_clinvar "without clinvar/clinvar_chr.vcf.gz"
lacks_arg --gnomad_mito "without gnomad_mito/"

# Both installed
seed_clinvar "$G"
mkdir -p "${G}/gnomad_mito"
printf 'placeholder\n' > "$GM"
run_rc both "${SCRIPTS}/run-all.sh" sample1 male
has_args --mito_clinvar "$MC"
has_args --gnomad_mito "$GM"

# ClinVar without its index
rm -f "${MC}.tbi"
run_rc noindex "${SCRIPTS}/run-all.sh" sample1 male
lacks_arg --mito_clinvar "without clinvar_chr.vcf.gz.tbi"
has_args --gnomad_mito "$GM"
echo "run-all.sh passes --mito_clinvar and --gnomad_mito when both are installed, and neither when not"
