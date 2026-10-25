#!/usr/bin/env bash
# 44-ancient.sh — Where your genome falls among ancient genomes (AADR, smartpca)
# Usage: ./scripts/44-ancient.sh <sample_name>
#
# Places the sample on a PCA of present-day West Eurasians computed with
# smartpca (EIGENSOFT), onto which the ancient people of the Allen Ancient DNA
# Resource (AADR) are projected, and lists the ancient groups nearest to it
# with their dates and places. Similarity on a PCA is not descent; qpAdm
# models are not part of this step.
#
# Needs the panel `scripts/setup.sh --ancient-panel <genome_dir>` builds
# (reference/ancient/aadr_<release>, ANCIENT_PANEL_DIR to use another, none to
# skip). Without it the step says so in one line and exits 0.
#
# Reads the gVCF of step 3 (vcf/<sample>.g.vcf.gz), so a site where the
# sample matches the reference is a real genotype and a site without reads is
# missing, never taken for the reference. The same commands run in the
# Nextflow ANCIENT_SITES, ANCIENT_GENO, ANCIENT_PCA and ANCIENT_REPORT processes.
#
# Output in ancient/:
#   <sample>_ancient.tsv       the 20 nearest ancient groups (rank, group, n,
#                              distance, pc1, pc2, date_bp, dates, country, locality)
#   <sample>_ancient_info.tsv  release, sites used, people, the sample's PC1 and PC2
#   <sample>_ancient.txt       the same in words, with the caveat
#   <sample>_ancient_pca.svg   PC1 against PC2, the sample marked (dark background)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image BCFTOOLS_IMAGE PYTHON_IMAGE EIGENSOFT_IMAGE

PANEL=${ANCIENT_PANEL_DIR:-${GENOME_DIR}/reference/ancient/aadr_${AADR_VERSION}}
if [ "$PANEL" = none ] || [ ! -s "${PANEL}/DONE" ]; then
  echo "Step 44 skipped: no ancient populations panel$([ "$PANEL" = none ] && echo " (ANCIENT_PANEL_DIR=none)" || echo " at ${PANEL}"); install it with scripts/setup.sh --ancient-panel ${GENOME_DIR}"
  exit 0
fi
S="${GENOME_DIR}/${SAMPLE}"
GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
OUTDIR="${S}/ancient"
WORK="${OUTDIR}/work"
RELEASE=$(awk -F'\t' '$1 == "release" {print $2}' "${PANEL}/DONE")

echo "============================================"
echo "  Step 44: Ancient populations (AADR ${RELEASE}, smartpca)"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${GVCF}"
echo "  Panel:  ${PANEL}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
for f in "$GVCF" "${GVCF}.tbi" "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}. Step 44 reads the gVCF of step 3 (a variant-only VCF cannot tell a reference match from a site without reads)." >&2; exit 1; }
done
P=$(cpath "$PANEL") W=$(cpath "$WORK") O=$(cpath "$OUTDIR")
mkdir -p "$WORK"
rm -f "${OUTDIR}/${SAMPLE}"_ancient*

echo "[1/4] Genotyping the panel's sites from the gVCF..."
# gvcf2vcf expands each reference block over a site into a record with the
# reference base, -T keeps the panel's sites, --trim-alt-alleles drops <*>,
# and a no-call is dropped, so a site without reads has no line.
# shellcheck disable=SC2016  # $1 to $4 belong to the inner bash
run_in --cpus 2 --memory 4g "$BCFTOOLS_IMAGE" \
  bash -euo pipefail -c '
    n=$(bcftools query -l "$3" | wc -l | tr -d " ")
    [ "$n" -eq 1 ] || { echo "ERROR: $3 holds $n samples; step 44 reads a one-sample gVCF (keep one with bcftools view -s <sample>)." >&2; exit 1; }
    awk -F"\t" -v OFS="\t" "NR > 1 {print \$2, \$3}" "$2/sites_grch38.tsv" > "$4/sites.tsv"
    bcftools convert --gvcf2vcf -f "$1" -R "$4/sites.tsv" -Ou "$3" \
      | bcftools view -T "$4/sites.tsv" --trim-alt-alleles -i "GT!=\"mis\"" -Ou \
      | bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t[%TGT]\n" > "$4/genotypes.tsv"
    rm -f "$4/sites.tsv"' \
  _ "$REF_FASTA_C" "$P" "$(cpath "$GVCF")" "$W" \
  || { echo "ERROR: step 44 could not genotype ${GVCF} at the panel's sites (the line above says why)." >&2; exit 1; }

echo "[2/4] Adding ${SAMPLE} to the present-day panel..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/ancient_report.py person-geno --sites "${P}/sites_grch38.tsv" --genotypes "${W}/genotypes.tsv" \
    --snp "${P}/moderns.snp" --moderns-geno "${P}/moderns.geno" --out "${W}/merged.geno" --stats "${W}/stats.tsv"
{ cat "${PANEL}/moderns.ind"; printf '%s\tU\tPerson\n' "$SAMPLE"; } > "${WORK}/merged.ind"

echo "[3/4] smartpca: the present-day axes, ${SAMPLE} projected (lsqproject)..."
# The parameters of setup.sh's ancient_smartpca_par: axes from the panel's
# present-day groups, the sample projected, no outlier removal.
printf '%s\n' "genotypename: ${W}/merged.geno" "snpname: ${P}/moderns.snp" "indivname: ${W}/merged.ind" \
  "poplistname: ${P}/poplist.txt" "evecoutname: ${W}/run.evec" "evaloutname: ${W}/run.eval" "lsqproject: YES" \
  "numoutevec: 10" "numoutlieriter: 0" "hashcheck: NO" "numthreads: ${THREADS}" > "${WORK}/run.par"
run_in --cpus "$THREADS" --memory 8g "$EIGENSOFT_IMAGE" smartpca -p "${W}/run.par" > "${WORK}/smartpca.log" 2>&1 \
  || { tail -n 20 "${WORK}/smartpca.log" >&2; echo "ERROR: smartpca failed (log: ${WORK}/smartpca.log)." >&2; exit 1; }

echo "[4/4] The nearest ancient groups..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/ancient_report.py report --panel-evec "${P}/panel.evec" --run-evec "${W}/run.evec" \
    --ancients "${P}/ancients.tsv" --sample "$SAMPLE" --release "$RELEASE" --stats "${W}/stats.tsv" \
    --out-tsv "${O}/${SAMPLE}_ancient.tsv" --out-info "${O}/${SAMPLE}_ancient_info.tsv" \
    --out-txt "${O}/${SAMPLE}_ancient.txt" --out-svg "${O}/${SAMPLE}_ancient_pca.svg"
# The merged genotypes are as large as the panel's (about 1 GB): not kept.
rm -f "${WORK}/merged.geno"
echo ""
echo "Results: ${OUTDIR}/ (the plot: ${SAMPLE}_ancient_pca.svg)"
