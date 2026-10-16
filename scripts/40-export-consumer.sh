#!/usr/bin/env bash
# 40-export-consumer.sh — a consumer raw-data file (23andMe v5 layout) from the gVCF
# Usage: ./scripts/40-export-consumer.sh <sample_name> [male|female]
#
# Opt-in. Genealogy and health services take the raw-data file a chip
# company gives its customers. bin/export_consumer.py genotypes every site of
# the 23andMe v5 chip list (the site table `setup.sh --consumer-sites`
# builds) from step 03's gVCF: a site inside a reference block is
# homozygous reference, a site with no call is '--'. It writes the 23andMe
# layout and the same sites in the AncestryDNA layout, on GRCh37 positions,
# reads both back with its strict reader, and zips them.
#
# Sex: with `male`, chrX and chrY sites outside the pseudoautosomal regions
#   (assets/par_grch38.bed) get one letter, as on a chip file, and a
#   heterozygous call there is '--'; with `female`, chrY sites are '--'.
#   Without it each site keeps the ploidy of its call. chrM sites always get
#   one letter ('--' for a heteroplasmic call).
#
# The step writes files and sends nothing anywhere. Uploading one to a
# service shares your genotypes with it: docs/40-export-consumer.md says what
# that means.
#
# Input:  vcf/<sample>.g.vcf.gz (step 03)
# Output: consumer_export/<sample>_23andme_v5.zip, <sample>_ancestrydna_layout.zip
#         and <sample>_consumer_export.tsv (what was written, for the report)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
SEX=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE
case "$SEX" in
  ''|male|female) ;;
  *) echo "ERROR: sex must be 'male' or 'female' (or left out), got '${SEX}'." >&2; exit 2 ;;
esac

S="${GENOME_DIR}/${SAMPLE}"
GVCF="${S}/vcf/${SAMPLE}.g.vcf.gz"
SITES="${GENOME_DIR}/reference/consumer/23andme_v5_sites_dbsnp${CONSUMER_DBSNP_BUILD}.tsv.gz"
OUTDIR="${S}/consumer_export"

echo "=== Raw-data export (23andMe v5 layout): ${SAMPLE} ==="
if [ ! -s "$SITES" ]; then
  echo "ERROR: the site table is not built (${SITES}): run scripts/setup.sh --consumer-sites ${GENOME_DIR}" >&2
  exit 1
fi
if [ ! -s "$GVCF" ]; then
  echo "ERROR: no gVCF for ${SAMPLE} (${GVCF}): run scripts/03-deepvariant.sh ${SAMPLE} first." >&2
  echo "  A VCF alone has no reference calls, so most chip sites would be '--'." >&2
  exit 1
fi
echo "Input:  ${GVCF}"
echo "Sites:  ${SITES}"
echo "Sex:    ${SEX:-not given (each site keeps the ploidy of its call)}"

mkdir -p "$OUTDIR"
rm -f "${OUTDIR}/${SAMPLE}_23andme_v5.zip" "${OUTDIR}/${SAMPLE}_ancestrydna_layout.zip" "${OUTDIR}/${SAMPLE}_consumer_export.tsv"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 4g "$PYTHON_IMAGE" \
  python3 /pgp-bin/export_consumer.py export \
    --gvcf "$(cpath "$GVCF")" --sites "$(cpath "$SITES")" --sample "$SAMPLE" --sex "$SEX" \
    --par /pgp-assets/par_grch38.bed --out-dir "$(cpath "$OUTDIR")"

echo ""
echo "Results: ${OUTDIR}/${SAMPLE}_23andme_v5.zip"
echo "         ${OUTDIR}/${SAMPLE}_ancestrydna_layout.zip"
echo "Uploading a file to a service shares these genotypes with it; docs/40-export-consumer.md says what to weigh first."
