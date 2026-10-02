#!/usr/bin/env bash
# VEP — Ensembl Variant Effect Predictor
# Full functional annotation: consequence, SIFT, PolyPhen, regulatory, etc.
# Requires: VEP cache (~26 GB download, one-time)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF_DIR=${VCF_DIR:-vcf}
VCF="${GENOME_DIR}/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz"
CACHE_DIR="${GENOME_DIR}/vep_cache"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/vep"

echo "=== VEP Annotation: ${SAMPLE} ==="

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# The cache must be the release of VEP_IMAGE (VEP_CACHE_RELEASE in
# versions.env). Another release in the same directory, such as the one CPSR
# uses, does not count.
if [ ! -f "${CACHE_DIR}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  echo "VEP ${VEP_CACHE_RELEASE} cache not found in ${CACHE_DIR}. Installing (26 GB download)..."
  install_vep_cache "$CACHE_DIR" "$VEP_CACHE_RELEASE"
  echo "Cache installed at ${CACHE_DIR}/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/"
fi

# Run VEP. The cache stays writable as before: VEP builds an index for a
# FASTA it finds there without one, and no CI run shows it never writes.
run_in \
  --cpus 4 --memory 8g \
  -v "${CACHE_DIR}:/opt/vep/.vep" \
  "${VEP_IMAGE}" \
  vep \
    --input_file "/genome/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz" \
    --output_file "/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf" \
    --vcf \
    --cache \
    --cache_version "${VEP_CACHE_RELEASE}" \
    --dir_cache /opt/vep/.vep \
    --offline \
    --assembly GRCh38 \
    --everything \
    --force_overwrite \
    --fork 4

echo "=== VEP complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_vep.vcf"
echo ""
echo "Filter HIGH impact variants:"
echo "  grep 'HIGH' ${OUTPUT_DIR}/${SAMPLE}_vep.vcf | head"
