#!/usr/bin/env bash
# VEP — Ensembl Variant Effect Predictor
# Full functional annotation: consequence, SIFT, PolyPhen, regulatory, etc.
# Requires: VEP cache (~17GB download, one-time)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
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

# Check if cache exists
if [ ! -d "${CACHE_DIR}/homo_sapiens" ]; then
  echo "VEP cache not found. Installing..."
  echo "Step 1: Download cache (17GB, takes ~10-20 min)"

  # Manual download is more reliable than INSTALL.pl
  mkdir -p "${CACHE_DIR}/tmp"
  wget -c https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz \
    -O "${CACHE_DIR}/tmp/homo_sapiens_vep_112_GRCh38.tar.gz"

  echo "Step 2: Extract cache..."
  cd "$CACHE_DIR" && tar xzf tmp/homo_sapiens_vep_112_GRCh38.tar.gz
  echo "Cache installed at ${CACHE_DIR}/homo_sapiens/"
fi

# Run VEP
docker run --rm \
  --cpus 4 --memory 8g \
  --user root \
  -v "${GENOME_DIR}:/genome" \
  -v "${CACHE_DIR}:/opt/vep/.vep" \
  "${VEP_IMAGE}" \
  vep \
    --input_file "/genome/${SAMPLE}/${VCF_DIR}/${SAMPLE}.vcf.gz" \
    --output_file "/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf" \
    --vcf \
    --cache \
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
