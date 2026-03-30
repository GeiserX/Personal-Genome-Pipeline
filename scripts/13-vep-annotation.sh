#!/usr/bin/env bash
# VEP — Ensembl Variant Effect Predictor
# Full functional annotation: consequence, SIFT, PolyPhen, regulatory, etc.
# Requires: VEP cache (~17GB download, one-time)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}
VCF="${GENOMA_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
CACHE_DIR="${GENOMA_DIR}/vep_cache"
OUTPUT_DIR="${GENOMA_DIR}/${SAMPLE}/vep"

echo "=== VEP Annotation: ${SAMPLE} ==="
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
  -v "${GENOMA_DIR}:/genoma" \
  -v "${CACHE_DIR}:/opt/vep/.vep" \
  ensemblorg/ensembl-vep:release_112.0 \
  vep \
    --input_file "/genoma/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --output_file "/genoma/${SAMPLE}/vep/${SAMPLE}_vep.vcf" \
    --vcf \
    --cache \
    --offline \
    --assembly GRCh38 \
    --sift b \
    --polyphen b \
    --regulatory \
    --symbol \
    --canonical \
    --biotype \
    --af_gnomade \
    --max_af \
    --force_overwrite \
    --fork 4

echo "=== VEP complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_vep.vcf"
echo ""
echo "Filter HIGH impact variants:"
echo "  grep 'HIGH' ${OUTPUT_DIR}/${SAMPLE}_vep.vcf | head"
