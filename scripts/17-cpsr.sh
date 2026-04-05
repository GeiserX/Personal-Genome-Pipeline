#!/usr/bin/env bash
# CPSR — Cancer Predisposition Sequencing Reporter (ACMG SF v3.2)
# Input: Germline VCF + PCGR 2.x data bundle + VEP cache
# Output: HTML report + classified variant TSV
# Requires: ~5GB ref data bundle + VEP cache (download once, reuse for all samples)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
VCF="${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz"
VEP_DIR="${GENOME_DIR}/vep_cache"
REFDATA_DIR="${GENOME_DIR}/pcgr_data/20250314"
OUTPUT_DIR="${SAMPLE_DIR}/cpsr"

echo "=== CPSR Cancer Predisposition: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "VEP cache: ${VEP_DIR}"
echo "Ref data bundle: ${REFDATA_DIR}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}" >&2
  exit 1
fi

if [ ! -d "${VEP_DIR}" ]; then
  echo "ERROR: VEP cache not found at ${VEP_DIR}/" >&2
  echo "Download and extract it first:" >&2
  echo "  mkdir -p ${VEP_DIR}" >&2
  echo "  wget -c -P ${VEP_DIR} https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz" >&2
  echo "  tar xzf ${VEP_DIR}/homo_sapiens_vep_112_GRCh38.tar.gz -C ${VEP_DIR}" >&2
  exit 1
fi

if [ ! -d "${REFDATA_DIR}/data" ]; then
  echo "ERROR: PCGR 2.x ref data bundle not found at ${REFDATA_DIR}/data/" >&2
  echo "Download and extract it first:" >&2
  echo "  cd ${GENOME_DIR}/pcgr_data" >&2
  echo "  wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz" >&2
  echo "  tar xzf pcgr_ref_data.20250314.grch38.tgz" >&2
  echo "  mkdir -p 20250314 && mv data/ 20250314/" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${VEP_DIR}:/mnt/.vep" \
  -v "${REFDATA_DIR}:/mnt/bundle" \
  -v "${SAMPLE_DIR}/vcf:/mnt/inputs" \
  -v "${SAMPLE_DIR}/cpsr:/mnt/outputs" \
  sigven/pcgr:2.2.5 \
  cpsr \
    --input_vcf "/mnt/inputs/${SAMPLE}.vcf.gz" \
    --vep_dir /mnt/.vep \
    --refdata_dir /mnt/bundle \
    --output_dir /mnt/outputs \
    --genome_assembly grch38 \
    --sample_id "${SAMPLE}" \
    --panel_id 0 \
    --classify_all \
    --force_overwrite

echo "=== CPSR complete ==="
echo "HTML report: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.html"
echo "Variant table: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.snvs_indels.tiers.tsv"
