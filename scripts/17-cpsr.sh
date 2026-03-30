#!/usr/bin/env bash
# CPSR — Cancer Predisposition Sequencing Reporter (ACMG SF v3.2)
# Input: Germline VCF + PCGR data bundle
# Output: HTML report + classified variant TSV
# Requires: 21GB data bundle (download once, reuse for all samples)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
VCF="${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz"
DATA_DIR="${GENOME_DIR}/pcgr_data"
OUTPUT_DIR="${SAMPLE_DIR}/cpsr"

echo "=== CPSR Cancer Predisposition: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "Data bundle: ${DATA_DIR}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}" >&2
  exit 1
fi

if [ ! -d "${DATA_DIR}/data" ]; then
  echo "ERROR: PCGR data bundle not found at ${DATA_DIR}/data/" >&2
  echo "Download and extract it first:" >&2
  echo "  mkdir -p ${GENOME_DIR}/pcgr_data && cd ${GENOME_DIR}/pcgr_data" >&2
  echo "  wget -c http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20220203.tgz" >&2
  echo "  tar xzf pcgr.databundle.grch38.20220203.tgz" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  sigven/pcgr:1.4.1 \
  cpsr \
    --input_vcf "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --pcgr_dir /genome/pcgr_data \
    --output_dir "/genome/${SAMPLE}/cpsr" \
    --genome_assembly grch38 \
    --sample_id "${SAMPLE}" \
    --panel_id 0 \
    --classify_all \
    --force_overwrite

echo "=== CPSR complete ==="
echo "HTML report: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.html"
echo "Variant table: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.snvs_indels.tiers.tsv"
