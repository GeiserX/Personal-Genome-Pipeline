#!/usr/bin/env bash
# MToolBox — Mitochondrial DNA analysis (heteroplasmy + disease variants)
# Input: Sorted BAM
# Output: Prioritized variants, heteroplasmy levels, haplogroup
# Runtime: ~15-30 minutes
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/mtoolbox"

echo "=== MToolBox Mitochondrial Analysis: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${OUTPUT_DIR}/output"

SAMTOOLS_IMAGE="staphb/samtools:1.21"
MTOOLBOX_IMAGE="robertopreste/mtoolbox:latest"

echo "[1/2] Extracting chrM reads to FASTQ..."
docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "$SAMTOOLS_IMAGE" \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM | \
    samtools sort -n - | \
    samtools fastq \
      -1 /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_R1.fastq.gz \
      -2 /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_R2.fastq.gz \
      -s /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_singleton.fastq.gz -
  "

echo "[2/2] Running MToolBox analysis..."
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${OUTPUT_DIR}:/input" \
  -v "${OUTPUT_DIR}/output:/output" \
  "$MTOOLBOX_IMAGE" \
  MToolBox.sh \
    -i "/input/${SAMPLE}_chrM_R1.fastq.gz" \
    -I "/input/${SAMPLE}_chrM_R2.fastq.gz" \
    -o /output \
    -m "-t 4"

echo "=== MToolBox complete ==="
echo "Prioritized variants: ${OUTPUT_DIR}/output/prioritized_variants.txt"
echo "Full annotation: ${OUTPUT_DIR}/output/annotation.csv"
echo "Haplogroup: ${OUTPUT_DIR}/output/mt_classification_best_results.csv"
echo ""
echo "View top pathogenic variants:"
echo "  head -20 ${OUTPUT_DIR}/output/prioritized_variants.txt"
