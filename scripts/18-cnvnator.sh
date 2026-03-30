#!/usr/bin/env bash
# CNVnator — Depth-based CNV detection (orthogonal to Manta)
# Input: Sorted BAM + reference FASTA
# Output: CNV calls (deletions + duplications)
# Runtime: ~2-4 hours per 30X genome
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/cnvnator"
BIN_SIZE=1000  # 1000bp bins for 30X WGS

echo "=== CNVnator: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Bin size: ${BIN_SIZE} bp"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

IMAGE="quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11"

echo "[1/5] Extracting read mapping..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$IMAGE" \
  cnvnator \
    -root "/genome/${SAMPLE}/cnvnator/${SAMPLE}.root" \
    -tree "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

echo "[2/5] Generating read-depth histogram..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$IMAGE" \
  cnvnator \
    -root "/genome/${SAMPLE}/cnvnator/${SAMPLE}.root" \
    -his "$BIN_SIZE" \
    -fasta /genome/reference/Homo_sapiens_assembly38.fasta

echo "[3/5] Computing statistics..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$IMAGE" \
  cnvnator \
    -root "/genome/${SAMPLE}/cnvnator/${SAMPLE}.root" \
    -stat "$BIN_SIZE"

echo "[4/5] Partitioning..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$IMAGE" \
  cnvnator \
    -root "/genome/${SAMPLE}/cnvnator/${SAMPLE}.root" \
    -partition "$BIN_SIZE"

echo "[5/5] Calling CNVs..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$IMAGE" \
  cnvnator \
    -root "/genome/${SAMPLE}/cnvnator/${SAMPLE}.root" \
    -call "$BIN_SIZE" \
  > "${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"

CNV_COUNT=$(wc -l < "${OUTPUT_DIR}/${SAMPLE}_cnvs.txt")
echo "=== CNVnator complete ==="
echo "Total CNVs called: ${CNV_COUNT}"
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"
echo ""
echo "Filter significant CNVs (e-value < 0.01):"
echo "  awk '\$5 < 0.01' ${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"
