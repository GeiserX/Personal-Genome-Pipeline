#!/usr/bin/env bash
# Alignment — minimap2 + samtools sort (FASTQ to sorted BAM)
# Input: paired-end FASTQ files + GRCh38 reference
# Output: sorted BAM + BAI index in $GENOME_DIR/<sample>/aligned/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
R1="${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
MMI="${GENOME_DIR}/reference/GRCh38.mmi"
OUTPUT_DIR="${SAMPLE_DIR}/aligned"

echo "=== Alignment: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Reference: ${REF}"

# Validate inputs
for f in "$R1" "$R2" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Step 1: Build minimap2 index (one-time, ~30 min)
if [ ! -f "$MMI" ]; then
  echo "Building minimap2 index (one-time, ~30 min)..."
  docker run --rm \
    --cpus 8 --memory 16g \
    -v "${GENOME_DIR}:/genome" \
    quay.io/biocontainers/minimap2:2.28--he4a0461_0 \
    minimap2 -d /genome/reference/GRCh38.mmi \
      /genome/reference/Homo_sapiens_assembly38.fasta
fi

# Step 2: Align + sort (1-2 hours for 30X WGS)
# minimap2 runs in its own container, pipes SAM to samtools for sorting.
# The -i flag on the samtools container keeps stdin open for the pipe.
echo "Aligning reads (this takes 1-2 hours for 30X WGS)..."
docker run --rm \
  --cpus "${THREADS}" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/minimap2:2.28--he4a0461_0 \
  minimap2 -t "${THREADS}" -a -x sr \
    /genome/reference/GRCh38.mmi \
    "/genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" \
    "/genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" \
| docker run --rm -i \
  --cpus "${THREADS}" --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.20 \
  samtools sort -@ "${THREADS}" \
    -o "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

# Step 3: Index BAM
echo "Indexing BAM..."
docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.20 \
  samtools index "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

echo "=== Alignment complete ==="
echo "BAM: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
echo "Index: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai"
ls -lh "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" 2>/dev/null || true
