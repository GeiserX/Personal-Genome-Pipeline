#!/usr/bin/env bash
# BWA-MEM2 — Alternative aligner (FASTQ to sorted BAM)
# Alternative to step 02 (minimap2). Outputs to aligned_bwamem2/ to avoid conflicts.
# Input: paired-end FASTQ files + GRCh38 reference
# Output: sorted BAM + BAI index in $GENOME_DIR/<sample>/aligned_bwamem2/
# Note: BWA-MEM2 produces XS (suboptimal alignment score) tags that some callers
#       (especially Strelka2) depend on. minimap2 does not produce XS tags.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
R1="${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
BWA_INDEX="${REF}.bwt.2bit.64"
OUTPUT_DIR="${SAMPLE_DIR}/aligned_bwamem2"
SAM_TEMP="${OUTPUT_DIR}/${SAMPLE}.sam"

echo "=== BWA-MEM2 Alignment: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
for f in "$R1" "$R2" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Step 1: Build BWA-MEM2 index if not present (one-time, ~1 hour)
if [ ! -f "$BWA_INDEX" ]; then
  echo "=== Building BWA-MEM2 index (one-time, ~1 hour) ==="
  docker run --rm \
    --user root \
    --cpus 8 --memory 24g \
    -v "${GENOME_DIR}:/genome" \
    quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5 \
    bwa-mem2 index /genome/reference/Homo_sapiens_assembly38.fasta
  echo "BWA-MEM2 index built."
else
  echo "BWA-MEM2 index found, skipping build."
fi

# Step 2: Align with BWA-MEM2 (4-8 hours for 30X WGS, writes SAM to temp file)
echo "=== Aligning reads with BWA-MEM2 (this takes 4-8 hours for 30X WGS) ==="
docker run --rm \
  --user root \
  --cpus 8 --memory 24g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/bwa-mem2:2.2.1--hd03093a_5 \
  bwa-mem2 mem -t "${THREADS}" \
    /genome/reference/Homo_sapiens_assembly38.fasta \
    "/genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" \
    "/genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" \
    -o "/genome/${SAMPLE}/aligned_bwamem2/${SAMPLE}.sam"

# Step 3: Sort + compress with samtools
echo "=== Sorting and compressing BAM ==="
docker run --rm \
  --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.20 \
  samtools sort -@ 2 \
    -o "/genome/${SAMPLE}/aligned_bwamem2/${SAMPLE}_sorted.bam" \
    "/genome/${SAMPLE}/aligned_bwamem2/${SAMPLE}.sam"

# Remove temporary SAM file (can be 200+ GB for 30X WGS)
echo "Removing temporary SAM file..."
rm -f "$SAM_TEMP"

# Step 4: Index BAM
echo "=== Indexing BAM ==="
docker run --rm \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.20 \
  samtools index "/genome/${SAMPLE}/aligned_bwamem2/${SAMPLE}_sorted.bam"

echo "=== BWA-MEM2 Alignment complete ==="
echo "BAM: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
echo "Index: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai"
ls -lh "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" 2>/dev/null || true
