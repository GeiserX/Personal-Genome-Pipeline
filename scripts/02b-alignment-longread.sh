#!/usr/bin/env bash
# Long-read Alignment — minimap2 + samtools sort (FASTQ/BAM to sorted BAM)
# Supports Oxford Nanopore (ONT) and PacBio HiFi long-read data.
# Alternative to step 02 (short-read alignment). Outputs to aligned_longread/ to avoid conflicts.
# Input: single FASTQ (.fastq.gz) or unaligned BAM (.bam) + GRCh38 reference
# Output: sorted BAM + BAI index in $GENOME_DIR/<sample>/aligned_longread/
#
# Long reads are single-end (no R1/R2 pairs). Set PLATFORM to select the minimap2 preset:
#   PLATFORM=ont   -> Oxford Nanopore (minimap2 preset: map-ont)
#   PLATFORM=hifi  -> PacBio HiFi/CCS (minimap2 preset: map-hifi)
#
# Runtime: ~1-3 hours for 30X long-read WGS depending on read length and throughput.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
PLATFORM=${PLATFORM:?Set PLATFORM to ont or hifi}
THREADS=${THREADS:-8}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/aligned_longread"

# Select minimap2 preset based on platform
case "$PLATFORM" in
  ont)
    MM2_PRESET="map-ont"
    ;;
  hifi)
    MM2_PRESET="map-hifi"
    ;;
  *)
    echo "ERROR: PLATFORM must be 'ont' or 'hifi', got '${PLATFORM}'" >&2
    exit 1
    ;;
esac

# Find input file: look for FASTQ first, then unaligned BAM
# Long-read data is typically a single file (not paired-end)
INPUT_FILE=""
for candidate in \
  "${SAMPLE_DIR}/fastq/${SAMPLE}.fastq.gz" \
  "${SAMPLE_DIR}/fastq/${SAMPLE}_lr.fastq.gz" \
  "${SAMPLE_DIR}/fastq/${SAMPLE}.fq.gz" \
  "${SAMPLE_DIR}/fastq/${SAMPLE}.bam"; do
  if [ -f "$candidate" ]; then
    INPUT_FILE="$candidate"
    break
  fi
done

# Allow explicit override via INPUT env var
INPUT_FILE="${INPUT:-${INPUT_FILE}}"

# Validate INPUT is physically inside GENOME_DIR (Docker only mounts GENOME_DIR).
# Resolve symlinks so a link under GENOME_DIR pointing outside still fails.
if [ -n "${INPUT_FILE}" ] && [ -f "${INPUT_FILE}" ]; then
  REAL_INPUT=$(python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$INPUT_FILE" 2>/dev/null \
    || readlink -f "$INPUT_FILE" 2>/dev/null \
    || echo "$INPUT_FILE")
  REAL_GENOME=$(python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$GENOME_DIR" 2>/dev/null \
    || readlink -f "$GENOME_DIR" 2>/dev/null \
    || echo "$GENOME_DIR")
  case "$REAL_INPUT" in
    "${REAL_GENOME}/"*)
      ;;
    *)
      echo "ERROR: INPUT path must be physically inside GENOME_DIR (${GENOME_DIR})." >&2
      echo "  Resolved INPUT: ${REAL_INPUT}" >&2
      echo "  Resolved GENOME_DIR: ${REAL_GENOME}" >&2
      echo "  The Docker container only mounts GENOME_DIR. Copy or move your file:" >&2
      echo "  cp \"${INPUT_FILE}\" \"${SAMPLE_DIR}/fastq/\"" >&2
      exit 1
      ;;
  esac
fi

if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
  echo "ERROR: No long-read input file found." >&2
  echo "Looked for:" >&2
  echo "  ${SAMPLE_DIR}/fastq/${SAMPLE}.fastq.gz" >&2
  echo "  ${SAMPLE_DIR}/fastq/${SAMPLE}_lr.fastq.gz" >&2
  echo "  ${SAMPLE_DIR}/fastq/${SAMPLE}.fq.gz" >&2
  echo "  ${SAMPLE_DIR}/fastq/${SAMPLE}.bam" >&2
  echo "Set INPUT=/path/to/reads to override." >&2
  exit 1
fi

echo "=== Long-read Alignment: ${SAMPLE} ==="
echo "Platform: ${PLATFORM} (minimap2 preset: ${MM2_PRESET})"
echo "Input: ${INPUT_FILE}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/"
echo "Threads: ${THREADS}"

# Validate reference exists
if [ ! -f "$REF" ]; then
  echo "ERROR: Reference not found: ${REF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

MINIMAP2_IMAGE="quay.io/biocontainers/minimap2:2.28--he4a0461_0"
SAMTOOLS_IMAGE="staphb/samtools:1.20"

# Compute container-relative input path
INPUT_RELPATH="${INPUT_FILE#"${GENOME_DIR}/"}"

# Align + sort
# Long-read minimap2 does NOT use a pre-built .mmi index — the preset-specific index
# differs from the short-read one. minimap2 builds it on the fly from the FASTA.
echo "[1/2] Aligning long reads with minimap2 (preset: ${MM2_PRESET})..."
echo "       This takes 1-3 hours for 30X long-read WGS."
docker run --rm \
  --cpus "${THREADS}" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  "$MINIMAP2_IMAGE" \
  minimap2 -t "${THREADS}" -a -x "${MM2_PRESET}" \
    --MD -Y \
    /genome/reference/Homo_sapiens_assembly38.fasta \
    "/genome/${INPUT_RELPATH}" \
| docker run --rm -i \
  --cpus "${THREADS}" --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$SAMTOOLS_IMAGE" \
  samtools sort -@ 4 -m 1G \
    -o "/genome/${SAMPLE}/aligned_longread/${SAMPLE}_sorted.bam"

# Index BAM
echo "[2/2] Indexing BAM..."
docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "$SAMTOOLS_IMAGE" \
  samtools index "/genome/${SAMPLE}/aligned_longread/${SAMPLE}_sorted.bam"

echo "=== Long-read Alignment complete ==="
echo "BAM: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
echo "Index: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai"
ls -lh "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" 2>/dev/null || true
echo ""
echo "Next steps:"
echo "  - Variant calling: PLATFORM=${PLATFORM} ./scripts/03e-clair3.sh ${SAMPLE}"
echo "  - SV calling:      ALIGN_DIR=aligned_longread ./scripts/04c-sniffles2.sh ${SAMPLE}"
echo "  - Or use DeepVariant with --model_type=$([ "$PLATFORM" = "ont" ] && echo "ONT_R104" || echo "PACBIO")"
