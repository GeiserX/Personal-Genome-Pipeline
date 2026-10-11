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
#
# The BAM is written under a temporary name, indexed and checked with samtools
# quickcheck, and only then renamed, as in step 02: a killed run leaves no
# <sample>_sorted.bam behind for steps 03e and 04c to trust. THREADS (default
# 8) sets the CPUs of minimap2 and of samtools sort.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
PLATFORM=${PLATFORM:?Set PLATFORM to ont or hifi}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/aligned_longread"

# Select minimap2 preset based on platform
case "$PLATFORM" in
  ont)
    MM2_PRESET="map-ont"
    RG_PLATFORM="ONT"
    ;;
  hifi)
    MM2_PRESET="map-hifi"
    RG_PLATFORM="PACBIO"
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

# Resolve GENOME_DIR once (Docker resolves symlinks for bind mounts, so the
# container-relative path must be computed from the real filesystem path).
_resolve() {
  python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$1" 2>/dev/null \
    || readlink -f "$1" 2>/dev/null \
    || echo "$1"
}
REAL_GENOME=$(_resolve "$GENOME_DIR")

# Validate INPUT is physically inside GENOME_DIR (Docker only mounts GENOME_DIR).
# Resolve symlinks so a link under GENOME_DIR pointing outside still fails.
if [ -n "${INPUT_FILE}" ] && [ -f "${INPUT_FILE}" ]; then
  REAL_INPUT=$(_resolve "$INPUT_FILE")
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

BAM="${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
TMP_BAM="${OUTPUT_DIR}/${SAMPLE}_sorted.tmp.bam"
SORT_TMP="${OUTPUT_DIR}/${SAMPLE}.sort_tmp"
# Whatever this run leaves half written goes when it exits, finished or not.
cleanup() { rm -rf "$TMP_BAM" "${TMP_BAM}.bai" "$SORT_TMP"; }
trap cleanup EXIT

# Compute container-relative input path from resolved paths.
# Docker resolves symlinks on bind mounts, so /genome/ maps to REAL_GENOME.
REAL_INPUT=${REAL_INPUT:-$(_resolve "$INPUT_FILE")}
INPUT_RELPATH="${REAL_INPUT#"${REAL_GENOME}/"}"

# Align + sort
# Long-read minimap2 does NOT use a pre-built .mmi index — the preset-specific index
# differs from the short-read one. minimap2 builds it on the fly from the FASTA.
# $1 = reads path inside the container, or - to read FASTQ from stdin.
# Any further arguments are extra minimap2 options.
# samtools sort spills to SORT_TMP (aligned_longread/<sample>.sort_tmp), a
# directory of its own that the trap removes, not to the container's own disk;
# -m is per thread, so the container gets THREADS + 4 GB.
SORT_MEM_GB=$((THREADS + 4))
_align_and_sort() {
  local reads="$1"
  shift
  # Only the uBAM path pipes reads in; -i on file input would make docker
  # read this script's stdin.
  local stdin_flag=()
  [ "$reads" = "-" ] && stdin_flag=(-i)
  run_in ${stdin_flag[@]+"${stdin_flag[@]}"} \
    --cpus "${THREADS}" --memory 16g \
    "$MINIMAP2_IMAGE" \
    minimap2 -t "${THREADS}" -a -x "${MM2_PRESET}" \
      --MD -Y \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${RG_PLATFORM}\tLB:${SAMPLE}" \
      "$@" \
      "${REF_FASTA_C}" \
      "$reads" \
  | run_in -i \
    --cpus "${THREADS}" --memory "${SORT_MEM_GB}g" \
    "$SAMTOOLS_IMAGE" \
    samtools sort -@ "${THREADS}" -m 1G -T "$(cpath "$SORT_TMP")/sort" \
      -o "$(cpath "$TMP_BAM")"
}

mkdir -p "$SORT_TMP"
echo "[1/2] Aligning long reads with minimap2 (preset: ${MM2_PRESET})..."
echo "       This takes 1-3 hours for 30X long-read WGS."
if [[ "$INPUT_RELPATH" == *.bam ]]; then
  # minimap2 reads FASTA/FASTQ only. An unaligned BAM (the usual PacBio HiFi
  # delivery) is streamed through samtools fastq; -T MM,ML puts the base
  # modification tags in the read comment and minimap2 -y copies them back.
  echo "       Unaligned BAM input: converting to FASTQ on the fly (MM/ML tags kept)."
  run_in \
    --cpus 2 --memory 4g \
    "$SAMTOOLS_IMAGE" \
    samtools fastq -T MM,ML "/genome/${INPUT_RELPATH}" \
  | _align_and_sort - -y
else
  _align_and_sort "/genome/${INPUT_RELPATH}"
fi

# Index and check, then rename. The old index goes first, so an index never
# sits next to a BAM it was not built from.
echo "[2/2] Indexing and checking BAM..."
run_in \
  --cpus "${THREADS}" --memory 2g \
  "$SAMTOOLS_IMAGE" \
  samtools index -@ "${THREADS}" "$(cpath "$TMP_BAM")"
run_in \
  --cpus 1 --memory 1g \
  "$SAMTOOLS_IMAGE" \
  samtools quickcheck -v "$(cpath "$TMP_BAM")"
rm -f "${BAM}.bai"
mv -f "$TMP_BAM" "$BAM"
mv -f "${TMP_BAM}.bai" "${BAM}.bai"

echo "=== Long-read Alignment complete ==="
echo "BAM: ${BAM}"
echo "Index: ${BAM}.bai"
ls -lh "$BAM" 2>/dev/null || true
echo ""
echo "Next steps:"
echo "  - Variant calling: PLATFORM=${PLATFORM} ./scripts/03e-clair3.sh ${SAMPLE} [male|female]"
echo "  - SV calling:      ALIGN_DIR=aligned_longread ./scripts/04c-sniffles2.sh ${SAMPLE}"
echo "  - Or DeepVariant: ALIGN_DIR=aligned_longread MODEL_TYPE=$([ "$PLATFORM" = "ont" ] && echo "ONT_R104" || echo "PACBIO") ./scripts/03-deepvariant.sh ${SAMPLE} [male|female]"
