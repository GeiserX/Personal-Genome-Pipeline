#!/usr/bin/env bash
# fastp QC + adapter trimming — pre-alignment quality control
# Input: paired-end FASTQ files from $GENOME_DIR/<sample>/fastq/
# Output: trimmed FASTQs + QC reports in $GENOME_DIR/<sample>/fastq_trimmed/
#
# Performs: adapter auto-detection (Illumina + BGI/MGI built-in), quality
# trimming (Q20 sliding window), polyG tail removal, length filtering (>=36bp),
# and per-read overlap-based PE adapter detection.
#
# Skip with: SKIP_TRIM=true
# The alignment step (02-alignment.sh) auto-detects trimmed FASTQs.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
R1="${SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"
OUTPUT_DIR="${SAMPLE_DIR}/fastq_trimmed"

# Skip if requested
if [ "${SKIP_TRIM:-false}" = "true" ]; then
  echo "=== fastp QC: SKIPPED (SKIP_TRIM=true) ==="
  echo "Alignment will use raw FASTQs."
  exit 0
fi

echo "=== fastp QC: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"

# Validate inputs
for f in "$R1" "$R2"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Skip if trimmed output already exists
if [ -f "${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz" ] && [ -f "${OUTPUT_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
  echo "Trimmed FASTQs already exist in ${OUTPUT_DIR}/, skipping."
  echo "Delete them to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

mkdir -p "$OUTPUT_DIR"

# fastp: adapter trimming + QC
# Flags:
#   --detect_adapter_for_pe  Enable PE overlap-based adapter detection (off by default for PE)
#   --qualified_quality_phred 20  Mark bases below Q20 as unqualified
#   --cut_front / --cut_tail  Sliding window quality trimming from both ends
#   --cut_mean_quality 20    Window mean quality threshold
#   --length_required 36     Discard reads shorter than 36bp after trimming
#   -g / --trim_poly_g       Enable polyG tail trimming (NovaSeq/NextSeq two-color chemistry)
#   -R                       Report title (used by MultiQC for sample naming)
#   -w                       Worker threads (default 3, max 16 effective for I/O-bound work)
echo "Running fastp (adapter trimming + quality filtering)..."
docker run --rm --user root \
  --cpus "${THREADS}" --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/fastp:1.3.1--h43da1c4_0 \
  fastp \
    -i "/genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz" \
    -I "/genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz" \
    -o "/genome/${SAMPLE}/fastq_trimmed/${SAMPLE}_R1.fastq.gz" \
    -O "/genome/${SAMPLE}/fastq_trimmed/${SAMPLE}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    --length_required 36 \
    -g \
    -R "${SAMPLE}" \
    -j "/genome/${SAMPLE}/fastq_trimmed/${SAMPLE}_fastp.json" \
    -h "/genome/${SAMPLE}/fastq_trimmed/${SAMPLE}_fastp.html" \
    -w "${THREADS}"

echo "=== fastp QC complete ==="
echo "Trimmed R1:   ${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz"
echo "Trimmed R2:   ${OUTPUT_DIR}/${SAMPLE}_R2.fastq.gz"
echo "JSON report:  ${OUTPUT_DIR}/${SAMPLE}_fastp.json"
echo "HTML report:  ${OUTPUT_DIR}/${SAMPLE}_fastp.html"
ls -lh "${OUTPUT_DIR}/${SAMPLE}_R1.fastq.gz" 2>/dev/null || true
