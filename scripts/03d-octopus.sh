#!/usr/bin/env bash
# Octopus — Haplotype-aware Bayesian variant caller
# Input: sorted BAM + GRCh38 reference FASTA
# Output: VCF in $GENOME_DIR/<sample>/vcf_octopus/
#
# Octopus uses a haplotype-aware, reference-panel-free Bayesian genotype model.
# Designed as a 5th caller for benchmarking alongside DeepVariant, GATK, FreeBayes,
# and Strelka2, or as a potential FreeBayes replacement (faster, lower memory).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_octopus"

echo "=== Octopus: ${SAMPLE} ==="
echo "BAM: ${BAM}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Skip if output already exists
if [ -f "${OUTPUT_DIR}/${SAMPLE}.vcf.gz" ]; then
  echo "Octopus output already exists, skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

mkdir -p "$OUTPUT_DIR"

# Restrict to specific regions if INTERVALS is set (e.g., INTERVALS=chr22 for testing)
REGION_ARGS=()
if [ -n "${INTERVALS:-}" ]; then
  echo "Restricting to region: ${INTERVALS}"
  REGION_ARGS=(--regions "${INTERVALS}")
fi

# Octopus in germline mode (default)
# Flags:
#   -R           Reference FASTA
#   -I           Input BAM
#   -o           Output VCF
#   --threads    Worker threads
#   --regions    Restrict to regions (optional, for testing)
echo "Running Octopus (this takes 2-4 hours for 30X WGS)..."
docker run --rm --user root \
  --cpus "${THREADS}" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  dancooke/octopus:0.7.4 \
  octopus \
    -R /genome/reference/Homo_sapiens_assembly38.fasta \
    -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -o "/genome/${SAMPLE}/vcf_octopus/${SAMPLE}.vcf.gz" \
    --threads "${THREADS}" \
    "${REGION_ARGS[@]}"

echo "=== Octopus complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
ls -lh "${OUTPUT_DIR}/${SAMPLE}.vcf.gz" 2>/dev/null || true
