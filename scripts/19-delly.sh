#!/usr/bin/env bash
# Delly — Structural variant caller (paired-end + split-read + read-depth)
# Input: Sorted BAM + reference FASTA
# Output: SV VCF (DEL, DUP, INV, BND, INS)
# Runtime: ~2-4 hours per 30X genome
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/delly"

echo "=== Delly SV calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

DELLY_IMAGE="${DELLY_IMAGE}"
BCFTOOLS_IMAGE="${BCFTOOLS_IMAGE}"

echo "[1/3] Calling structural variants..."
run_in  --cpus 4 --memory 8g \
  "$DELLY_IMAGE" \
  delly call \
    -g "${REF_FASTA_C}" \
    -o "/genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf" \
    "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

echo "[2/3] Converting BCF to VCF..."
run_in  "$BCFTOOLS_IMAGE" \
  bcftools view \
    "/genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf" \
    -Oz -o "/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz"

echo "[3/3] Indexing VCF..."
run_in  "$BCFTOOLS_IMAGE" \
  bcftools index -t \
    "/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz"

echo "=== Delly complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz"
echo ""
echo "View PASS variants only:"
echo "  bcftools view -f PASS ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz"
echo ""
echo "Count by SV type:"
echo "  bcftools query -f '%INFO/SVTYPE\n' ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz | sort | uniq -c"
