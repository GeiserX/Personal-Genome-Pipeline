#!/usr/bin/env bash
# duphold — Annotate structural variants with depth-based quality metrics
# Input: Manta diploidSV.vcf.gz + sorted BAM + reference FASTA
# Output: SV VCF with DHBFC/DHFFC annotations (filter false positives)
# Very fast (~20 minutes)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
MANTA_VCF="${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz"
# Fall back to manta2/ if a second Manta run was used
[ ! -f "$MANTA_VCF" ] && MANTA_VCF="${SAMPLE_DIR}/manta2/results/variants/diploidSV.vcf.gz"
OUTPUT_DIR="${SAMPLE_DIR}/duphold"

echo "=== duphold: ${SAMPLE} ==="
echo "Input SV VCF: ${MANTA_VCF}"
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf"

# Validate inputs
for f in "$MANTA_VCF" "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  brentp/duphold:v0.2.3 \
  duphold \
    -v "/genome/${SAMPLE}/$(echo "$MANTA_VCF" | sed "s|${SAMPLE_DIR}/||")" \
    -b "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    -f /genome/reference/Homo_sapiens_assembly38.fasta \
    -o "/genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf"

echo "=== duphold complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf"
echo ""
echo "Filter high-confidence DELs (DHFFC < 0.7):"
echo "  grep -v '^#' ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf | awk '\$8 ~ /DHFFC=/ && \$5 ~ /DEL/'"
echo ""
echo "Filter high-confidence DUPs (DHBFC > 1.3):"
echo "  grep -v '^#' ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf | awk '\$8 ~ /DHBFC=/ && \$5 ~ /DUP/'"
