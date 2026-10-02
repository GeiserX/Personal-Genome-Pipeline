#!/usr/bin/env bash
# Mitochondrial analysis — Heteroplasmy detection + variant calling with GATK Mutect2
# Input: Sorted BAM
# Output: Mitochondrial VCF with heteroplasmy fractions
# Runtime: ~15-30 minutes
# Note: Originally planned for MToolBox, but no working Docker image exists.
#       GATK Mutect2 in mitochondrial mode is the standard alternative.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/mito"

echo "=== Mitochondrial Analysis (GATK Mutect2): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"


echo "[1/4] Extracting chrM reads..."
run_in --cpus 2 --memory 4g \
  "$SAMTOOLS_IMAGE" \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM \
      > /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam && \
    samtools index /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam
  "

echo "[2/4] Checking sequence dictionary..."
if [ ! -f "${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict" ]; then
  echo "  Creating sequence dictionary..."
  run_in --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk CreateSequenceDictionary \
      -R "${REF_FASTA_C}" \
      -O /genome/reference/Homo_sapiens_assembly38.dict
else
  echo "  Sequence dictionary already exists, skipping."
fi

echo "[3/4] Running Mutect2 in mitochondrial mode..."
run_in --cpus 4 --memory 8g \
  "$GATK_IMAGE" \
  gatk Mutect2 \
    -R "${REF_FASTA_C}" \
    -I "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" \
    -L chrM \
    --mitochondria-mode \
    --max-mnp-distance 0 \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz"

echo "[4/4] Filtering variants..."
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk FilterMutectCalls \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz" \
    --mitochondria-mode \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"

echo "=== Mitochondrial analysis complete ==="
echo "Raw calls: ${OUTPUT_DIR}/${SAMPLE}_chrM_mutect2.vcf.gz"
echo "Filtered: ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz"
echo ""
echo "View heteroplasmic variants (AF < 0.95):"
echo "  bcftools query -f '%POS\t%REF\t%ALT\t[%AF]\n' ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz | awk '\$4 < 0.95'"
