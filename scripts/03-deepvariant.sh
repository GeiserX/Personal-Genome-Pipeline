#!/usr/bin/env bash
# DeepVariant — Variant calling (BAM to VCF)
# Input: sorted BAM + GRCh38 reference
# Output: VCF.gz with SNPs and small indels (~5.5M variants per 30X WGS)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/vcf"

# Select DeepVariant model type: WGS (default), WES, or PACBIO/ONT_R104
# WES uses a model trained on exome depth profiles and capture boundaries.
MODEL_TYPE=${MODEL_TYPE:-WGS}
case "$MODEL_TYPE" in
  WGS|WES|PACBIO|ONT_R104) ;;
  *) echo "ERROR: MODEL_TYPE must be WGS, WES, PACBIO, or ONT_R104, got '${MODEL_TYPE}'" >&2; exit 1 ;;
esac

echo "=== DeepVariant: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Model type: ${MODEL_TYPE}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 8 --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type="${MODEL_TYPE}" \
    --ref="/genome/reference/Homo_sapiens_assembly38.fasta" \
    --reads="/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --output_vcf="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --num_shards=8

echo "=== DeepVariant complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
