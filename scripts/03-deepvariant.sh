#!/usr/bin/env bash
# DeepVariant — Variant calling (BAM to VCF)
# Input: sorted BAM + GRCh38 reference
# Output: VCF.gz with SNPs and small indels (~5.5M variants per 30X WGS)
# Optional: INTERVALS="chr20:10000001-10500000 chr22:1-50818468" calls only
#   those regions (space-separated region literals, passed to --regions).
#   Unset means the whole genome.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
INTERVALS=${INTERVALS:-}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
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
if [ -n "$INTERVALS" ]; then
  echo "Regions: ${INTERVALS}"
fi
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

DV_ARGS=(
  --model_type="${MODEL_TYPE}"
  --ref="${REF_FASTA_C}"
  --reads="/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  --output_vcf="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
  --sample_name="${SAMPLE}"
  --num_shards=8
)
if [ -n "$INTERVALS" ]; then
  DV_ARGS+=(--regions "$INTERVALS")
fi

run_in \
  --cpus 8 --memory 32g \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant "${DV_ARGS[@]}"

echo "=== DeepVariant complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "${BCFTOOLS_IMAGE}" bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
