#!/usr/bin/env bash
# GATK HaplotypeCaller — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_gatk/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (with .dict and .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_gatk/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
INTERVALS=${INTERVALS:-""}

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
REF_DICT="${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_gatk"

GATK_IMAGE="${GATK_IMAGE}"
BCFTOOLS_IMAGE="${BCFTOOLS_IMAGE}"

echo "=== GATK HaplotypeCaller: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"
if [ -n "$INTERVALS" ]; then
  echo "Intervals: ${INTERVALS}"
fi
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" "$REF_DICT"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Build GATK command
GATK_CMD=(
  gatk HaplotypeCaller
  -R "${REF_FASTA_C}"
  -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  -O "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
  --native-pair-hmm-threads "$THREADS"
  -ERC NONE
)

if [ -n "$INTERVALS" ]; then
  GATK_CMD+=(--intervals "$INTERVALS")
fi

echo "=== [1/3] Running GATK HaplotypeCaller ==="
run_in --cpus "$THREADS" --memory 32g \
  "$GATK_IMAGE" \
  "${GATK_CMD[@]}"

echo "=== [2/3] Indexing VCF with bcftools ==="
run_in --cpus 2 --memory 2g \
  "$BCFTOOLS_IMAGE" \
  bcftools index -ft "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"

echo "=== [3/3] Variant statistics ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "$BCFTOOLS_IMAGE" bcftools stats "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"

echo "=== GATK HaplotypeCaller complete ==="
