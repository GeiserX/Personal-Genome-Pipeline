#!/usr/bin/env bash
# GATK HaplotypeCaller — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_gatk/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (with .dict and .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_gatk/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
INTERVALS=${INTERVALS:-""}

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
REF_DICT="${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_gatk"

GATK_IMAGE="broadinstitute/gatk:4.6.1.0"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

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
  -R /genome/reference/Homo_sapiens_assembly38.fasta
  -I "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
  -O "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
  --native-pair-hmm-threads "$THREADS"
  -ERC NONE
)

if [ -n "$INTERVALS" ]; then
  GATK_CMD+=(--intervals "$INTERVALS")
fi

echo "=== [1/3] Running GATK HaplotypeCaller ==="
docker run --rm --user root \
  --cpus "$THREADS" --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  "$GATK_IMAGE" \
  "${GATK_CMD[@]}"

echo "=== [2/3] Indexing VCF with bcftools ==="
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools index -ft "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"

echo "=== [3/3] Variant statistics ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(docker run --rm -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" bcftools stats "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"

echo "=== GATK HaplotypeCaller complete ==="
