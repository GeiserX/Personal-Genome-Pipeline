#!/usr/bin/env bash
# FreeBayes — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_freebayes/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (.fasta + .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_freebayes/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_freebayes"
INTERVALS=${INTERVALS:-""}

echo "=== FreeBayes: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
if [ -n "$INTERVALS" ]; then
  echo "Region: ${INTERVALS}"
fi

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Step 1: Run FreeBayes (single-threaded, outputs unsorted VCF)
echo "Running FreeBayes (single-threaded, this may take several hours for 30X WGS)..."
FREEBAYES_ARGS=(-f /genome/reference/Homo_sapiens_assembly38.fasta)
if [ -n "$INTERVALS" ]; then
  FREEBAYES_ARGS+=(--region "$INTERVALS")
fi
FREEBAYES_ARGS+=("/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam")

docker run --rm \
  --cpus 4 --memory 16g \
  --user root \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/freebayes:1.3.6--hbfe0e7f_2 \
  freebayes "${FREEBAYES_ARGS[@]}" \
  > "${OUTPUT_DIR}/${SAMPLE}_raw.vcf"

# Step 2: Sort, compress, and index with bcftools
echo "Sorting and compressing VCF..."
docker run --rm \
  --cpus 4 --memory 4g \
  --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools sort /genome/${SAMPLE}/vcf_freebayes/${SAMPLE}_raw.vcf \
    | bcftools view -Oz -o /genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"

echo "Indexing VCF..."
docker run --rm \
  --cpus 1 --memory 1g \
  --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"

# Clean up raw unsorted VCF
rm -f "${OUTPUT_DIR}/${SAMPLE}_raw.vcf"

echo "=== FreeBayes complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 bcftools stats "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
echo ""
echo "NOTE: FreeBayes tends to call more variants than DeepVariant (higher sensitivity, more false positives)."
echo "Consider running bcftools filter or vcffilter for quality filtering."
