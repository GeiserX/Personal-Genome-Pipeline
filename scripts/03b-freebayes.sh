#!/usr/bin/env bash
# FreeBayes — Alternative variant caller (SNPs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_freebayes/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference (.fasta + .fai)
# Output: VCF.gz in $GENOME_DIR/<sample>/vcf_freebayes/
# Runtime: ~9 hours single-threaded for 30X WGS
# Memory: peaks at ~13 GB for full genome; needs 32 GB allocation for safety margin
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
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
FREEBAYES_ARGS=(-f "${REF_FASTA_C}")
if [ -n "$INTERVALS" ]; then
  FREEBAYES_ARGS+=(--region "$INTERVALS")
fi
FREEBAYES_ARGS+=("/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam")

# Through a temporary name: a FreeBayes that fails leaves no raw VCF behind.
atomic_out "${OUTPUT_DIR}/${SAMPLE}_raw.vcf" run_in \
  --cpus 4 --memory 32g \
  "${FREEBAYES_IMAGE}" \
  freebayes "${FREEBAYES_ARGS[@]}"

# Step 2: Sort and compress in one bcftools call (no pipe whose first half can
# fail unseen), with its temporary files in the output directory, then index.
# The raw VCF is removed only after both succeeded.
echo "Sorting and compressing VCF..."
run_in \
  --cpus 4 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bcftools sort -Oz -o "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz" \
    -T "/genome/${SAMPLE}/vcf_freebayes/sort-tmp" \
    "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}_raw.vcf"

echo "Indexing VCF..."
run_in \
  --cpus 1 --memory 1g \
  "${BCFTOOLS_IMAGE}" \
  bcftools index -f -t "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"

# Clean up raw unsorted VCF
rm -f "${OUTPUT_DIR}/${SAMPLE}_raw.vcf"

echo "=== FreeBayes complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "${BCFTOOLS_IMAGE}" bcftools stats "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
echo ""
echo "NOTE: FreeBayes tends to call more variants than DeepVariant (higher sensitivity, more false positives)."
echo "Consider running bcftools filter or vcffilter for quality filtering."
