#!/usr/bin/env bash
# Strelka2 — Alternative germline SV + small variant caller
# Alternative to step 04 (Manta). Outputs to sv_strelka2/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference
# Output: germline variants VCF in $GENOME_DIR/<sample>/sv_strelka2/
# Runtime: ~2-4 hours per 30X genome
#
# NOTE: Strelka2's SomaticEVS scoring model is optimized for BWA-MEM alignments.
# When used with minimap2-aligned BAMs, SNP precision may be reduced due to:
#   - Missing XS (suboptimal alignment score) tags that BWA-MEM produces
#   - Different AS (alignment score) scaling between minimap2 and BWA-MEM
# Indel and SV calls are less affected. Consider this when interpreting results
# from minimap2 alignments.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/sv_strelka2"
THREADS=8

echo "=== Strelka2 Germline Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

STRELKA_IMAGE="quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# Step 1: Configure Strelka2 germline workflow
echo "[1/2] Configuring Strelka2 germline workflow..."
docker run --rm --user root \
  --cpus "$THREADS" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  "$STRELKA_IMAGE" \
  configureStrelkaGermlineWorkflow.py \
    --bam "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    --referenceFasta /genome/reference/Homo_sapiens_assembly38.fasta \
    --runDir "/genome/${SAMPLE}/sv_strelka2"

# Step 2: Run the workflow
echo "[2/2] Running Strelka2 (this takes 2-4 hours for 30X WGS)..."
docker run --rm --user root \
  --cpus "$THREADS" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  "$STRELKA_IMAGE" \
  "/genome/${SAMPLE}/sv_strelka2/runWorkflow.py" \
    -m local \
    -j "$THREADS"

VARIANT_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view -H "/genome/${SAMPLE}/sv_strelka2/results/variants/variants.vcf.gz" | wc -l)

echo "=== Strelka2 complete ==="
echo "Total variants called: ${VARIANT_COUNT}"
echo "Results: ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "View PASS variants only:"
echo "  bcftools view -f PASS ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "Count variant types:"
echo "  bcftools stats ${OUTPUT_DIR}/results/variants/variants.vcf.gz | grep 'number of records'"
