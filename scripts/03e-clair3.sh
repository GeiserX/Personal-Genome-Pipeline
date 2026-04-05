#!/usr/bin/env bash
# Clair3 — Long-read variant calling (BAM to VCF)
# Alternative to step 03 (DeepVariant) for long-read data.
# Supports Oxford Nanopore (ONT) and PacBio HiFi platforms.
# Input: sorted BAM from long-read alignment + GRCh38 reference
# Output: VCF.gz with SNPs and small indels in $GENOME_DIR/<sample>/vcf_clair3/
#
# Set PLATFORM to select the appropriate model:
#   PLATFORM=ont   -> ONT R10.4.1 model (r1041_e82_400bps_sup_v500)
#   PLATFORM=hifi  -> PacBio HiFi/Revio model (hifi_revio)
#
# Runtime: ~2-4 hours for 30X long-read WGS
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
PLATFORM=${PLATFORM:?Set PLATFORM to ont or hifi}
THREADS=${THREADS:-8}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned_longread}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_clair3"

CLAIR3_IMAGE="hkubal/clair3:v2.0.0"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# Select model path based on platform
case "$PLATFORM" in
  ont)
    MODEL_PATH="/opt/models/r1041_e82_400bps_sup_v500"
    ;;
  hifi)
    MODEL_PATH="/opt/models/hifi_revio"
    ;;
  *)
    echo "ERROR: PLATFORM must be 'ont' or 'hifi', got '${PLATFORM}'" >&2
    exit 1
    ;;
esac

echo "=== Clair3 Variant Calling: ${SAMPLE} ==="
echo "Platform: ${PLATFORM}"
echo "Model: ${MODEL_PATH}"
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/"
echo "Threads: ${THREADS}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

echo "[1/1] Running Clair3 (this takes 2-4 hours for 30X long-read WGS)..."
docker run --rm \
  --cpus "${THREADS}" --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  "$CLAIR3_IMAGE" \
  /opt/bin/run_clair3.sh \
    --bam_fn="/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --ref_fn="/genome/reference/Homo_sapiens_assembly38.fasta" \
    --platform="${PLATFORM}" \
    --model_path="${MODEL_PATH}" \
    --output="/genome/${SAMPLE}/vcf_clair3" \
    --threads="${THREADS}" \
    --sample_name="${SAMPLE}"

# Clair3 outputs merge_output.vcf.gz as the final merged VCF
# Rename to match pipeline conventions
CLAIR3_VCF="${OUTPUT_DIR}/merge_output.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

if [ -f "$CLAIR3_VCF" ] && [ "$CLAIR3_VCF" != "$FINAL_VCF" ]; then
  echo "Renaming output to match pipeline conventions..."
  cp "$CLAIR3_VCF" "$FINAL_VCF"
  cp "${CLAIR3_VCF}.tbi" "${FINAL_VCF}.tbi" 2>/dev/null || true
fi

echo "=== Clair3 complete ==="
echo "VCF: ${FINAL_VCF}"
echo ""
echo "Quick stats:"
VARIANT_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/vcf_clair3/${SAMPLE}.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
VARIANT_COUNT=${VARIANT_COUNT:-unknown}
echo "  Total variants: ${VARIANT_COUNT}"
PASS_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view -f PASS "/genome/${SAMPLE}/vcf_clair3/${SAMPLE}.vcf.gz" \
  | grep -vc '^#' || echo "unknown")
echo "  PASS variants: ${PASS_COUNT}"
echo ""
echo "Next steps:"
echo "  - Run downstream VCF steps (ClinVar, PharmCAT, VEP) pointing at vcf_clair3/"
echo "  - Compare with DeepVariant: ./scripts/benchmark-variants.sh ${SAMPLE}"
