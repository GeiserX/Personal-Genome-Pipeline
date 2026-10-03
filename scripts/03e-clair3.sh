#!/usr/bin/env bash
# Clair3 — Long-read variant calling (BAM to VCF)
# Alternative to step 03 (DeepVariant) for long-read data.
# Supports Oxford Nanopore (ONT) and PacBio HiFi platforms.
# Input: sorted BAM from long-read alignment + GRCh38 reference
# Output: VCF.gz with SNPs and small indels in $GENOME_DIR/<sample>/vcf_clair3/
#
# Usage: PLATFORM=ont|hifi ./scripts/03e-clair3.sh <sample_name> [male|female]
#
# Set PLATFORM to select the appropriate model:
#   PLATFORM=ont   -> ONT R10.4.1 model (r1041_e82_400bps_sup_v500)
#   PLATFORM=hifi  -> PacBio HiFi/Revio model (hifi_revio)
# CLAIR3_MODEL overrides the model path inside the image (e.g.
#   /opt/models/r1041_e82_400bps_sup_v520 for newer Dorado basecalls).
# Sex: `male` calls chrX and chrY haploid outside the pseudoautosomal regions
#   (Clair3 --gender male with assets/par_grch38.bed); `female` leaves chrY
#   out and calls chrX diploid; no sex calls both diploid.
#
# Runtime: ~2-4 hours for 30X long-read WGS
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
SEX=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
case "$SEX" in
  ''|male|female) ;;
  *) echo "ERROR: sex must be 'male' or 'female' (or left out), got '${SEX}'." >&2; exit 2 ;;
esac
PAR_BED="${PGP_ROOT}/assets/par_grch38.bed"
PLATFORM=${PLATFORM:?Set PLATFORM to ont or hifi}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned_longread}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/vcf_clair3"


# Select model path based on platform
case "$PLATFORM" in
  ont)
    MODEL_PATH="${CLAIR3_MODEL:-/opt/models/r1041_e82_400bps_sup_v500}"
    ;;
  hifi)
    MODEL_PATH="${CLAIR3_MODEL:-/opt/models/hifi_revio}"
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
echo "Sex: ${SEX:-not given (chrX and chrY called diploid)}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# --gender exists from Clair3 v2.0.3 on; male also gets the PAR BED, so the
# pseudoautosomal regions stay diploid.
SEX_ARGS=()
case "$SEX" in
  male) SEX_ARGS=(--gender=male --par_regions_bed=/pgp/par_grch38.bed) ;;
  female) SEX_ARGS=(--gender=female) ;;
esac

echo "[1/1] Running Clair3 (this takes 2-4 hours for 30X long-read WGS)..."
run_in \
  --cpus "${THREADS}" --memory 32g \
  -v "${PAR_BED}:/pgp/par_grch38.bed:ro" \
  "$CLAIR3_IMAGE" \
  /opt/bin/run_clair3.sh \
    --bam_fn="/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --ref_fn="${REF_FASTA_C}" \
    --platform="${PLATFORM}" \
    --model_path="${MODEL_PATH}" \
    --output="/genome/${SAMPLE}/vcf_clair3" \
    --threads="${THREADS}" \
    --sample_name="${SAMPLE}" \
    ${SEX_ARGS[@]+"${SEX_ARGS[@]}"}

# Clair3 writes merge_output.vcf.gz (and its .tbi) as the final merged VCF;
# it is moved to the pipeline's name. Without it the run failed, whatever
# Clair3's exit code said, and an older <sample>.vcf.gz is not a result.
CLAIR3_VCF="${OUTPUT_DIR}/merge_output.vcf.gz"
FINAL_VCF="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"

if ! wrote_vcf "$CLAIR3_VCF"; then
  echo "ERROR: Clair3 left no ${CLAIR3_VCF}; see its log in ${OUTPUT_DIR}/." >&2
  exit 1
fi
echo "Renaming output to match pipeline conventions..."
rm -f "${FINAL_VCF}.tbi"
mv -f "$CLAIR3_VCF" "$FINAL_VCF"
if [ -f "${CLAIR3_VCF}.tbi" ]; then
  mv -f "${CLAIR3_VCF}.tbi" "${FINAL_VCF}.tbi"
else
  run_in "$BCFTOOLS_IMAGE" bcftools index -f -t "/genome/${SAMPLE}/vcf_clair3/${SAMPLE}.vcf.gz"
fi

echo "=== Clair3 complete ==="
echo "VCF: ${FINAL_VCF}"
echo ""
echo "Quick stats:"
VARIANT_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/vcf_clair3/${SAMPLE}.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
VARIANT_COUNT=${VARIANT_COUNT:-unknown}
echo "  Total variants: ${VARIANT_COUNT}"
PASS_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools view -f PASS "/genome/${SAMPLE}/vcf_clair3/${SAMPLE}.vcf.gz" \
  | grep -vc '^#' || echo "unknown")
echo "  PASS variants: ${PASS_COUNT}"
echo ""
echo "Next steps:"
echo "  - Run downstream VCF steps (ClinVar, PharmCAT, VEP) pointing at vcf_clair3/"
echo "  - Compare with DeepVariant: ./scripts/benchmark-variants.sh ${SAMPLE}"
