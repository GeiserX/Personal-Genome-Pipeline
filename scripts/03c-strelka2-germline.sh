#!/usr/bin/env bash
# Strelka2 — Alternative germline small variant caller (SNVs + indels)
# Alternative to step 03 (DeepVariant). Outputs to vcf_strelka2/ to avoid conflicts.
# Strelka2 is a SMALL VARIANT caller (SNVs + indels up to ~49 bp), NOT a structural
# variant caller. It complements Manta, which handles SVs (>50 bp).
# Input: sorted BAM + GRCh38 reference
# Output: germline variants VCF in $GENOME_DIR/<sample>/vcf_strelka2/
# Runtime: ~1-2 hours per 30X genome
#
# NOTE: Strelka2's scoring model is optimized for BWA-MEM alignments.
# When used with minimap2-aligned BAMs, SNP precision may be reduced due to:
#   - Missing XS (suboptimal alignment score) tags that BWA-MEM produces
#   - Different AS (alignment score) scaling between minimap2 and BWA-MEM
# For best results, use BWA-MEM2 alignments (scripts/02a-alignment-bwamem2.sh).
#
# Reference: Kim et al. Strelka2: fast and accurate calling of germline and
# somatic variants. Nature Methods (2018). https://doi.org/10.1038/s41592-018-0051-x
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
OUTPUT_DIR="${SAMPLE_DIR}/vcf_strelka2"
THREADS=${THREADS:-8}

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

STRELKA_IMAGE="${STRELKA_IMAGE}"
BCFTOOLS_IMAGE="${BCFTOOLS_IMAGE}"

# A finished run is not repeated: the configure script refuses a runDir
# that already holds a workflow, so a second run used to fail.
if [ -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz" ] && [ -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz.tbi" ]; then
  echo "Strelka2 already done for ${SAMPLE}: ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
  echo "To run it again, delete ${OUTPUT_DIR}/ first."
  exit 0
fi

# Step 1: Configure Strelka2 germline workflow (skipped when an interrupted
# run can be resumed)
if [ -f "${OUTPUT_DIR}/runWorkflow.py" ]; then
  echo "[1/2] Found an unfinished Strelka2 run in ${OUTPUT_DIR}/; resuming it."
else
  if [ -d "$OUTPUT_DIR" ]; then
    # Files written by the container belong to root, so remove them from a container
    echo "Removing leftover ${OUTPUT_DIR}/ (no workflow and no results)..."
    run_in      "$STRELKA_IMAGE" \
      rm -rf "/genome/${SAMPLE}/vcf_strelka2"
  fi
  mkdir -p "$OUTPUT_DIR"
  echo "[1/2] Configuring Strelka2 germline workflow..."
  run_in    --cpus "$THREADS" --memory 16g \
    "$STRELKA_IMAGE" \
    configureStrelkaGermlineWorkflow.py \
      --bam "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
      --referenceFasta "${REF_FASTA_C}" \
      --runDir "/genome/${SAMPLE}/vcf_strelka2"
fi

# Step 2: Run the workflow
echo "[2/2] Running Strelka2 (this takes 1-2 hours for 30X WGS)..."
run_in  --cpus "$THREADS" --memory 16g \
  "$STRELKA_IMAGE" \
  "/genome/${SAMPLE}/vcf_strelka2/runWorkflow.py" \
    -m local \
    -j "$THREADS"

if [ ! -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz" ] || [ ! -f "${OUTPUT_DIR}/results/variants/variants.vcf.gz.tbi" ]; then
  echo "ERROR: Strelka2 finished without ${OUTPUT_DIR}/results/variants/variants.vcf.gz(.tbi)." >&2
  echo "  Delete ${OUTPUT_DIR}/ and run this step again." >&2
  exit 1
fi

VARIANT_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
VARIANT_COUNT=${VARIANT_COUNT:-unknown}

echo "=== Strelka2 complete ==="
echo "Total variants called: ${VARIANT_COUNT}"
echo "Results: ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "View PASS variants only:"
echo "  bcftools view -f PASS ${OUTPUT_DIR}/results/variants/variants.vcf.gz"
echo ""
echo "Count variant types:"
echo "  bcftools stats ${OUTPUT_DIR}/results/variants/variants.vcf.gz | grep 'number of records'"
