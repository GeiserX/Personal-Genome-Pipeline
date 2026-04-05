#!/usr/bin/env bash
# TIDDIT — Alternative SV caller (large structural variants)
# Alternative to step 04 (Manta). Outputs to sv_tiddit/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference
# Output: SV VCF in $GENOME_DIR/<sample>/sv_tiddit/
# Runtime: ~30-60 minutes per 30X genome (with --skip_assembly)
# NOTE: TIDDIT >=3.9 requires BWA index for local assembly. Since the default
# pipeline uses minimap2, we use --skip_assembly. If using BWA-MEM2 alignment
# (02a-alignment-bwamem2.sh), you can remove --skip_assembly for better breakpoint resolution.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/sv_tiddit"

echo "=== TIDDIT SV Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

TIDDIT_IMAGE="quay.io/biocontainers/tiddit:3.9.5--py312h6e8b409_0"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# Detect BWA index — if present, enable local assembly for better breakpoint resolution
BWA_INDEX="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta.bwt.2bit.64"
TIDDIT_EXTRA_ARGS=()
if [ -f "$BWA_INDEX" ]; then
  echo "BWA index detected — enabling local assembly for breakpoint refinement."
else
  echo "No BWA index found — using --skip_assembly (minimap2 alignment)."
  TIDDIT_EXTRA_ARGS+=(--skip_assembly)
fi

echo "[1/3] Running TIDDIT SV caller..."
docker run --rm --user root \
  --cpus "$THREADS" --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$TIDDIT_IMAGE" \
  tiddit --sv \
    --bam "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --ref /genome/reference/Homo_sapiens_assembly38.fasta \
    --threads "$THREADS" \
    "${TIDDIT_EXTRA_ARGS[@]}" \
    -o "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}"

echo "[2/3] Compressing VCF with bcftools..."
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view \
    "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}.vcf" \
    -Oz -o "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz"

echo "[3/3] Indexing VCF..."
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools index -t \
    "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz"

SV_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
SV_COUNT=${SV_COUNT:-unknown}

echo "=== TIDDIT complete ==="
echo "Total SVs called: ${SV_COUNT}"
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz"
echo "Auxiliary files: ${OUTPUT_DIR}/${SAMPLE}.ploidies.tab, ${OUTPUT_DIR}/${SAMPLE}.signals.tab"
echo ""
echo "Count by SV type:"
printf '  bcftools query -f "%%INFO/SVTYPE\\n" %s/%s_sv.vcf.gz | sort | uniq -c\n' "${OUTPUT_DIR}" "${SAMPLE}"
