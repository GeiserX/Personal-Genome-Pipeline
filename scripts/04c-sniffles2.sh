#!/usr/bin/env bash
# Sniffles2 — Long-read structural variant caller
# Designed for ONT and PacBio long-read data. Detects SVs (>50bp): deletions,
# duplications, insertions, inversions, translocations, and nested SVs.
# Alternative to step 04 (Manta, short-read only). Outputs to sv_sniffles/ to avoid conflicts.
# Input: sorted BAM from long-read alignment + GRCh38 reference
# Output: SV VCF in $GENOME_DIR/<sample>/sv_sniffles/
# Runtime: ~30-90 minutes per 30X long-read genome
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned_longread}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/sv_sniffles"

SNIFFLES_IMAGE="quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

echo "=== Sniffles2 SV Calling: ${SAMPLE} ==="
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

# Run Sniffles2
echo "[1/3] Running Sniffles2 SV caller..."
echo "       This takes 30-90 minutes for 30X long-read WGS."
docker run --rm --user root \
  --cpus "$THREADS" --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  "$SNIFFLES_IMAGE" \
  sniffles \
    -i "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -v "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv_raw.vcf" \
    --reference /genome/reference/Homo_sapiens_assembly38.fasta \
    --threads "$THREADS" \
    --sample-id "${SAMPLE}"

# Compress and index with bcftools
echo "[2/3] Compressing VCF..."
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view \
    "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv_raw.vcf" \
    -Oz -o "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz"

echo "[3/3] Indexing VCF..."
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools index -t \
    "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz"

# Clean up raw VCF
rm -f "${OUTPUT_DIR}/${SAMPLE}_sv_raw.vcf"

SV_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
SV_COUNT=${SV_COUNT:-unknown}

echo "=== Sniffles2 complete ==="
echo "Total SVs called: ${SV_COUNT}"
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz"
echo ""
echo "Count by SV type:"
docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools query -f '%INFO/SVTYPE\n' \
    "/genome/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz" \
  | sort | uniq -c | sort -rn || echo "  (run bcftools query manually)"
echo ""
echo "Next steps:"
echo "  - Annotate SVs: ./scripts/05-annotsv.sh ${SAMPLE} (set SV_VCF to point at sv_sniffles)"
echo "  - Merge with short-read SVs: ./scripts/22-survivor-merge.sh ${SAMPLE}"
echo "  - Filter quality: duphold can re-annotate long-read SVs with GQ scores"
