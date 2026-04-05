#!/usr/bin/env bash
# [EXPERIMENTAL] Somatic variant calling — Mutect2 tumor-only mode
# Finds mutations acquired during life (not inherited), such as clonal hematopoiesis
# Input: Sorted BAM + GRCh38 reference
# Output: Filtered somatic VCF in $GENOME_DIR/<sample>/somatic/
# Runtime: ~2-6 hours (full genome), ~15-30 min per chromosome with INTERVALS
#
# WARNING: Tumor-only mode (no matched normal) has a HIGH false positive rate.
# Many germline variants will be called as somatic. Use gnomAD and PoN resources
# to reduce false positives, and treat results as exploratory only.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
INTERVALS=${INTERVALS:-""}

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
REF_DICT="${GENOME_DIR}/reference/Homo_sapiens_assembly38.dict"
OUTPUT_DIR="${SAMPLE_DIR}/somatic"

GATK_IMAGE="broadinstitute/gatk:4.6.1.0"
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# Optional resources (improve filtering if present)
GNOMAD_VCF="${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz"
PON_VCF="${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz"

echo "=== [EXPERIMENTAL] Somatic Variant Calling (Mutect2 Tumor-Only): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"
if [ -n "$INTERVALS" ]; then
  echo "Intervals: ${INTERVALS}"
fi
echo "Output: ${OUTPUT_DIR}/"
echo ""

# Check for idempotent skip
FINAL_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz"
if [ -f "$FINAL_OUTPUT" ]; then
  echo "Output already exists: ${FINAL_OUTPUT}"
  echo "Skipping. Delete the file to re-run."
  exit 0
fi

# Validate required inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# GATK needs .dict file
if [ ! -f "$REF_DICT" ]; then
  echo "ERROR: Sequence dictionary not found: ${REF_DICT}" >&2
  echo "Generate it with: docker run --rm --user root -v \${GENOME_DIR}:/genome ${GATK_IMAGE} gatk CreateSequenceDictionary -R /genome/reference/Homo_sapiens_assembly38.fasta" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Build Mutect2 command
MUTECT2_CMD=(
  gatk Mutect2
  -R /genome/reference/Homo_sapiens_assembly38.fasta
  -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  -O "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_unfiltered.vcf.gz"
  --native-pair-hmm-threads "$THREADS"
  --max-mnp-distance 0
)

# Add gnomAD germline resource if available (reduces germline false positives)
if [ -f "$GNOMAD_VCF" ] && [ -f "${GNOMAD_VCF}.tbi" ]; then
  echo "Using gnomAD germline resource: ${GNOMAD_VCF}"
  MUTECT2_CMD+=(--germline-resource "/genome/somatic/af-only-gnomad.hg38.vcf.gz")
else
  echo "WARNING: gnomAD AF-only VCF not found at ${GNOMAD_VCF}"
  echo "  Without gnomAD, many common germline variants will appear as somatic calls."
  echo "  See docs/29-mutect2-somatic.md for download instructions."
  echo ""
fi

# Add Panel of Normals if available (reduces recurrent technical artifacts)
if [ -f "$PON_VCF" ] && [ -f "${PON_VCF}.tbi" ]; then
  echo "Using Panel of Normals: ${PON_VCF}"
  MUTECT2_CMD+=(-pon "/genome/somatic/1000g_pon.hg38.vcf.gz")
else
  echo "INFO: Panel of Normals not found at ${PON_VCF} (optional, reduces artifacts)."
  echo ""
fi

# Add interval restriction if specified
if [ -n "$INTERVALS" ]; then
  MUTECT2_CMD+=(--intervals "$INTERVALS")
fi

echo "=== [1/3] Running Mutect2 in tumor-only mode ==="
echo "  This is the slowest step. Full genome takes ~2-6 hours."
echo "  For quick testing, set INTERVALS=chr22"
echo ""
docker run --rm --user root \
  --cpus "$THREADS" --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  "$GATK_IMAGE" \
  "${MUTECT2_CMD[@]}"

echo ""
echo "=== [2/3] Filtering somatic calls (FilterMutectCalls) ==="
FILTER_CMD=(
  gatk FilterMutectCalls
  -R /genome/reference/Homo_sapiens_assembly38.fasta
  -V "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_unfiltered.vcf.gz"
  -O "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz"
)

docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "$GATK_IMAGE" \
  "${FILTER_CMD[@]}"

echo ""
echo "=== [3/3] Somatic variant statistics ==="

# Count PASS variants
PASS_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view -f PASS "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz" \
  2>/dev/null | grep -c "^[^#]" || echo "0")

TOTAL_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "$BCFTOOLS_IMAGE" \
  bcftools view "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz" \
  2>/dev/null | grep -c "^[^#]" || echo "0")

echo "  Total calls: ${TOTAL_COUNT}"
echo "  PASS calls:  ${PASS_COUNT}"
echo ""

echo "=== [EXPERIMENTAL] Somatic variant calling complete ==="
echo ""
echo "Output files:"
echo "  Unfiltered: ${OUTPUT_DIR}/${SAMPLE}_somatic_unfiltered.vcf.gz"
echo "  Filtered:   ${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz"
echo "  Stats:      ${OUTPUT_DIR}/${SAMPLE}_somatic_unfiltered.vcf.gz.stats"
echo ""
echo "IMPORTANT: Tumor-only mode produces many false positives."
echo "  - Most PASS variants in a healthy individual are germline, not somatic."
echo "  - True somatic variants (e.g., clonal hematopoiesis) typically have low AF (<0.1)."
echo "  - Cross-reference with ClinVar and gnomAD before interpreting any variant."
echo ""
echo "View low-AF PASS variants (potential somatic):"
echo "  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' -i 'FILTER=\"PASS\"' \\"
echo "    ${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz | awk '\$5 < 0.1'"
