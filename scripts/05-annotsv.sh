#!/usr/bin/env bash
# AnnotSV — Annotate structural variants with ACMG pathogenicity classification
# Input: Manta diploidSV.vcf.gz
# Output: *_sv_annotated.tsv (ACMG class 1-5 for each SV)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
MANTA_VCF="${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz"
# Fall back to manta2/ if a second Manta run was used
[ ! -f "$MANTA_VCF" ] && MANTA_VCF="${SAMPLE_DIR}/manta2/results/variants/diploidSV.vcf.gz"
OUTPUT_DIR="${SAMPLE_DIR}/annotsv"

echo "=== AnnotSV: ${SAMPLE} ==="
echo "Input: ${MANTA_VCF}"
echo "Output: ${OUTPUT_DIR}/"

if [ ! -f "$MANTA_VCF" ]; then
  echo "ERROR: Manta VCF not found: ${MANTA_VCF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Determine relative path of Manta VCF within SAMPLE_DIR
MANTA_REL=$(echo "$MANTA_VCF" | sed "s|${GENOME_DIR}/||")

docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  getwilds/annotsv:latest \
  AnnotSV \
    -SVinputFile "/genome/${MANTA_REL}" \
    -outputFile "/genome/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv" \
    -genomeBuild GRCh38 \
    -annotationMode both

echo "=== AnnotSV complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv_annotated.tsv"
echo ""
echo "Filter pathogenic (class 4-5, <5MB):"
echo "  awk -F'\t' 'NR==1 || (\$120==4 || \$120==5) && \$16==\"full\"' ${OUTPUT_DIR}/${SAMPLE}_sv_annotated.tsv"
