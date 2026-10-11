#!/usr/bin/env bash
# AnnotSV — Annotate structural variants with ACMG pathogenicity classification
# Input: step 15's duphold-filtered calls (duphold/<sample>_sv_filtered.vcf.gz),
#   the input of the Nextflow ANNOTSV; without them, Manta's diploidSV.vcf.gz,
#   with a note. SV_VCF overrides both.
# Output: *_sv_annotated.tsv (ACMG class 1-5 for each SV)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
# Allow SV_VCF env var override (e.g., for Sniffles2 long-read SVs)
if [ -n "${SV_VCF:-}" ]; then
  MANTA_VCF="$SV_VCF"
else
  MANTA_VCF="${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz"
  # Fall back to manta2/ if a second Manta run was used
  [ ! -f "$MANTA_VCF" ] && MANTA_VCF="${SAMPLE_DIR}/manta2/results/variants/diploidSV.vcf.gz"
  # Step 15's depth-filtered calls, unless Manta's calls are newer than them
  FILTERED="${SAMPLE_DIR}/duphold/${SAMPLE}_sv_filtered.vcf.gz"
  if [ -f "$FILTERED" ] && [ -f "${FILTERED}.tbi" ] && ! [ "$MANTA_VCF" -nt "$FILTERED" ]; then
    MANTA_VCF="$FILTERED"
  elif [ -f "$FILTERED" ]; then
    echo "NOTE: ${FILTERED} is unfinished or older than Manta's calls; annotating Manta's calls."
    echo "  Run scripts/15-duphold.sh ${SAMPLE} first for the depth-filtered list the pipeline annotates."
  else
    echo "NOTE: no duphold-filtered calls (scripts/15-duphold.sh); annotating Manta's calls without the depth filter."
  fi
fi
OUTPUT_DIR="${SAMPLE_DIR}/annotsv"
# The AnnotSV image holds code only; its annotation data is a separate
# download that setup.sh puts here.
ANNOTATIONS_DIR="${GENOME_DIR}/annotsv_annotations"

echo "=== AnnotSV: ${SAMPLE} ==="
echo "Input: ${MANTA_VCF}"
echo "Annotations: ${ANNOTATIONS_DIR}/"
echo "Output: ${OUTPUT_DIR}/"

if [ ! -d "${ANNOTATIONS_DIR}/Annotations_Human/Genes/GRCh38" ]; then
  echo "ERROR: AnnotSV annotation data not found in ${ANNOTATIONS_DIR}/Annotations_Human/" >&2
  echo "  Run ./scripts/setup.sh ${GENOME_DIR} to download it (~5.3 GB), or see docs/05-annotsv.md." >&2
  exit 1
fi

if [ ! -f "$MANTA_VCF" ]; then
  echo "ERROR: SV VCF not found: ${MANTA_VCF}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Determine relative path of Manta VCF within SAMPLE_DIR
MANTA_REL=$(echo "$MANTA_VCF" | sed "s|${GENOME_DIR}/||")

# AnnotSV builds sorted copies of its annotation files inside the annotations
# directory the first time it runs, so that directory is writable here.
run_in --rw "$ANNOTATIONS_DIR" --cpus 4 --memory 8g \
  "${ANNOTSV_IMAGE}" \
  AnnotSV \
    -SVinputFile "/genome/${MANTA_REL}" \
    -outputFile "/genome/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv" \
    -genomeBuild GRCh38 \
    -annotationMode both \
    -annotationsDir /genome/annotsv_annotations

echo "=== AnnotSV complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv_annotated.tsv"
echo ""
echo "Filter pathogenic (class 4-5, <5MB):"
echo "  awk -F'\t' 'NR==1 || (\$120==4 || \$120==5) && \$16==\"full\"' ${OUTPUT_DIR}/${SAMPLE}_sv_annotated.tsv"
