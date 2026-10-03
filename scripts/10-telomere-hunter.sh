#!/usr/bin/env bash
# TelomereHunter — Estimate telomere length from WGS BAM
# Output: tel_content metric (GC-corrected telomeric reads per million)
# Higher values = longer telomeres. Provides biological age baseline.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/telomere/${SAMPLE}"

echo "=== TelomereHunter: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"
echo "WARNING: This reads the entire BAM (~30-40GB). Takes 30-60 minutes."

for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# TelomereHunter sorts reads into intratelomeric, subtelomeric and junction
# classes by chromosome band, and without -b it uses its own hg19 bands
# (`telomerehunter --help`: "If no banding file is specified, the banding
# information of hg19 will be used"). The BAM is GRCh38, so UCSC's GRCh38
# bands (installed by setup.sh) are passed.
BAND_ARGS=()
if BANDS=$(data_file cytoband); then
  echo "Chromosome bands: ${BANDS}"
  BAND_ARGS=(-b "$(cpath "$BANDS")")
else
  echo "WARNING: the GRCh38 chromosome bands are not installed (${BANDS})."
  echo "  TelomereHunter falls back to its hg19 bands, so the subtelomeric and"
  echo "  junction read classes use hg19 band ends on GRCh38 positions."
  echo "  Install the bands with: ./scripts/setup.sh ${GENOME_DIR}"
fi

# --root: this image has not been shown to run as an unprivileged user.
# TelomereHunter has no thread option; THREADS caps the container.
run_in --root \
  --cpus "$THREADS" --memory 4g \
  "${TELOMEREHUNTER_IMAGE}" \
  telomerehunter \
    -ibt "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -o "/genome/${SAMPLE}/telomere/${SAMPLE}" \
    -p "$SAMPLE" \
    ${BAND_ARGS[@]+"${BAND_ARGS[@]}"}

echo "=== TelomereHunter complete ==="
SUMMARY="${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_summary.tsv"
if [ -f "$SUMMARY" ]; then
  echo "Summary: $SUMMARY"
  TEL_CONTENT=$(awk -F'\t' 'NR==2 {print $11}' "$SUMMARY")
  echo "Telomere content: ${TEL_CONTENT}"
fi
