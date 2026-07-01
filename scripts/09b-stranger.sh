#!/usr/bin/env bash
# Stranger — Annotate ExpansionHunter STR VCF with clinical pathogenicity status
# Adds STR_STATUS (normal/pre_mutation/full_mutation), disease name, OMIM number,
# inheritance mode, and normal/pathogenic repeat ranges to each locus.
# Input:  ExpansionHunter VCF produced by step 09
# Output: Annotated VCF in $GENOME_DIR/<sample>/expansion_hunter/
#
# This step exits cleanly (exit 0) when the ExpansionHunter VCF does not exist
# so the pipeline can treat it as optional without failing run-all.sh.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"
EH_VCF="${OUTPUT_DIR}/${SAMPLE}_eh.vcf"
OUT_VCF="${OUTPUT_DIR}/${SAMPLE}_eh_stranger.vcf"

echo "=== Stranger STR Annotation: ${SAMPLE} ==="

# Exit cleanly if the ExpansionHunter VCF does not exist.
# Step 09 may have been skipped or run separately — this is not an error.
if [ ! -f "$EH_VCF" ]; then
  echo "INFO: ExpansionHunter VCF not found: ${EH_VCF}"
  echo "INFO: Run scripts/09-expansion-hunter.sh first, or skip this annotation step."
  exit 0
fi

# Skip if output already exists
if [ -f "$OUT_VCF" ]; then
  echo "Stranger output already exists: ${OUT_VCF}"
  echo "Delete to re-run: rm ${OUT_VCF}"
  exit 0
fi

echo "Input:  ${EH_VCF}"

# Stranger annotates each STR locus with clinical pathogenicity thresholds from its
# bundled repeat catalog (derived from ClinGen/OMIM data).
# An optional custom catalog (TSV) can be supplied via STRANGER_REPEATS.
if [ -n "${STRANGER_REPEATS:-}" ]; then
  if [ ! -f "${STRANGER_REPEATS}" ]; then
    echo "ERROR: STRANGER_REPEATS catalog not found: ${STRANGER_REPEATS}" >&2
    exit 1
  fi
  # Compute container-relative path for the catalog file
  STRANGER_REPEATS_REL="${STRANGER_REPEATS#"${GENOME_DIR}/"}"
  echo "Repeat catalog: ${STRANGER_REPEATS} (custom)"
  docker run --rm \
    --cpus 1 --memory 1g \
    -v "${GENOME_DIR}:/genome" \
    "${STRANGER_IMAGE}" \
    stranger \
      --repeats-file "/genome/${STRANGER_REPEATS_REL}" \
      "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf" \
    > "${OUT_VCF}"
else
  echo "Repeat catalog: bundled clinical database (default)"
  docker run --rm \
    --cpus 1 --memory 1g \
    -v "${GENOME_DIR}:/genome" \
    "${STRANGER_IMAGE}" \
    stranger \
      "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf" \
    > "${OUT_VCF}"
fi

echo "=== Stranger complete ==="
echo "Annotated VCF: ${OUT_VCF}"
echo ""
echo "Each locus now has:"
echo "  STR_STATUS: normal / pre_mutation / full_mutation"
echo "  Disease name, OMIM number, inheritance mode, and repeat size ranges"
echo "  See docs/09b-stranger.md for interpretation guidance."
