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

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
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

# Skip only a finished output: written through a temporary name, so an empty
# or cut-short file from a failed run is never taken for one.
if have_output "$OUT_VCF"; then
  echo "Stranger output already exists: ${OUT_VCF}"
  echo "Delete to re-run: rm ${OUT_VCF}"
  exit 0
fi

echo "Input:  ${EH_VCF}"

# Stranger annotates each STR locus with clinical pathogenicity thresholds from a
# repeat catalog. Its default is the GRCh37 catalog it bundles, while step 09
# calls the GRCh38 loci, so the bundled GRCh38 catalog is passed explicitly
# (found next to the installed package). A custom catalog (TSV or JSON) can be
# supplied via STRANGER_REPEATS.
if [ -n "${STRANGER_REPEATS:-}" ]; then
  if [ ! -f "${STRANGER_REPEATS}" ]; then
    echo "ERROR: STRANGER_REPEATS catalog not found: ${STRANGER_REPEATS}" >&2
    exit 1
  fi
  REPEATS_C=$(cpath "$STRANGER_REPEATS") || exit 2
  echo "Repeat catalog: ${STRANGER_REPEATS} (custom)"
  atomic_out "$OUT_VCF" run_in \
    --cpus 1 --memory 1g \
    "${STRANGER_IMAGE}" \
    stranger \
      --repeats-file "$REPEATS_C" \
      "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
else
  echo "Repeat catalog: bundled GRCh38 clinical catalog (variant_catalog_grch38.json)"
  # shellcheck disable=SC2016  # the dollar signs are for the container's sh
  atomic_out "$OUT_VCF" run_in \
    --cpus 1 --memory 1g \
    "${STRANGER_IMAGE}" \
    sh -c 'catalog=$(python -c "import os, stranger; print(os.path.join(os.path.dirname(stranger.__file__), \"resources\", \"variant_catalog_grch38.json\"))") &&
      [ -f "$catalog" ] || { echo "ERROR: no GRCh38 catalog in the Stranger image (${catalog:-not found})" >&2; exit 1; }
      exec stranger --repeats-file "$catalog" "$1"' \
    stranger "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf"
fi

echo "=== Stranger complete ==="
echo "Annotated VCF: ${OUT_VCF}"
echo ""
echo "Each locus now has:"
echo "  STR_STATUS: normal / pre_mutation / full_mutation"
echo "  Disease name, OMIM number, inheritance mode, and repeat size ranges"
echo "  See docs/09b-stranger.md for interpretation guidance."
