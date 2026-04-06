#!/usr/bin/env bash
# 21-cyrius.sh — [EXPERIMENTAL] CYP2D6 star allele calling using Cyrius
# Usage: ./scripts/21-cyrius.sh <sample_name>
#
# CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
# and complex structural variants (gene deletions, duplications, hybrids).
# Cyrius uses depth-based analysis specifically designed for CYP2D6.
#
# EXPERIMENTAL: Cyrius is installed at runtime via pip (unpinned version) and
# may return "None" for complex CYP2D6 arrangements. Verify results against
# PharmCAT or clinical lab calls before acting on them.
#
# Requires: Sorted BAM with index
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

# shellcheck source=versions.env
source "$(dirname "$0")/versions.env"

BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
BAI="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai"
OUTDIR="${GENOME_DIR}/${SAMPLE}/cyrius"
mkdir -p "$OUTDIR"

# Validate inputs
for FILE in "$BAM" "$BAI"; do
  if [ ! -f "$FILE" ]; then
    echo "ERROR: Required file not found: ${FILE}"
    exit 1
  fi
done

echo "============================================"
echo "  Step 21: CYP2D6 Star Allele Calling"
echo "  Tool: Cyrius (Illumina)"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${BAM}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Cyrius is a Python tool. We use a Python container and install it on the fly.
# This avoids dependency on a specific Cyrius Docker image that may not exist.
# The manifest file (list of BAM paths) is created inside the container.
echo "[1/2] Running Cyrius CYP2D6 caller..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  -w /tmp \
  "${PYTHON_IMAGE}" \
  bash -c "
    pip install -q cyrius 2>/dev/null &&
    echo '/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam' > /tmp/manifest.txt &&
    star_caller \
      --manifest /tmp/manifest.txt \
      --genome 38 \
      --prefix ${SAMPLE}_cyp2d6 \
      --outDir /genome/${SAMPLE}/cyrius/ \
      --threads 4
  "

echo ""
echo "[2/2] Parsing results..."

RESULT_FILE="${OUTDIR}/${SAMPLE}_cyp2d6.tsv"
if [ -f "$RESULT_FILE" ]; then
  echo ""
  echo "  CYP2D6 Results:"
  echo "  ─────────────────"
  # Display results
  column -t "$RESULT_FILE" 2>/dev/null || cat "$RESULT_FILE"
  echo ""

  # Extract diplotype and phenotype
  DIPLOTYPE=$(awk -F'\t' 'NR==2 {print $2}' "$RESULT_FILE" 2>/dev/null || echo "N/A")
  echo "  Diplotype: ${DIPLOTYPE}"
  echo ""
  echo "  Interpret this with PharmGKB: https://www.pharmgkb.org/gene/PA128"
else
  echo "  WARNING: No output file found. Cyrius may have failed."
  echo "  Check the error messages above."
fi

echo ""
echo "============================================"
echo "  CYP2D6 calling complete: ${SAMPLE}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
