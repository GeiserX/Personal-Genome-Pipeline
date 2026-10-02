#!/usr/bin/env bash
# 21-cyrius.sh — [EXPERIMENTAL] CYP2D6 star allele calling using Cyrius
# Usage: ./scripts/21-cyrius.sh <sample_name>
#
# CYP2D6 is the hardest pharmacogene to call because of its pseudogene (CYP2D7)
# and complex structural variants (gene deletions, duplications, hybrids).
# Cyrius uses depth-based analysis specifically designed for CYP2D6.
#
# EXPERIMENTAL: Cyrius 1.1.1 is installed at runtime with pip (its dependencies
# pinned by scripts/cyrius-constraints.txt, so it needs network access) and
# may return "None" for complex CYP2D6 arrangements. Verify results against
# PharmCAT or clinical lab calls before acting on them.
#
# Requires: Sorted BAM with index
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

CONSTRAINTS="${SCRIPT_DIR}/cyrius-constraints.txt"
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
BAI="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai"
OUTDIR="${GENOME_DIR}/${SAMPLE}/cyrius"
mkdir -p "$OUTDIR"

# Validate inputs
for FILE in "$BAM" "$BAI" "$CONSTRAINTS"; do
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

# Cyrius is a Python tool with no maintained image. We install the pinned
# release in the Python container, its dependencies held to the versions in
# cyrius-constraints.txt, and keep pip's own errors in the log. The package's
# console script is `cyrius` (there is no `star_caller` command).
# The manifest file (list of BAM paths) is created inside the container.
echo "[1/2] Running Cyrius CYP2D6 caller..."
# --net: pip downloads Cyrius and its dependencies from PyPI. pip installs
# them for the calling user under HOME (/tmp in the container), so the step
# needs no root and its outputs belong to the caller.
run_in --net --cpus 4 --memory 8g \
  -v "${CONSTRAINTS}:/constraints.txt:ro" \
  -w /tmp \
  "${PYTHON_IMAGE}" \
  bash -c "
    pip install --user --no-cache-dir --disable-pip-version-check -q -c /constraints.txt 'cyrius==${CYRIUS_VERSION}' &&
    export PATH=\"\$HOME/.local/bin:\$PATH\" &&
    echo '/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam' > /tmp/manifest.txt &&
    cyrius \
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
  echo "ERROR: Cyrius finished but wrote no ${RESULT_FILE}. Check the messages above." >&2
  exit 1
fi

echo ""
echo "============================================"
echo "  CYP2D6 calling complete: ${SAMPLE}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
