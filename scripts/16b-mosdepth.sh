#!/usr/bin/env bash
# mosdepth — fast per-base and per-region coverage statistics from BAM
# Input: sorted BAM + BAI from $GENOME_DIR/<sample>/aligned/
# Output: coverage distributions, thresholds, summary in $GENOME_DIR/<sample>/mosdepth/
#
# Complements indexcov (step 16): indexcov reads only the BAM index for a quick
# estimate; mosdepth reads actual alignments for precise per-base depth.
#
# A second pass measures how well each gene of assets/callability_genes.tsv
# (CALLABILITY_GENES names another table) could be read: for each gene on
# chr1-22, X and Y, the bases with 10 or more reads of mapping quality 20 or
# more, in <sample>_callability.bed.gz. The report lists the genes where
# fewer than 90% of the bases reach that (docs/16b-mosdepth.md).
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
OUTPUT_DIR="${SAMPLE_DIR}/mosdepth"
CALLABILITY_GENES=${CALLABILITY_GENES:-${PGP_ROOT}/assets/callability_genes.tsv}
CALLABLE="${OUTPUT_DIR}/${SAMPLE}_callability.bed.gz"

echo "=== mosdepth Coverage: ${SAMPLE} ==="
echo "BAM: ${BAM}"

# Validate inputs
if [ ! -f "$BAM" ]; then
  echo "ERROR: BAM not found: ${BAM}" >&2
  exit 1
fi
if [ ! -f "${BAM}.bai" ]; then
  echo "ERROR: BAM index not found: ${BAM}.bai" >&2
  echo "Create it with: samtools index ${BAM}" >&2
  exit 1
fi

if [ ! -s "$CALLABILITY_GENES" ]; then
  echo "ERROR: gene table not found or empty: ${CALLABILITY_GENES}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# callability_pass: mosdepth over the genes of the table, reads of MAPQ 20 or
# more only, the bases at 10x or more per gene (the same command as the
# MOSDEPTH module). The genes with no position on chr1-22, X or Y (ALT-only)
# are left out of the BED; the report reads them from the table.
callability_pass() {
  local work="${OUTPUT_DIR}/callability.part"
  rm -rf "$work"
  mkdir -p "$work"
  awk -F'\t' -v OFS='\t' '!/^#/ && $1 != "gene" && $4 ~ /^[0-9]+$/ {print $3, $4, $5, $1}' \
    "$CALLABILITY_GENES" > "${work}/genes.bed"
  echo "Gene callability (MAPQ >= 20, 10x) over $(grep -c . "${work}/genes.bed") genes of ${CALLABILITY_GENES}..."
  run_in --cpus "${THREADS}" --memory 4g \
    "${MOSDEPTH_IMAGE}" \
    mosdepth \
      --by "/genome/${SAMPLE}/mosdepth/callability.part/genes.bed" \
      --mapq 20 \
      --fast-mode \
      --no-per-base \
      --threads "${THREADS}" \
      --thresholds 10 \
      "/genome/${SAMPLE}/mosdepth/callability.part/c" \
      "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  if [ -f "${work}/c.thresholds.bed.gz" ]; then
    mv "${work}/c.thresholds.bed.gz" "$CALLABLE"
  else
    echo "WARNING: mosdepth wrote no gene table; the report says the gene callability was not measured." >&2
  fi
  rm -rf "$work"
}

# Skip if output already exists. mosdepth writes the summary last, so a
# non-empty summary means the run finished. The gene pass runs on its own
# when only its table is missing (a folder from before it existed).
if have_output "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"; then
  echo "mosdepth output already exists in ${OUTPUT_DIR}/, skipping."
  if ! have_output "$CALLABLE"; then
    callability_pass
  fi
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

# mosdepth flags:
#   --by 500        Window size for per-region coverage (500bp bins; WGS default)
#   --by <BED>      Capture BED for on-target coverage (WES; set CAPTURE_BED env var)
#   --fast-mode     Use faster algorithm (does NOT skip per-base output by itself)
#   --no-per-base   Omit the large per-base BED (saves ~3-5 GB per WGS sample)
#   --threads       Decompression threads (mosdepth uses 1 main + N decompression)
#   --thresholds    Report fraction of bases at these coverage thresholds

# Use capture BED for WES on-target coverage, or 500bp bins for WGS
if [ -n "${CAPTURE_BED:-}" ]; then
  if [ ! -f "${CAPTURE_BED}" ]; then
    echo "ERROR: CAPTURE_BED not found: ${CAPTURE_BED}" >&2
    exit 1
  fi
  # Compute container-relative path for the BED file
  CAPTURE_BED_REL="${CAPTURE_BED#"${GENOME_DIR}/"}"
  BY_FLAG="/genome/${CAPTURE_BED_REL}"
  echo "Using capture BED for WES on-target coverage: ${CAPTURE_BED}"
else
  BY_FLAG="500"
fi

echo "Computing coverage statistics..."
run_in --cpus "${THREADS}" --memory 4g \
  "${MOSDEPTH_IMAGE}" \
  mosdepth \
    --by "${BY_FLAG}" \
    --fast-mode \
    --no-per-base \
    --threads "${THREADS}" \
    --thresholds 1,5,10,15,20,30,50 \
    "/genome/${SAMPLE}/mosdepth/${SAMPLE}" \
    "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"

rm -f "$CALLABLE"
callability_pass

echo "=== mosdepth complete ==="
echo "Summary:        ${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
echo "Distribution:   ${OUTPUT_DIR}/${SAMPLE}.mosdepth.global.dist.txt"
echo "Region coverage: ${OUTPUT_DIR}/${SAMPLE}.regions.bed.gz"
echo "Thresholds:     ${OUTPUT_DIR}/${SAMPLE}.thresholds.bed.gz"
echo "Gene callability: ${CALLABLE}"

# Print summary
if [ -f "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt" ]; then
  echo ""
  echo "Coverage summary:"
  head -5 "${OUTPUT_DIR}/${SAMPLE}.mosdepth.summary.txt"
fi
