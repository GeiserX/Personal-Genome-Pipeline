#!/usr/bin/env bash
# 37-y-haplogroup.sh — Y-chromosome haplogroup (paternal line) with Yleaf
# Usage: ./scripts/37-y-haplogroup.sh <sample_name>
#
# Opt-in. Yleaf (YLEAF_IMAGE, 3.2.1: the only biocontainer; upstream is 4.x)
# reads the BAM's pileup at its Y-chromosome markers and predicts the
# haplogroup they support. The image lacks the samtools Yleaf calls and its
# marker tables (setup.sh --yleaf-data installs them), so the pileup is made
# in SAMTOOLS_IMAGE and bin/yleaf_run.py runs Yleaf on it. Runs only when step 16 (indexcov) infers male from
# the reads: the sex comes from indexcov/indexcov-indexcov.ped, not from an
# argument. A female sample, or one with no Y to read, is skipped with one line.
#
# Output: y_haplogroup/<sample>_y_haplogroup.txt, Yleaf's prediction table
# (Hg, Hg_marker, Total_reads, Valid_markers, QC-score). Hg NA means too few
# markers had reads to call one (the step says so).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image YLEAF_IMAGE SAMTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
PED="${S}/indexcov/indexcov-indexcov.ped"
OUTDIR="${S}/y_haplogroup"
OUT="${OUTDIR}/${SAMPLE}_y_haplogroup.txt"
YDATA="${GENOME_DIR}/reference/yleaf-${YLEAF_DATA_VERSION}/data"

echo "=== Y haplogroup (Yleaf): ${SAMPLE} ==="
if [ ! -s "${YDATA}/hg38/new_positions.txt" ]; then
  echo "ERROR: Yleaf's marker tables are not installed (${YDATA}): run scripts/setup.sh --yleaf-data ${GENOME_DIR}" >&2
  exit 1
fi
for f in "$BAM" "${BAM}.bai" "$REF_FASTA"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done
if [ ! -s "$PED" ]; then
  echo "ERROR: no sex check for ${SAMPLE} (${PED}): run scripts/16-indexcov.sh ${SAMPLE} first; this step reads the sex it infers." >&2
  exit 1
fi
# goleft's .ped: sex in PED coding (1 male, 2 female) under the header name 'sex'.
SEX=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == "sex") c = i }; next }
  c { s = $c } END { print (s == 1 ? "male" : s == 2 ? "female" : "unknown") }' "$PED")
if [ "$SEX" != male ]; then
  echo "Step 37 skipped: indexcov infers ${SEX} from the reads (${PED#"$S"/}); a Y haplogroup needs a Y chromosome."
  exit 0
fi

mkdir -p "$OUTDIR"
rm -rf "${OUTDIR}/yleaf" "$OUT"
# The Yleaf image has no samtools, and Yleaf downloads the whole hg38 FASTA
# unless its read-only config names one: bin/yleaf_run.py writes Yleaf's
# marker positions, SAMTOOLS_IMAGE makes the pileup Yleaf would make, and
# Yleaf then runs on it, offline (bin/yleaf_run.py says why and how).
O=$(cpath "$OUTDIR")
echo "[1/3] Yleaf's Y marker positions..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 2g "$YLEAF_IMAGE" \
  python3 /pgp-bin/yleaf_run.py positions --data "$(cpath "$YDATA")" "${O}/positions.txt"
echo "[2/3] Pileup at those positions (samtools)..."
run_in --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c \
  'samtools idxstats "$1" > "$2/idxstats.txt" && samtools mpileup -l "$2/positions.txt" -AQ20q1 "$1" > "$2/pileup.txt"' \
  _ "$(cpath "$BAM")" "$O"
echo "[3/3] Yleaf..."
echo "Reference: ${REF_FASTA}"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" --cpus 1 --memory 4g "$YLEAF_IMAGE" \
  python3 /pgp-bin/yleaf_run.py predict --data "$(cpath "$YDATA")" --bam "$(cpath "$BAM")" --reference "$REF_FASTA_C" \
    --idxstats "${O}/idxstats.txt" --pileup "${O}/pileup.txt" --out "${O}/yleaf"
rm -f "${OUTDIR}/idxstats.txt" "${OUTDIR}/pileup.txt"

PRED="${OUTDIR}/yleaf/hg_prediction.hg"
if [ ! -s "$PRED" ] || [ "$(grep -c . "$PRED")" -lt 2 ]; then
  echo "ERROR: Yleaf wrote no prediction (${PRED})." >&2
  exit 1
fi
cp "$PRED" "$OUT"
read -r HG MARKERS QC < <(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
  NR == 2 { print $c["Hg"], $c["Valid_markers"], $c["QC-score"] }' "$OUT")
echo ""
if [ "$HG" = NA ]; then
  echo "Y haplogroup: insufficient markers (Yleaf found ${MARKERS} Y markers with haplogroup information; none passes its quality threshold)"
else
  echo "Y haplogroup: ${HG} (${MARKERS} markers, QC-score ${QC})"
fi
echo "Results: ${OUT}"
