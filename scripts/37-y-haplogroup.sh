#!/usr/bin/env bash
# 37-y-haplogroup.sh — Y-chromosome haplogroup (paternal line) with Yleaf
# Usage: ./scripts/37-y-haplogroup.sh <sample_name>
#
# Opt-in. Yleaf (YLEAF_IMAGE, 3.2.1: the only biocontainer; upstream is 4.x)
# reads the BAM's pileup at its Y-chromosome markers and predicts the
# haplogroup they support. Runs only when step 16 (indexcov) infers male from
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
require_image YLEAF_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
PED="${S}/indexcov/indexcov-indexcov.ped"
OUTDIR="${S}/y_haplogroup"
OUT="${OUTDIR}/${SAMPLE}_y_haplogroup.txt"

echo "=== Y haplogroup (Yleaf): ${SAMPLE} ==="
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
# Yleaf downloads the whole hg38 FASTA on its first run unless its config
# names one, and the image's config is read-only: its constant is pointed at
# the reference before it starts (a BAM never needs the sequence itself).
# Its multiprocessing pools are replaced by a serial map: a failing samtools
# call inside a pool worker raises SystemExit, the worker dies and the pool
# waits for its result forever; serial, Yleaf stops with the error.
echo "Reference: ${REF_FASTA}"
run_in --cpus 1 --memory 4g "$YLEAF_IMAGE" python3 -c 'import sys, multiprocessing; from pathlib import Path; multiprocessing.Pool = type("SerialPool", (), {"__init__": lambda s, *a, **k: None, "__enter__": lambda s: s, "__exit__": lambda s, *a: False, "map": lambda s, f, xs: list(map(f, xs))}); from yleaf import yleaf_constants; yleaf_constants.HG38_FULL_GENOME = Path(sys.argv[1]); from yleaf import Yleaf; sys.argv = ["Yleaf"] + sys.argv[2:]; Yleaf.main()' "$REF_FASTA_C" -bam "$(cpath "$BAM")" -o "$(cpath "${OUTDIR}/yleaf")" -rg hg38 -force -t 1

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
