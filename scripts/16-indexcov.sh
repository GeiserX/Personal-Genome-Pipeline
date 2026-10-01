#!/usr/bin/env bash
# indexcov (goleft) — Rapid coverage QC and sex chromosome check from BAM index
# Input: sorted BAM with .bai index
# Output: HTML report with per-chromosome coverage + sex check
# Instant (~5 seconds, reads only the BAM index)
#
# Usage: ./scripts/16-indexcov.sh <sample_name> [declared_sex: male|female]
# With a declared sex, the step exits non-zero when the sex inferred from X/Y
# coverage differs (a sample swap or a sex-chromosome aneuploidy). SEX_CHECK=warn
# prints the mismatch and exits 0.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [declared_sex: male|female]}
DECLARED_SEX=${2:-}
SEX_CHECK=${SEX_CHECK:-fail}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
OUTPUT_DIR="${SAMPLE_DIR}/indexcov"

echo "=== indexcov: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

case "$DECLARED_SEX" in
  ""|male|female) ;;
  *) echo "ERROR: declared sex must be 'male' or 'female', got '${DECLARED_SEX}'" >&2; exit 1 ;;
esac

mkdir -p "$OUTPUT_DIR"
PED="${OUTPUT_DIR}/indexcov-indexcov.ped"
# Never read a .ped left by an earlier run
rm -f "$PED"

docker run --rm \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/goleft:0.2.6--he881be0_1 \
  goleft indexcov \
    --directory "/genome/${SAMPLE}/indexcov" \
    "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

echo "=== indexcov complete ==="
echo "Results: ${OUTPUT_DIR}/"
echo ""
echo "Key outputs:"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.html  — interactive coverage plot"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.ped   — sex chromosome inference"
echo "  ${OUTPUT_DIR}/indexcov-indexcov.roc   — coverage uniformity"
echo ""
if [ ! -f "$PED" ]; then
  echo "ERROR: goleft wrote no ${PED}" >&2
  exit 1
fi

# goleft's .ped columns: #family_id sample_id paternal_id maternal_id sex phenotype CNchrX CNchrY ...
# sex is PED coding (1 male, 2 female, anything else unknown); phenotype is always -9.
# Columns are found by header name, so a reordered file cannot be misread.
echo "Sex check (from X/Y coverage):"
echo "  .ped row: $(grep -v '^#' "$PED" | tail -1)"
read -r SEX_CODE CN_X CN_Y < <(awk '
  /^#/ { for(i=1;i<=NF;i++) { h=$i; sub(/^#/,"",h); col[h]=i } next }
  { split($0,f) }
  END {
    if(!("sex" in col)) { print "NOCOL NA NA"; exit }
    x=("CNchrX" in col) ? f[col["CNchrX"]] : "NA"
    y=("CNchrY" in col) ? f[col["CNchrY"]] : "NA"
    print f[col["sex"]], x, y
  }' "$PED")
if [ "$SEX_CODE" = "NOCOL" ]; then
  echo "ERROR: ${PED} has no 'sex' column in its header" >&2
  exit 1
fi
case "$SEX_CODE" in
  1) INFERRED_SEX=male ;;
  2) INFERRED_SEX=female ;;
  *) INFERRED_SEX=unknown ;;
esac
echo "  CNchrX: ${CN_X}  CNchrY: ${CN_Y}"
echo "  Predicted sex: ${INFERRED_SEX}"

if [ -n "$DECLARED_SEX" ]; then
  echo "  Declared sex:  ${DECLARED_SEX}"
  if [ "$INFERRED_SEX" = "$DECLARED_SEX" ]; then
    echo "  Sex check: OK"
  else
    echo "" >&2
    echo "!!! SEX CHECK MISMATCH: declared ${DECLARED_SEX}, inferred ${INFERRED_SEX} (CNchrX=${CN_X}, CNchrY=${CN_Y})" >&2
    echo "!!! Either the sample is not who you think it is, the declared sex is wrong, or the" >&2
    echo "!!! sample has a sex-chromosome aneuploidy. Steps that use the sex (ExpansionHunter) would get the wrong value." >&2
    if [ "$SEX_CHECK" = "warn" ]; then
      echo "!!! SEX_CHECK=warn: continuing." >&2
    else
      echo "!!! Set SEX_CHECK=warn to continue anyway." >&2
      exit 1
    fi
  fi
fi
