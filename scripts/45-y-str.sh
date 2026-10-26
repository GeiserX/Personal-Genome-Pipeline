#!/usr/bin/env bash
# 45-y-str.sh — Y-STR repeat counts (FamilyTreeDNA's Y-111 markers) with HipSTR
# Usage: ./scripts/45-y-str.sh <sample_name>
#
# Opt-in. HipSTR (HIPSTR_IMAGE) genotypes the Y-STR regions of
# assets/ystr_panel.tsv on chrY from the BAM, and bin/ystr_report.py, run with
# the image's own Python, turns its calls into repeat counts in
# FamilyTreeDNA's marker order (Y-12 to Y-111), with a confidence per marker.
# Multi-copy markers (DYS385, DYS459, DYS464, YCAII, CDY, DYF395S1, DYS413)
# and repeats longer than 150 bp reads can span are listed as not callable.
# The same commands as the Y_STR Nextflow process.
# Runs only when step 16 (indexcov) infers male from the reads: the sex comes
# from indexcov/indexcov-indexcov.ped, as for step 37. A female sample, or one
# with no Y to read, is skipped with one line.
#
# Output: y_str/<sample>_y_str.tsv     one row per marker: count, confidence, reads, note
#         y_str/<sample>_y_str.csv     marker,value,confidence in FamilyTreeDNA's order
#         y_str/<sample>_y_str.vcf.gz  HipSTR's calls; y_str/<sample>_y_str.log its log
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image HIPSTR_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
PED="${S}/indexcov/indexcov-indexcov.ped"
OUTDIR="${S}/y_str"
TSV="${OUTDIR}/${SAMPLE}_y_str.tsv"

echo "=== Y-STR markers (HipSTR): ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai"; do
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
  echo "Step 45 skipped: indexcov infers ${SEX} from the reads (${PED#"$S"/}); Y-STRs need a Y chromosome."
  # A table from an earlier run of this sample would reach the reports as current.
  rm -f "$TSV" "${OUTDIR}/${SAMPLE}_y_str."{csv,vcf.gz,log}
  exit 0
fi

echo "BAM: ${BAM}"
echo "Reference: ${REF_FASTA}"
mkdir -p "$OUTDIR"
rm -f "$TSV" "${OUTDIR}/${SAMPLE}_y_str."{csv,vcf.gz,log} "${OUTDIR}/regions.bed"
O=$(cpath "$OUTDIR")
# One container: the region file from the panel, HipSTR, then the table. The
# flags and why: bin/ystr_report.py and docs/45-y-str.md.
# shellcheck disable=SC2016  # $1 to $5 belong to the inner bash
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$HIPSTR_IMAGE" \
  bash -euo pipefail -c '
    python3 /pgp-bin/ystr_report.py bed --panel /pgp-assets/ystr_panel.tsv --out "$1/regions.bed"
    HipSTR --bams "$2" --fasta "$3" --regions "$1/regions.bed" --str-vcf "$1/$4_y_str.vcf.gz" \
      --bam-samps "$4" --bam-libs "$4" --haploid-chrs chrY --def-stutter-model --min-reads 5 \
      --max-str-len 150 --output-filters --log "$1/$4_y_str.log"
    python3 /pgp-bin/ystr_report.py report --panel /pgp-assets/ystr_panel.tsv --vcf "$1/$4_y_str.vcf.gz" \
      --log "$1/$4_y_str.log" --tsv "$1/$4_y_str.tsv" --csv "$1/$4_y_str.csv"' \
  _ "$O" "$(cpath "$BAM")" "$REF_FASTA_C" "$SAMPLE"
rm -f "${OUTDIR}/regions.bed"

[ -s "$TSV" ] || { echo "ERROR: no table written (${TSV})" >&2; exit 1; }
echo ""
echo "Table: ${TSV}"
echo "CSV in FamilyTreeDNA's order: ${OUTDIR}/${SAMPLE}_y_str.csv"
