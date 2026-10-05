#!/usr/bin/env bash
# 35-paralogs.sh — [OPT-IN] SMN1 and SMN2 copy number with Parascopy
# Usage: ./scripts/35-paralogs.sh <sample_name>
#
# SMN1 and SMN2 are near-identical copies on chr5: short reads cannot be
# placed on one of them, so no VCF-based step sees an SMN1 copy-number loss,
# the usual cause of SMA carrier status. Parascopy (MIT licence) estimates
# the copy number of the pair (agCN) and, where the sequence differences
# between the copies allow, of each copy (psCN), each with a Phred quality
# (20 = 99% likely right). It needs a BAM aligned to a reference without ALT
# contigs, the default, and the homology table and 1000 Genomes models that
# `setup.sh --parascopy-data` installs.
#
# OPT-IN: a default run leaves this step out (run-all.sh runs it with
# TOOLS=...,parascopy).
#
# Environment:
#   PARASCOPY_POPULATION  model of AFR, AMR, EAS, EUR or SAS (default EUR)
#   PARASCOPY_DEPTH_BED   background windows (one size) for a BAM that covers
#                         part of the genome; default Parascopy's GRCh38 set
#
# Output: ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_smn_copy_number.tsv
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_parascopy/ (Parascopy's own files)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${GENOME_DIR}/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
OUTDIR="${GENOME_DIR}/${SAMPLE}/paralogs"
DATA="${GENOME_DIR}/reference/parascopy-${PARASCOPY_DATA_VERSION}"
POP=${PARASCOPY_POPULATION:-EUR}
MODEL="${DATA}/models_GRCh38_1KGP/${POP}/SMN1.gz"
TABLE="${DATA}/homology_table/GRCh38.bed.gz"

echo "=== Parascopy SMN1/SMN2 copy number: ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done
if [ ! -f "$TABLE" ] || [ ! -f "$MODEL" ]; then
  echo "ERROR: Parascopy ${PARASCOPY_DATA_VERSION}'s data is not installed (${TABLE}, ${MODEL})." >&2
  echo "  Install it once (~50 MB): ./scripts/setup.sh --parascopy-data ${GENOME_DIR}" >&2
  echo "  Population models: AFR AMR EAS EUR SAS (PARASCOPY_POPULATION, now ${POP})." >&2
  exit 1
fi
BACKGROUND=(-g GRCh38)
if [ -n "${PARASCOPY_DEPTH_BED:-}" ]; then
  BACKGROUND=(-b "$(cpath "$PARASCOPY_DEPTH_BED")")
  echo "Background windows: ${PARASCOPY_DEPTH_BED}"
fi
echo "Model: ${POP} (${MODEL})"

mkdir -p "$OUTDIR"
WORK="${OUTDIR}/${SAMPLE}_parascopy.part"
rm -rf "$WORK" "${OUTDIR}/depth.part"
# Background depth, then the copy number of the SMN1/SMN2 locus
run_in --cpus "${THREADS}" --memory 8g "${PARASCOPY_IMAGE}" \
  parascopy depth -i "$(cpath "$BAM")::${SAMPLE}" -f "$REF_FASTA_C" "${BACKGROUND[@]}" \
    -o "$(cpath "${OUTDIR}/depth.part")" -@ "${THREADS}"
run_in --cpus "${THREADS}" --memory 8g "${PARASCOPY_IMAGE}" \
  parascopy cn-using "$(cpath "$MODEL")" \
    -i "$(cpath "$BAM")::${SAMPLE}" \
    -f "$REF_FASTA_C" \
    -t "$(cpath "$TABLE")" \
    -d "$(cpath "${OUTDIR}/depth.part")" \
    -o "$(cpath "$WORK")" \
    -@ "${THREADS}"
rm -rf "${OUTDIR}/depth.part" "${OUTDIR}/${SAMPLE}_parascopy"
mv "$WORK" "${OUTDIR}/${SAMPLE}_parascopy"

# One row per region of the SMN1/SMN2 profile: aggregate and paralog-specific
# copy number, each with its filter and quality (as the PARASCOPY module).
OUT="${OUTDIR}/${SAMPLE}_smn_copy_number.tsv"
{
  printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n'
  gzip -dc "${OUTDIR}/${SAMPLE}_parascopy/res.samples.bed.gz" \
    | awk -F'\t' -v OFS='\t' '!/^#/ {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11, $13}'
} > "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"

echo ""
echo "=== Parascopy complete ==="
echo "Copy number: ${OUT}"
column -t -s $'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo ""
echo "SMN1 sits at chr5:70.92-70.95 Mb and SMN2 at chr5:70.05-70.08 Mb (GRCh38). psCN lists the copy"
echo "number of the row's region first, then of each region in homologous_regions. A quality under 20"
echo "or a filter other than PASS means the value is not reliable. See docs/35-paralogs.md."
