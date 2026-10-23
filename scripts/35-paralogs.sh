#!/usr/bin/env bash
# 35-paralogs.sh — [OPT-IN] SMN1 and SMN2 copy number with Parascopy, HBA1/HBA2 from depth
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
# HBA1 and HBA2 (alpha-thalassaemia deletions) have no Parascopy model, so the
# step counts them from read depth: the depth over both gene bodies against
# single-copy sequence on each side of the cluster gives the number of
# alpha-globin genes (bin/hba_copy_number.py). A failure there writes an HBA
# table that says not assessed, and does not stop the step.
#
# OPT-IN: a default run leaves this step out (run-all.sh runs it with
# TOOLS=...,parascopy).
#
# Environment:
#   PARASCOPY_POPULATION  model of AFR, AMR, EAS, EUR or SAS (default EUR)
#   PARASCOPY_DEPTH_BED   background windows (one size) for a BAM that covers
#                         part of the genome; default Parascopy's GRCh38 set.
#                         With it the depth is not stratified by GC content
#                         (--no-gc): a few regions rarely span the GC range
#                         Parascopy's GC model needs, and it stops on them.
#
# Output: ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_smn_copy_number.tsv
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_hba_copy_number.tsv
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
  BACKGROUND=(-b "$(cpath "$PARASCOPY_DEPTH_BED")" --no-gc)
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

# HBA1/HBA2 from depth (as the PARASCOPY module): samtools bedcov over the
# regions bin/hba_copy_number.py prints, all MAPQ, then its verdict.
HBA_OUT="${OUTDIR}/${SAMPLE}_hba_copy_number.tsv"
HBA_WORK="${OUTDIR}/hba_depth.part"
rm -rf "$HBA_OUT" "$HBA_WORK"
mkdir -p "$HBA_WORK"
if run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PARASCOPY_IMAGE}" \
     python3 /pgp-bin/hba_copy_number.py bed > "${HBA_WORK}/regions.bed" \
   && run_in --cpus 2 --memory 2g "${PARASCOPY_IMAGE}" \
     samtools bedcov "$(cpath "${HBA_WORK}/regions.bed")" "$(cpath "$BAM")" > "${HBA_WORK}/bedcov.tsv" \
   && run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PARASCOPY_IMAGE}" \
     python3 /pgp-bin/hba_copy_number.py check --bedcov "$(cpath "${HBA_WORK}/bedcov.tsv")" \
       --out "$(cpath "${HBA_OUT}.tmp")"; then
  mv "${HBA_OUT}.tmp" "$HBA_OUT"
else
  rm -f "${HBA_OUT}.tmp"
  echo "WARNING: the HBA1/HBA2 depth check failed; the report lists HBA1/HBA2 as not assessed." >&2
  printf 'metric\tvalue\nstatus\tnot assessed\nalpha_globin_genes\t.\nestimate\t.\nhba2_depth\t.\nhba1_depth\t.\nflank_left_depth\t.\nflank_right_depth\t.\nmessage\t%s\n' \
    'not assessed: the HBA1/HBA2 depth check failed; see the log of step 35' > "$HBA_OUT"
fi
rm -rf "$HBA_WORK"

echo ""
echo "=== Parascopy complete ==="
echo "Copy number: ${OUT}"
echo "HBA1/HBA2:   ${HBA_OUT}"
column -t -s $'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo ""
echo "SMN1 sits at chr5:70.92-70.95 Mb and SMN2 at chr5:70.05-70.08 Mb (GRCh38). psCN lists the copy"
echo "number of the row's region first, then of each region in homologous_regions. A quality under 20"
echo "or a filter other than PASS means the value is not reliable. HBA1/HBA2 is a count of"
echo "alpha-globin genes from depth, 4 without a deletion. See docs/35-paralogs.md."
