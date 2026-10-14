#!/usr/bin/env bash
# 41-bloodtype.sh — Blood group antigens (red cell and platelet) predicted from the genome
# Usage: ./scripts/41-bloodtype.sh <sample_name>
#
# Opt-in. Reads the sample's genotype at the defining variants of ABO (A1, A2,
# B, O), RhCE, Kell, Duffy (with the GATA box), Kidd, MNS, Lutheran, Diego,
# Dombrock, Colton and the platelet antigens HPA-1 to -5 and -15
# (assets/bloodtype_markers.tsv), and RhD from read depth: the common
# RhD-negative allele deletes the whole RHD gene, which a VCF does not show.
# bin/bloodtype_report.py applies the rules; it is the same code as the
# BLOODTYPE Nextflow processes. A prediction from the DNA, not a serological
# typing.
#
# Input:  vcf/<sample>.vcf.gz; aligned/<sample>_sorted.bam with its .bai
#         (ALIGN_DIR as for the other BAM steps), optional: without it RhD is
#         not determined and sites the VCF has no record for are read as
#         reference with confidence medium.
# Output: bloodtype/<sample>_bloodtype.tsv           one row per system
#         bloodtype/<sample>_bloodtype_evidence.tsv  what the VCF and the reads said per site
#         bloodtype/<sample>_bloodtype_depth.tsv     samtools bedcov over the table (with a BAM)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
require_image PYTHON_IMAGE SAMTOOLS_IMAGE

S="${GENOME_DIR}/${SAMPLE}"
VCF="${S}/vcf/${SAMPLE}.vcf.gz"
BAM="${S}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
OUTDIR="${S}/bloodtype"
TABLE="${OUTDIR}/${SAMPLE}_bloodtype.tsv"
EVIDENCE="${OUTDIR}/${SAMPLE}_bloodtype_evidence.tsv"
DEPTH="${OUTDIR}/${SAMPLE}_bloodtype_depth.tsv"

echo "=== Blood groups (predicted from the genome): ${SAMPLE} ==="
[ -f "$VCF" ] || { echo "ERROR: File not found: ${VCF}" >&2; exit 1; }
echo "VCF: ${VCF}"
mkdir -p "$OUTDIR"
rm -f "$TABLE" "$EVIDENCE" "$DEPTH" "${OUTDIR}/targets.bed"
O=$(cpath "$OUTDIR")
ARGS=()
if [ -f "$BAM" ] && [ -f "${BAM}.bai" ]; then
  echo "BAM: ${BAM}"
  echo "[1/2] Read depth at the marker sites and over RHD (samtools bedcov, MAPQ >= 20)..."
  # The marker table is a BED with extra columns; bedcov reads it once the
  # comment and header lines are gone.
  run_in -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$SAMTOOLS_IMAGE" bash -euo pipefail -c \
    'grep -v "^#" "$1" | tail -n +2 > "$2/targets.bed" && samtools bedcov -Q 20 -j "$2/targets.bed" "$3" > "$4"' \
    _ /pgp-assets/bloodtype_markers.tsv "$O" "$(cpath "$BAM")" "$(cpath "$DEPTH")"
  rm -f "${OUTDIR}/targets.bed"
  ARGS+=(--depth "$(cpath "$DEPTH")")
else
  echo "[1/2] No BAM (${BAM}): RhD needs read depth and is reported as not determined."
fi
echo "[2/2] Blood group rules (bin/bloodtype_report.py)..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" -v "${PGP_ROOT}/assets:/pgp-assets:ro" --cpus 1 --memory 2g "$PYTHON_IMAGE" \
  python3 /pgp-bin/bloodtype_report.py --sample "$SAMPLE" --vcf "$(cpath "$VCF")" \
    --markers /pgp-assets/bloodtype_markers.tsv ${ARGS[@]+"${ARGS[@]}"} \
    --out "$(cpath "$TABLE")" --evidence "$(cpath "$EVIDENCE")"

[ -s "$TABLE" ] || { echo "ERROR: no table written (${TABLE})" >&2; exit 1; }
echo ""
echo "Table: ${TABLE}"
echo "Evidence per site: ${EVIDENCE}"
