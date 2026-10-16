#!/usr/bin/env bash
# 35b-paraphase.sh — [OPT-IN] paralog genes from PacBio HiFi reads with Paraphase
# Usage: ./scripts/35b-paraphase.sh <sample_name>
#
# Near-identical gene copies (SMN1/SMN2, HBA1/HBA2, CYP21A2 in the RCCX
# module, GBA1/GBAP1, PMS2/PMS2CL, STRC/STRCP1, NCF1, OTOA and the rest of
# Paraphase's 162 GRCh38 regions) hide copy-number changes and variants from
# a normal alignment. Paraphase gathers the reads of each gene family,
# phases them into haplotypes and calls each gene's copy number, the HBA
# 3.7 kb and 4.2 kb deletions by name, and phased variants.
#
# PacBio HiFi only. Paraphase is written for HiFi reads, and its licence (the
# PacBio Software License Agreement) allows it only on data from a PacBio
# instrument. The step reads the BAM header first and stops unless every read
# group says PL:PACBIO (bin/paraphase_report.py check-bam). The default BAM is
# the one step 02b writes with PLATFORM=hifi. For short reads, step 35
# (Parascopy) covers SMN1 and SMN2.
#
# OPT-IN: a default run leaves this step out (run-all.sh runs it with
# TOOLS=...,paraphase, on the BAM it is given).
#
# Environment:
#   ALIGN_DIR          folder of <sample>_sorted.bam (default aligned_longread)
#   PARAPHASE_REGIONS  Paraphase regions, comma separated (default: all 162),
#                      e.g. smn1,hba,rccx,GBA,pms2,strc,ncf1,OTOA
#
# Output: ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_paraphase_copy_number.tsv
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_paraphase_variants.tsv
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_paraphase/ (Paraphase's JSON, VCFs, BAM)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned_longread}
BAM="${GENOME_DIR}/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
OUTDIR="${GENOME_DIR}/${SAMPLE}/paralogs"
REGIONS=${PARAPHASE_REGIONS:-}

echo "=== Paraphase paralog genes: ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done
if [[ -n "$REGIONS" && ! "$REGIONS" =~ ^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$ ]]; then
  echo "ERROR: PARAPHASE_REGIONS must be region names separated by commas (e.g. smn1,hba), got '${REGIONS}'." >&2
  exit 1
fi
echo "Input: ${BAM}"
echo "Regions: ${REGIONS:-all GRCh38 regions of Paraphase}"

# [1/3] PacBio reads only: the header's read groups must all say PL:PACBIO
echo "[1/3] Checking the read groups..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/paraphase_report.py check-bam "$(cpath "$BAM")"

# [2/3] Paraphase, into a temporary folder renamed when it finishes
echo "[2/3] Running Paraphase..."
mkdir -p "$OUTDIR"
WORK="${OUTDIR}/${SAMPLE}_paraphase.part"
FINAL="${OUTDIR}/${SAMPLE}_paraphase"
rm -rf "$WORK"
run_in --cpus "${THREADS}" --memory 8g "${PARAPHASE_IMAGE}" \
  paraphase -b "$(cpath "$BAM")" -r "$REF_FASTA_C" -o "$(cpath "$WORK")" -p "$SAMPLE" \
    -t "${THREADS}" ${REGIONS:+-g "$REGIONS"}

# [3/3] One row per gene and one per variant (bin/paraphase_report.py, as the PARAPHASE_REPORT module).
# The tables are written from WORK before it replaces FINAL, so a failed report leaves the
# previous results and their tables together.
echo "[3/3] Writing the tables..."
CN="${OUTDIR}/${SAMPLE}_paraphase_copy_number.tsv"
VARS="${OUTDIR}/${SAMPLE}_paraphase_variants.tsv"
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/paraphase_report.py table --sample "$SAMPLE" \
    --json "$(cpath "${WORK}/${SAMPLE}.paraphase.json")" \
    --vcf-dir "$(cpath "${WORK}/${SAMPLE}_paraphase_vcfs")" \
    --out-cn "$(cpath "$CN")" --out-variants "$(cpath "$VARS")"
rm -rf "$FINAL"
mv "$WORK" "$FINAL"

echo ""
echo "=== Paraphase complete ==="
echo "Copy number: ${CN}"
echo "Variants:    ${VARS}"
column -t -s $'\t' "$CN" 2>/dev/null | cut -c1-160 || cat "$CN"
echo ""
echo "A 'no call' means Paraphase saw one haplotype and the depth did not settle one copy or two."
echo "This is a research estimate, not a diagnostic test. See docs/35-paralogs.md."
