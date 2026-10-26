#!/usr/bin/env bash
# 35-paralogs.sh — [OPT-IN] paralog copy number with Parascopy (SMN1 and SMN2 first), HBA1/HBA2 from depth
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
# The loci come from assets/paralog_loci.tsv: one row per Parascopy model of
# that data, run with one `parascopy cn-using` each. A row is not assessed,
# with the reason in its table, when its model is not installed, when the BAM
# has no read over the row's region, or when Parascopy fails on it; the other
# rows still run.
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
#   PARALOG_LOCI          the catalog (default assets/paralog_loci.tsv)
#
# Output: ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_<table>_copy_number.tsv
#           (one per catalog row; the SMN1 row's is ${SAMPLE}_smn_copy_number.tsv)
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_hba_copy_number.tsv
#         ${GENOME_DIR}/${SAMPLE}/paralogs/${SAMPLE}_parascopy/<table>/ (Parascopy's own files)
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
MODELS="${DATA}/models_GRCh38_1KGP/${POP}"
TABLE="${DATA}/homology_table/GRCh38.bed.gz"
LOCI=${PARALOG_LOCI:-${PGP_ROOT}/assets/paralog_loci.tsv}

echo "=== Parascopy paralog copy number: ${SAMPLE} ==="
for f in "$BAM" "${BAM}.bai" "$REF_FASTA" "${REF_FASTA}.fai" "$LOCI"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done
if [ ! -f "$TABLE" ] || [ ! -d "$MODELS" ]; then
  echo "ERROR: Parascopy ${PARASCOPY_DATA_VERSION}'s data is not installed (${TABLE}, ${MODELS}/)." >&2
  echo "  Install it once (~50 MB): ./scripts/setup.sh --parascopy-data ${GENOME_DIR}" >&2
  echo "  Population models: AFR AMR EAS EUR SAS (PARASCOPY_POPULATION, now ${POP})." >&2
  exit 1
fi
# One spec per catalog row, model|region|table (assets/paralog_loci.tsv; the
# Nextflow module reads the same columns). They reach file names and a shell
# loop, so a field with any other character stops the step.
SPECS=$(awk -F'\t' '/^#/ || !NF {next}
  !h {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
  {print $c["model"] "|" $c["region"] "|" $c["table"]}' "$LOCI")
if [ -z "$SPECS" ]; then
  echo "ERROR: ${LOCI} has no locus row." >&2
  exit 1
fi
while IFS= read -r spec; do
  if ! [[ "$spec" =~ ^[A-Za-z0-9._-]+\.gz\|chr[0-9XY]+:[0-9]+-[0-9]+\|[A-Za-z0-9._-]+$ ]]; then
    echo "ERROR: ${LOCI}: a row reads '${spec}' (model|region|table); see tests/test_paralog_loci.py." >&2
    exit 1
  fi
done <<<"$SPECS"
BACKGROUND=(-g GRCh38)
if [ -n "${PARASCOPY_DEPTH_BED:-}" ]; then
  BACKGROUND=(-b "$(cpath "$PARASCOPY_DEPTH_BED")" --no-gc)
  echo "Background windows: ${PARASCOPY_DEPTH_BED}"
fi
echo "Models: ${POP} (${MODELS}), loci: ${LOCI}"

mkdir -p "$OUTDIR"
WORK="${OUTDIR}/${SAMPLE}_parascopy.part"
rm -rf "$WORK" "${OUTDIR}/depth.part"
mkdir -p "$WORK"
# Background depth once, then the copy number of each catalog row
run_in --cpus "${THREADS}" --memory 8g "${PARASCOPY_IMAGE}" \
  parascopy depth -i "$(cpath "$BAM")::${SAMPLE}" -f "$REF_FASTA_C" "${BACKGROUND[@]}" \
    -o "$(cpath "${OUTDIR}/depth.part")" -@ "${THREADS}"
header() { printf 'chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n'; }
TABLES=()
# shellcheck disable=SC2086  # one word per row, checked above
for spec in $SPECS; do
  model=${spec%%|*} rest=${spec#*|}
  region=${rest%%|*} name=${rest#*|}
  out="${WORK}/${SAMPLE}_${name}_copy_number.tsv"
  TABLES+=("${SAMPLE}_${name}_copy_number.tsv")
  # A row is not assessed when its model is not installed, when the BAM has
  # no read over its region (a BAM of part of the genome), or when Parascopy
  # fails on it; the other rows still run.
  reason="" reads=0
  if [ -f "${MODELS}/${model}" ]; then
    # primary alignments only: no unmapped (0x4), secondary or supplementary
    reads=$(run_in --cpus 2 --memory 2g "${PARASCOPY_IMAGE}" \
      samtools view -c -F 0x904 "$(cpath "$BAM")" "$region") || reads=0
  fi
  if [ ! -f "${MODELS}/${model}" ]; then
    reason="no model ${model} in the installed Parascopy data (${POP})"
  elif ! [[ "$reads" =~ ^[0-9]+$ ]] || [ "$reads" -eq 0 ]; then
    reason="no read over ${region}"
  elif run_in --cpus "${THREADS}" --memory 8g "${PARASCOPY_IMAGE}" \
         parascopy cn-using "$(cpath "${MODELS}/${model}")" \
           -i "$(cpath "$BAM")::${SAMPLE}" \
           -f "$REF_FASTA_C" \
           -t "$(cpath "$TABLE")" \
           -d "$(cpath "${OUTDIR}/depth.part")" \
           -o "$(cpath "${WORK}/${name}")" \
           -@ "${THREADS}"; then
    # One row per region of the locus profile: aggregate and paralog-specific
    # copy number, each with its filter and quality (as the PARASCOPY module).
    {
      header
      gzip -dc "${WORK}/${name}/res.samples.bed.gz" \
        | awk -F'\t' -v OFS='\t' '!/^#/ {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11, $13}'
    } > "$out"
    continue
  else
    reason="parascopy cn-using failed; see the log of step 35"
  fi
  echo "WARNING: ${model}: not assessed: ${reason}" >&2
  rm -rf "${WORK:?}/${name}"
  { header; printf '.\t.\t.\t%s\tnot assessed: %s\t.\t.\t.\t.\t.\t.\n' "${model%.gz}" "$reason"; } > "$out"
done
for t in "${TABLES[@]}"; do mv "${WORK}/${t}" "${OUTDIR}/${t}"; done
rm -rf "${OUTDIR}/depth.part" "${OUTDIR}/${SAMPLE}_parascopy"
mv "$WORK" "${OUTDIR}/${SAMPLE}_parascopy"

# HBA1/HBA2 from depth (as the PARASCOPY module): samtools bedcov over the
# regions bin/hba_copy_number.py prints, all MAPQ and MAPQ 1 or more (a BAM
# with no MAPQ 0 read over the genes left out the reads that fit both), then
# its verdict.
HBA_OUT="${OUTDIR}/${SAMPLE}_hba_copy_number.tsv"
HBA_WORK="${OUTDIR}/hba_depth.part"
rm -rf "$HBA_OUT" "$HBA_WORK"
mkdir -p "$HBA_WORK"
if run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PARASCOPY_IMAGE}" \
     python3 /pgp-bin/hba_copy_number.py bed > "${HBA_WORK}/regions.bed" \
   && run_in --cpus 2 --memory 2g "${PARASCOPY_IMAGE}" \
     samtools bedcov "$(cpath "${HBA_WORK}/regions.bed")" "$(cpath "$BAM")" > "${HBA_WORK}/bedcov.tsv" \
   && run_in --cpus 2 --memory 2g "${PARASCOPY_IMAGE}" \
     samtools bedcov -Q 1 "$(cpath "${HBA_WORK}/regions.bed")" "$(cpath "$BAM")" > "${HBA_WORK}/bedcov_q1.tsv" \
   && run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PARASCOPY_IMAGE}" \
     python3 /pgp-bin/hba_copy_number.py check --bedcov "$(cpath "${HBA_WORK}/bedcov.tsv")" \
       --bedcov-mapq1 "$(cpath "${HBA_WORK}/bedcov_q1.tsv")" --out "$(cpath "${HBA_OUT}.tmp")"; then
  mv "${HBA_OUT}.tmp" "$HBA_OUT"
else
  rm -f "${HBA_OUT}.tmp"
  echo "WARNING: the HBA1/HBA2 depth check failed; the report lists HBA1/HBA2 as not assessed." >&2
  printf 'metric\tvalue\nstatus\tnot assessed\nalpha_globin_genes\t.\nestimate\t.\nhba2_depth\t.\nhba1_depth\t.\nflank_left_depth\t.\nflank_right_depth\t.\nmapq0_share\t.\nmessage\t%s\n' \
    'not assessed: the HBA1/HBA2 depth check failed; see the log of step 35' > "$HBA_OUT"
fi
rm -rf "$HBA_WORK"

echo ""
echo "=== Parascopy complete ==="
echo "HBA1/HBA2:   ${HBA_OUT}"
for t in "${TABLES[@]}"; do
  echo "Copy number: ${OUTDIR}/${t}"
  column -t -s $'\t' "${OUTDIR}/${t}" 2>/dev/null || cat "${OUTDIR}/${t}"
done
echo ""
echo "A row that says not assessed names the reason. SMN1 sits at chr5:70.92-70.95 Mb and SMN2 at chr5:70.05-70.08 Mb (GRCh38). psCN lists the copy"
echo "number of the row's region first, then of each region in homologous_regions. A quality under 20"
echo "or a filter other than PASS means the value is not reliable. HBA1/HBA2 is a count of"
echo "alpha-globin genes from depth, 4 without a deletion. See docs/35-paralogs.md."
