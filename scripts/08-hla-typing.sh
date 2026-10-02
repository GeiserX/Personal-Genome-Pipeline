#!/usr/bin/env bash
# HLA Typing — T1K (Class I + II, 4-digit resolution)
# Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
# Uses IPD-IMGT/HLA database aligned against GRCh38 reference
#
# IMPORTANT: The coordinate file MUST be built with the actual reference FASTA,
# NOT the .fai index. Using .fai produces -1 coordinates and empty results.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
IDX_DIR="${GENOME_DIR}/t1k_idx"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/hla_t1k"

echo "=== T1K HLA Typing: ${SAMPLE} ==="

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# Step 1: Build HLA reference index (one-time, ~5 min)
if [ ! -f "${IDX_DIR}/hlaidx/_dna_seq.fa" ]; then
  echo "Building HLA reference index..."
  mkdir -p "${IDX_DIR}"
  # --net: t1k-build.pl downloads the IPD-IMGT/HLA database (hla.dat).
  # The index is shared by every sample, so t1k_idx/ is writable here.
  run_in --net --rw "$IDX_DIR" --cpus 2 --memory 2g \
    "${T1K_IMAGE}" \
    t1k-build.pl -o /genome/t1k_idx/hlaidx --download IPD-IMGT/HLA
fi

# Step 2: Build coordinate file (one-time, ~30 min — reads entire 3.1GB reference)
# CRITICAL: Use the actual FASTA file, NOT the .fai index!
if [ ! -f "${IDX_DIR}/hlaidx_grch38/_dna_coord.fa" ]; then
  echo "Building coordinate file from reference genome (this takes ~30 min)..."
  run_in --rw "$IDX_DIR" --cpus 4 --memory 8g \
    "${T1K_IMAGE}" \
    t1k-build.pl \
      -d "/genome/t1k_idx/hlaidx/hla.dat" \
      -g "${REF_FASTA_C}" \
      -o /genome/t1k_idx/hlaidx_grch38
fi

# Step 3: Run HLA typing
echo "Running T1K genotyping..."
run_in \
  --cpus 4 --memory 8g \
  "${T1K_IMAGE}" \
  run-t1k \
    -b "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    -f "/genome/t1k_idx/hlaidx_grch38/_dna_seq.fa" \
    -c "/genome/t1k_idx/hlaidx_grch38/_dna_coord.fa" \
    --preset hla-wgs \
    -t 4 \
    --od "/genome/${SAMPLE}/hla_t1k/" \
    -o "${SAMPLE}_hla"

echo "=== T1K complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_hla_genotype.tsv"
