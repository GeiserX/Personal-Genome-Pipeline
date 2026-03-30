#!/usr/bin/env bash
# HLA Typing — T1K (Class I + II, 4-digit resolution)
# Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
# Uses IPD-IMGT/HLA database aligned against GRCh38 reference
#
# IMPORTANT: The coordinate file MUST be built with the actual reference FASTA,
# NOT the .fai index. Using .fai produces -1 coordinates and empty results.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}
BAM="${GENOMA_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOMA_DIR}/reference/Homo_sapiens_assembly38.fasta"
IDX_DIR="${GENOMA_DIR}/t1k_idx"
OUTPUT_DIR="${GENOMA_DIR}/${SAMPLE}/hla_t1k"

echo "=== T1K HLA Typing: ${SAMPLE} ==="
mkdir -p "$OUTPUT_DIR"

# Step 1: Build HLA reference index (one-time, ~5 min)
if [ ! -f "${IDX_DIR}/hlaidx/_dna_seq.fa" ]; then
  echo "Building HLA reference index..."
  mkdir -p "${IDX_DIR}"
  docker run --rm --cpus 2 --memory 2g \
    -v "${IDX_DIR}:/idx" \
    quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
    t1k-build.pl -o /idx/hlaidx --download IPD-IMGT/HLA
fi

# Step 2: Build coordinate file (one-time, ~30 min — reads entire 3.1GB reference)
# CRITICAL: Use the actual FASTA file, NOT the .fai index!
if [ ! -f "${IDX_DIR}/hlaidx_grch38/_dna_coord.fa" ]; then
  echo "Building coordinate file from reference genome (this takes ~30 min)..."
  docker run --rm --cpus 4 --memory 8g \
    -v "${GENOMA_DIR}:/genoma" \
    quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
    t1k-build.pl \
      -d "/genoma/t1k_idx/hlaidx/hla.dat" \
      -g "/genoma/reference/Homo_sapiens_assembly38.fasta" \
      -o /genoma/t1k_idx/hlaidx_grch38
fi

# Step 3: Run HLA typing
echo "Running T1K genotyping..."
docker run --rm \
  --cpus 4 --memory 8g \
  -v "${GENOMA_DIR}:/genoma" \
  quay.io/biocontainers/t1k:1.0.9--h5ca1c30_0 \
  run-t1k \
    -b "/genoma/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    -f "/genoma/t1k_idx/hlaidx_grch38/_dna_seq.fa" \
    -c "/genoma/t1k_idx/hlaidx_grch38/_dna_coord.fa" \
    --preset hla-wgs \
    -t 4 \
    --od "/genoma/${SAMPLE}/hla_t1k/" \
    -o "${SAMPLE}_hla"

echo "=== T1K complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_hla_genotype.tsv"
