#!/usr/bin/env bash
# ExpansionHunter — Screen for pathogenic short tandem repeat (STR) expansions
# Detects: Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, SCAs, myotonic dystrophy, etc.
# Input: sorted BAM + GRCh38 reference FASTA
# Output: JSON + VCF with repeat genotypes
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <male|female>}
SEX=${2:?Usage: $0 <sample_name> <male|female>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"

echo "=== ExpansionHunter: ${SAMPLE} (${SEX}) ==="

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  --entrypoint /ExpansionHunter/bin/ExpansionHunter \
  weisburd/expansionhunter:latest \
    --bam "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    --ref-fasta /genome/reference/Homo_sapiens_assembly38.fasta \
    --repeat-specs /pathogenic_repeats/GRCh38/ \
    --vcf "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf" \
    --json "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.json" \
    --log "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.log" \
    --sex "$SEX"

echo "=== ExpansionHunter complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_eh.json"
echo ""
echo "Key disease thresholds (alleles >= threshold = pathogenic):"
echo "  HTT(Huntington): 36  FMR1(FragileX): 55  FXN(Friedreich): 66"
echo "  C9ORF72(ALS/FTD): 30  DMPK(MyotonicDystrophy): 50  ATXN1(SCA1): 39"
