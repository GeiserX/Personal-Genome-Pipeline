#!/usr/bin/env bash
# ExpansionHunter — Screen for pathogenic short tandem repeat (STR) expansions
# Detects: Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, SCAs, myotonic dystrophy, etc.
# Input: sorted BAM + GRCh38 reference FASTA
# Output: JSON + VCF with repeat genotypes
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <male|female>}
SEX=${2:?Usage: $0 <sample_name> <male|female>}
GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}
SAMPLE_DIR="${GENOMA_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOMA_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"

echo "=== ExpansionHunter: ${SAMPLE} (${SEX}) ==="
mkdir -p "$OUTPUT_DIR"

docker run --rm \
  --cpus 4 --memory 4g \
  -v "${GENOMA_DIR}:/genoma" \
  --entrypoint /ExpansionHunter/bin/ExpansionHunter \
  weisburd/expansionhunter:latest \
    --bam "/genoma/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    --ref-fasta /genoma/reference/Homo_sapiens_assembly38.fasta \
    --repeat-specs /pathogenic_repeats/GRCh38/ \
    --vcf "/genoma/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.vcf" \
    --json "/genoma/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.json" \
    --log "/genoma/${SAMPLE}/expansion_hunter/${SAMPLE}_eh.log" \
    --sex "$SEX"

echo "=== ExpansionHunter complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_eh.json"
echo ""
echo "Key disease thresholds (alleles >= threshold = pathogenic):"
echo "  HTT(Huntington): 36  FMR1(FragileX): 55  FXN(Friedreich): 66"
echo "  C9ORF72(ALS/FTD): 30  DMPK(MyotonicDystrophy): 50  ATXN1(SCA1): 39"
