#!/usr/bin/env bash
# ExpansionHunter v5 — Screen for pathogenic short tandem repeat (STR) expansions
# Detects: Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, SCAs, myotonic dystrophy, etc.
# Input: sorted BAM + GRCh38 reference FASTA
# Output: JSON + VCF with repeat genotypes in $GENOME_DIR/<sample>/expansion_hunter/
#
# Upgraded from v2.5.5 to v5.0.0:
#   - --bam → --reads, --ref-fasta → --reference, --repeat-specs → --variant-catalog
#   - --vcf/--json/--log → --output-prefix (auto-generates .vcf, .json)
#   - Multithreading support (--threads)
#   - Bundled GRCh38 variant catalog (31 pathogenic loci) inside the container
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <male|female>}
SEX=${2:?Usage: $0 <sample_name> <male|female>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/expansion_hunter"

echo "=== ExpansionHunter v5: ${SAMPLE} (${SEX}) ==="

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# ExpansionHunter v5.0.0 via biocontainer
# The variant catalog (31 pathogenic GRCh38 loci) is bundled inside the container
# at /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json
docker run --rm \
  --cpus "${THREADS}" --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/expansionhunter:5.0.0--hc26b3af_5 \
  ExpansionHunter \
    --reads "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    --reference /genome/reference/Homo_sapiens_assembly38.fasta \
    --variant-catalog /usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json \
    --output-prefix "/genome/${SAMPLE}/expansion_hunter/${SAMPLE}_eh" \
    --threads "${THREADS}" \
    --sex "$SEX"

echo "=== ExpansionHunter complete ==="
echo "VCF:  ${OUTPUT_DIR}/${SAMPLE}_eh.vcf"
echo "JSON: ${OUTPUT_DIR}/${SAMPLE}_eh.json"
echo ""
echo "Key disease thresholds (alleles >= threshold = pathogenic):"
echo "  HTT(Huntington): 36  FMR1(FragileX): 55  FXN(Friedreich): 66"
echo "  C9ORF72(ALS/FTD): 30  DMPK(MyotonicDystrophy): 50  ATXN1(SCA1): 39"
