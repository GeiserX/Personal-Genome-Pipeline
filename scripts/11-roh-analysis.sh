#!/usr/bin/env bash
# ROH Analysis — Runs of homozygosity (consanguinity screening)
# Detects long stretches of homozygous DNA that indicate shared ancestry
# Centromeric regions (chr1 125-143MB, chr9 42-60MB, chr18 15-20MB) are known artifacts
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt"

echo "=== ROH Analysis: ${SAMPLE} ==="

docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  staphb/bcftools:1.21 \
  bcftools roh --AF-dflt 0.4 -o "/data/${SAMPLE}_roh.txt" "/data/${SAMPLE}.vcf.gz"

echo "=== ROH complete ==="
echo "Results: ${OUTPUT}"
echo ""
echo "Autosomal ROH >5MB (potential consanguinity signal):"
grep '^RG' "$OUTPUT" 2>/dev/null | awk '$3 !~ /chrX|chrY/ && $6 > 5000000 {printf "%s:%s-%s  %.1fMB\n", $3,$4,$5,$6/1e6}' || true
echo ""
echo "NOTE: Centromeric ROH (chr1:125-143MB, chr9:42-60MB, chr18:15-20MB) are technical artifacts, not real."
