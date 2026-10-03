#!/usr/bin/env bash
# Mitochondrial analysis — Heteroplasmy detection + variant calling with GATK Mutect2
# Input: Sorted BAM
# Output: Mitochondrial VCF with heteroplasmy fractions
# Runtime: ~15-30 minutes
# Note: Originally planned for MToolBox, but no working Docker image exists.
#       GATK Mutect2 in mitochondrial mode is the standard alternative; the file
#       name is historical.
#
# NuMTs (nuclear copies of mitochondrial DNA) put reads on chrM that come from
# the autosomes. GATK's NuMTFilterTool marks an allele as possible_numt when
# its depth is what such copies could explain at the sample's median autosomal
# coverage. That median is read from step 16b's mosdepth output, or taken from
# AUTOSOMAL_COVERAGE; without either the filter marks nothing.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/mito"

echo "=== Mitochondrial Analysis (GATK Mutect2): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}"

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"


echo "[1/5] Extracting chrM reads..."
run_in --cpus 2 --memory 4g \
  "$SAMTOOLS_IMAGE" \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM \
      > /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam && \
    samtools index /genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam
  "

echo "[2/5] Checking sequence dictionary..."
if [ ! -f "$REF_DICT" ]; then
  echo "  Creating sequence dictionary..."
  # The dictionary goes next to the FASTA, so its directory is writable here.
  run_in --rw "$(dirname "$REF_FASTA")" --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk CreateSequenceDictionary \
      -R "${REF_FASTA_C}" \
      -O "$(cpath "$REF_DICT")"
else
  echo "  Sequence dictionary already exists, skipping."
fi

echo "[3/5] Running Mutect2 in mitochondrial mode..."
run_in --cpus 4 --memory 8g \
  "$GATK_IMAGE" \
  gatk Mutect2 \
    -R "${REF_FASTA_C}" \
    -I "/genome/${SAMPLE}/mito/${SAMPLE}_chrM.bam" \
    -L chrM \
    --mitochondria-mode \
    --max-mnp-distance 0 \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz"

echo "[4/5] Filtering variants..."
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk FilterMutectCalls \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2.vcf.gz" \
    --mitochondria-mode \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz"

echo "[5/5] Marking possible NuMTs..."
# Median autosomal depth: the depth at which half of the chr1-22 bases are
# covered at least that deep, from mosdepth's per-chromosome distribution.
MOSDEPTH_PREFIX="${SAMPLE_DIR}/mosdepth/${SAMPLE}.mosdepth"
if [ -n "${AUTOSOMAL_COVERAGE:-}" ]; then
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (AUTOSOMAL_COVERAGE)"
elif [ -s "${MOSDEPTH_PREFIX}.summary.txt" ] && [ -s "${MOSDEPTH_PREFIX}.global.dist.txt" ]; then
  AUTOSOMAL_COVERAGE=$(awk -F'\t' '
    FNR == NR { if ($1 ~ /^chr[0-9]+$/) len[$1] = $2; next }
    ($1 in len) { at_least[$2] += len[$1] * $3 }
    END {
      for (c in len) total += len[c]
      best = 0
      if (total > 0) for (k in at_least) if (at_least[k] / total >= 0.5 && k + 0 > best) best = k + 0
      print best
    }' "${MOSDEPTH_PREFIX}.summary.txt" "${MOSDEPTH_PREFIX}.global.dist.txt")
  echo "  Median autosomal coverage: ${AUTOSOMAL_COVERAGE} (from ${MOSDEPTH_PREFIX}.global.dist.txt)"
else
  AUTOSOMAL_COVERAGE=0
  echo "  WARNING: no mosdepth output (run scripts/16b-mosdepth.sh first) and no AUTOSOMAL_COVERAGE:"
  echo "  the NuMT filter runs with coverage 0 and marks nothing."
fi
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk NuMTFilterTool \
    -R "${REF_FASTA_C}" \
    -V "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_mutect2_filtered.vcf.gz" \
    --autosomal-coverage "$AUTOSOMAL_COVERAGE" \
    -O "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz"

echo "=== Mitochondrial analysis complete ==="
echo "Raw calls: ${OUTPUT_DIR}/${SAMPLE}_chrM_mutect2.vcf.gz"
echo "Filtered (FilterMutectCalls, then possible_numt): ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz"
echo ""
echo "PASS variants below AF 0.95 (possible heteroplasmy; calls near the ends of the"
echo "control region, chrM:1-500 and 16000-16569, are less reliable):"
echo "  bcftools query -f '%POS\t%REF\t%ALT\t[%AF]\n' -i 'FILTER=\"PASS\"' ${OUTPUT_DIR}/${SAMPLE}_chrM_filtered.vcf.gz | awk '\$4 < 0.95'"
