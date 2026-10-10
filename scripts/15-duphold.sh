#!/usr/bin/env bash
# duphold — Annotate structural variants with depth-based quality metrics
# Input: Manta diploidSV.vcf.gz + sorted BAM + reference FASTA
# Output: SV VCF with DHBFC/DHFFC annotations, and the same calls after the
#   depth filter (<sample>_sv_filtered.vcf.gz), which step 05 annotates
# Very fast (~20 minutes)
#
# The filter is the Nextflow DUPHOLD_FILTER's: a deletion stays only when the
# depth drop against its flanks is real (DHFFC < 0.7), a duplication only when
# the gain against GC-matched bins is real (DHBFC > 1.3); other types, and
# records without a value, stay.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
MANTA_VCF="${SAMPLE_DIR}/manta/results/variants/diploidSV.vcf.gz"
# Fall back to manta2/ if a second Manta run was used
[ ! -f "$MANTA_VCF" ] && MANTA_VCF="${SAMPLE_DIR}/manta2/results/variants/diploidSV.vcf.gz"
OUTPUT_DIR="${SAMPLE_DIR}/duphold"

echo "=== duphold: ${SAMPLE} ==="
echo "Input SV VCF: ${MANTA_VCF}"
echo "Input BAM: ${BAM}"
echo "Output: ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf"

# Validate inputs
for f in "$MANTA_VCF" "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

run_in \
  --cpus 4 --memory 4g \
  "${DUPHOLD_IMAGE}" \
  duphold \
    -v "/genome/${SAMPLE}/$(echo "$MANTA_VCF" | sed "s|${SAMPLE_DIR}/||")" \
    -b "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
    -f "${REF_FASTA_C}" \
    -o "/genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf"

# The expression of DUPHOLD_FILTER (modules/local/duphold/main.nf), applied
# through a temporary name so a failed run leaves no half-written file.
FILTER_EXPR='(INFO/SVTYPE="DEL" && FMT/DHFFC[0] >= 0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0] <= 1.3)'
ANNOTATED_C="/genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf"
FILTERED="${OUTPUT_DIR}/${SAMPLE}_sv_filtered.vcf.gz"
FILTERED_C="/genome/${SAMPLE}/duphold/${SAMPLE}_sv_filtered.vcf.gz"
rm -f "${FILTERED}.tbi" "${FILTERED}.tmp"
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" \
  bcftools view -e "$FILTER_EXPR" -Oz -o "${FILTERED_C}.tmp" "$ANNOTATED_C"
mv -f "${FILTERED}.tmp" "$FILTERED"
run_in --cpus 1 --memory 2g "$BCFTOOLS_IMAGE" bcftools index -f -t "$FILTERED_C"
N_IN=$(run_in "$BCFTOOLS_IMAGE" bcftools view -H "$ANNOTATED_C" | wc -l | tr -d ' ')
N_OUT=$(run_in "$BCFTOOLS_IMAGE" bcftools view -H "$FILTERED_C" | wc -l | tr -d ' ')

echo "=== duphold complete ==="
echo "Annotated: ${OUTPUT_DIR}/${SAMPLE}_sv_duphold.vcf"
echo "Filtered:  ${FILTERED} (kept ${N_OUT} of ${N_IN} records; step 05 annotates this file)"
