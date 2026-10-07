#!/usr/bin/env bash
# TIDDIT — Alternative SV caller (large structural variants)
# Alternative to step 04 (Manta). Outputs to sv_tiddit/ to avoid conflicts.
# Input: sorted BAM + GRCh38 reference
# Output: SV VCF in $GENOME_DIR/<sample>/sv_tiddit/
# Runtime: ~30-60 minutes per 30X genome (with --skip_assembly)
# NOTE: TIDDIT's local assembly realigns contigs with classic `bwa mem`, so it
# needs the classic BWA index next to the reference (.amb .ann .bwt .pac .sa,
# the files GRIDSS needs too). BWA-MEM2's index (.bwt.2bit.64) does not count.
# Without the classic index the step runs with --skip_assembly.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/sv_tiddit"

echo "=== TIDDIT SV Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"


# Local assembly only with the classic BWA index: with BWA-MEM2's index alone
# TIDDIT's bwa call fails and the run stops with "file does not contain
# alignment data".
TIDDIT_EXTRA_ARGS=()
BWA_MISSING=""
for ext in amb ann bwt pac sa; do
  [ -f "${REF}.${ext}" ] || BWA_MISSING="${BWA_MISSING} .${ext}"
done
if [ -z "$BWA_MISSING" ]; then
  echo "Classic BWA index found: local assembly enabled for breakpoint refinement."
else
  echo "No classic BWA index (missing:${BWA_MISSING}): assembly skipped (--skip_assembly)."
  if [ -f "${REF}.bwt.2bit.64" ]; then
    echo "  The BWA-MEM2 index next to the reference is not one TIDDIT can use."
  fi
  TIDDIT_EXTRA_ARGS+=(--skip_assembly)
fi

# TIDDIT exits 0 when one of its own checks fails, so its output is checked
# below and its log kept for the error message.
TIDDIT_LOG="${OUTPUT_DIR}/${SAMPLE}_tiddit.log"
rm -f "${OUTPUT_DIR}/${SAMPLE}.vcf"
echo "[1/3] Running TIDDIT SV caller..."
run_in --cpus "$THREADS" --memory 8g \
  "$TIDDIT_IMAGE" \
  tiddit --sv \
    --bam "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    --ref "${REF_FASTA_C}" \
    --threads "$THREADS" \
    ${TIDDIT_EXTRA_ARGS[@]+"${TIDDIT_EXTRA_ARGS[@]}"} \
    -o "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}" 2>&1 | tee "$TIDDIT_LOG"

if ! have_output "${OUTPUT_DIR}/${SAMPLE}.vcf"; then
  echo "ERROR: TIDDIT wrote no VCF (${OUTPUT_DIR}/${SAMPLE}.vcf). The end of its log:" >&2
  tail -n 20 "$TIDDIT_LOG" >&2
  exit 1
fi

echo "[2/3] Compressing VCF with bcftools..."
run_in "$BCFTOOLS_IMAGE" \
  bcftools view \
    "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}.vcf" \
    -Oz -o "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz"

echo "[3/3] Indexing VCF..."
run_in "$BCFTOOLS_IMAGE" \
  bcftools index -f -t \
    "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz"

SV_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools stats "/genome/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz" \
  | grep '^SN' | grep 'number of records' | awk '{print $NF}')
SV_COUNT=${SV_COUNT:-unknown}

echo "=== TIDDIT complete ==="
echo "Total SVs called: ${SV_COUNT}"
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_sv.vcf.gz"
echo "Auxiliary files: ${OUTPUT_DIR}/${SAMPLE}.ploidies.tab, ${OUTPUT_DIR}/${SAMPLE}.signals.tab"
echo ""
echo "Count by SV type:"
printf '  bcftools query -f "%%INFO/SVTYPE\\n" %s/%s_sv.vcf.gz | sort | uniq -c\n' "${OUTPUT_DIR}" "${SAMPLE}"
