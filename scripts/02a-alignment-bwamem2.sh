#!/usr/bin/env bash
# BWA-MEM2 — Alternative aligner (FASTQ to sorted BAM)
# Alternative to step 02 (minimap2). Outputs to aligned_bwamem2/ to avoid conflicts.
# Input: paired-end FASTQ files + GRCh38 reference
# Output: sorted BAM + BAI index in $GENOME_DIR/<sample>/aligned_bwamem2/
# Note: BWA-MEM2 produces XS (suboptimal alignment score) tags that some callers
#       (especially Strelka2) depend on. minimap2 does not produce XS tags.
# Memory: the one-time index build needs about 28 GB of RAM per Gbp of
#   reference, so about 90 GB for GRCh38 (upstream's figure). The container
#   gets no memory cap; a host with less RAM kills the build (exit 137).
#   Alignment itself needs about 16 GB.
# Reads go bwa-mem2 -> samtools fixmate -m -> sort -> markdup in one pipe (no
# SAM on disk), and the BAM is renamed into place only once it is indexed and
# passes samtools quickcheck, as in step 02. THREADS (default 8) sets the CPUs.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"

# Detect trimmed FASTQs (same logic as 02-alignment.sh)
if [ -n "${FASTQ_SUBDIR:-}" ]; then
  echo "Using explicit FASTQ_SUBDIR=${FASTQ_SUBDIR}."
elif [ -f "${SAMPLE_DIR}/fastq_trimmed/${SAMPLE}_R1.fastq.gz" ] && \
     [ -f "${SAMPLE_DIR}/fastq_trimmed/${SAMPLE}_R2.fastq.gz" ]; then
  FASTQ_SUBDIR="fastq_trimmed"
else
  FASTQ_SUBDIR="fastq"
fi

R1="${SAMPLE_DIR}/${FASTQ_SUBDIR}/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/${FASTQ_SUBDIR}/${SAMPLE}_R2.fastq.gz"
REF="$REF_FASTA"
BWA_INDEX="${REF}.bwt.2bit.64"
OUTPUT_DIR="${SAMPLE_DIR}/aligned_bwamem2"

echo "=== BWA-MEM2 Alignment: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Reference: ${REF}"
echo "Output: ${OUTPUT_DIR}/"

# Validate inputs
for f in "$R1" "$R2" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

BAM="${OUTPUT_DIR}/${SAMPLE}_sorted.bam"
TMP_BAM="${OUTPUT_DIR}/${SAMPLE}_sorted.tmp.bam"
SORT_TMP="${OUTPUT_DIR}/${SAMPLE}.sort_tmp"
# bwa-mem2 index -p writes these five files; .bwt.2bit.64 is the one this
# script tests for, so it is renamed last.
IDX_TMP="${REF}.tmp.$$"
IDX_EXTS=(0123 amb ann pac bwt.2bit.64)
cleanup() {
  local e
  rm -rf "$TMP_BAM" "${TMP_BAM}.bai" "$SORT_TMP"
  for e in "${IDX_EXTS[@]}"; do rm -f "${IDX_TMP}.${e}"; done
}
trap cleanup EXIT

# Step 1: Build the BWA-MEM2 index if not present (one-time, ~1 hour), under a
# temporary prefix: an interrupted build leaves no .bwt.2bit.64 to trust.
if [ ! -f "$BWA_INDEX" ]; then
  echo "=== Building BWA-MEM2 index (one-time, ~1 hour, about 90 GB of RAM for GRCh38) ==="
  # The index files go next to the FASTA, so reference/ is writable here.
  rc=0
  run_in --rw "$(dirname "$REF_FASTA")" \
    --cpus "${THREADS}" \
    "${BWAMEM2_IMAGE}" \
    bwa-mem2 index -p "$(cpath "$IDX_TMP")" "${REF_FASTA_C}" || rc=$?
  if [ "$rc" -eq 137 ]; then
    echo "ERROR: bwa-mem2 index was killed (exit 137), almost always for lack of memory." >&2
    echo "  It needs about 28 GB of RAM per Gbp of reference: about 90 GB for GRCh38." >&2
    echo "  Build the index on a larger machine and copy the five ${REF##*/}.* index files next to the FASTA," >&2
    echo "  or use step 02 (minimap2), which needs about 20 GB of RAM for GRCh38." >&2
    exit 1
  elif [ "$rc" -ne 0 ]; then
    echo "ERROR: bwa-mem2 index failed (exit ${rc})." >&2
    exit "$rc"
  fi
  for e in "${IDX_EXTS[@]}"; do mv -f "${IDX_TMP}.${e}" "${REF}.${e}"; done
  echo "BWA-MEM2 index built."
else
  echo "BWA-MEM2 index found, skipping build."
fi

# Step 2: Align with BWA-MEM2, mark duplicates and sort (4-8 hours for 30X WGS).
# bwa-mem2 keeps the two reads of a pair together, which fixmate needs.
echo "=== Aligning reads with BWA-MEM2 and marking duplicates (this takes 4-8 hours for 30X WGS) ==="
SORT_MEM_GB=$((THREADS + 4))
mkdir -p "$SORT_TMP"
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
run_in \
  --cpus "${THREADS}" --memory 24g \
  "${BWAMEM2_IMAGE}" \
  bwa-mem2 mem -t "${THREADS}" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
    "${REF_FASTA_C}" \
    "/genome/${SAMPLE}/${FASTQ_SUBDIR}/${SAMPLE}_R1.fastq.gz" \
    "/genome/${SAMPLE}/${FASTQ_SUBDIR}/${SAMPLE}_R2.fastq.gz" \
| run_in -i \
  --cpus "${THREADS}" --memory "${SORT_MEM_GB}g" \
  "${SAMTOOLS_IMAGE}" \
  bash -euo pipefail -c '
    threads=$1 tmp=$2 out=$3
    samtools fixmate -u -m - - \
      | samtools sort -u -@ "$threads" -m 1G -T "${tmp}/sort" - \
      | samtools markdup -@ "$threads" -T "${tmp}/markdup" - "$out"' \
  _ "${THREADS}" "$(cpath "$SORT_TMP")" "$(cpath "$TMP_BAM")"

# Step 3: Index and check, then rename (the old index goes first).
echo "=== Indexing and checking BAM ==="
run_in \
  --cpus "${THREADS}" --memory 2g \
  "${SAMTOOLS_IMAGE}" \
  samtools index -@ "${THREADS}" "$(cpath "$TMP_BAM")"
run_in \
  --cpus 1 --memory 1g \
  "${SAMTOOLS_IMAGE}" \
  samtools quickcheck -v "$(cpath "$TMP_BAM")"
rm -f "${BAM}.bai"
mv -f "$TMP_BAM" "$BAM"
mv -f "${TMP_BAM}.bai" "${BAM}.bai"

echo "=== BWA-MEM2 Alignment complete ==="
echo "BAM: ${BAM}"
echo "Index: ${BAM}.bai"
ls -lh "$BAM" 2>/dev/null || true
echo ""
echo "Next step: call variants from this BAM into its own directory, so the"
echo "minimap2 VCF in vcf/ is kept:"
echo "  ALIGN_DIR=aligned_bwamem2 VCF_OUT_DIR=vcf_bwamem2 ./scripts/03-deepvariant.sh ${SAMPLE} [male|female]"
