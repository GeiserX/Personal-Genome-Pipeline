#!/usr/bin/env bash
# Alignment — minimap2 + samtools (FASTQ to sorted, duplicate-marked BAM)
# Input: paired-end FASTQ files + GRCh38 reference
# Output: sorted BAM + BAI index in $GENOME_DIR/<sample>/aligned/
#
# Reads go minimap2 -> samtools fixmate -m -> sort -> markdup, so PCR and
# optical duplicates carry the 0x400 flag: GATK, FreeBayes, Octopus and the SV
# and depth steps would otherwise count them as independent reads.
# The BAM is written under a temporary name, indexed and checked with
# samtools quickcheck, and only then renamed: a killed run leaves no
# <sample>_sorted.bam behind for the next run to trust.
# THREADS (default 8) sets the CPUs of the aligner and of samtools.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"

# Allow explicit override (e.g., FASTQ_SUBDIR=fastq to use raw reads even when trimmed exist)
if [ -n "${FASTQ_SUBDIR:-}" ]; then
  echo "Using explicit FASTQ_SUBDIR=${FASTQ_SUBDIR}."
elif [ -f "${SAMPLE_DIR}/fastq_trimmed/${SAMPLE}_R1.fastq.gz" ] && [ -f "${SAMPLE_DIR}/fastq_trimmed/${SAMPLE}_R2.fastq.gz" ]; then
  FASTQ_SUBDIR="fastq_trimmed"
  echo "Using trimmed FASTQs from fastp."
else
  FASTQ_SUBDIR="fastq"
fi
R1="${SAMPLE_DIR}/${FASTQ_SUBDIR}/${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE_DIR}/${FASTQ_SUBDIR}/${SAMPLE}_R2.fastq.gz"
REF="$REF_FASTA"
# The index is built with the same preset the reads are mapped with (-x sr:
# k21, w11; a plain `minimap2 -d` builds k15, w10) and is named after the
# reference, so another REF_FASTA gets its own index.
REF_BASE="${REF_FASTA%.gz}"
REF_BASE="${REF_BASE%.*}"
MMI="${REF_BASE}.sr.mmi"
OUTPUT_DIR="${SAMPLE_DIR}/aligned"

echo "=== Alignment: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Reference: ${REF}"

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
MMI_TMP="${MMI}.tmp.$$"
# Whatever this run leaves half written goes when it exits, finished or not.
cleanup() { rm -rf "$TMP_BAM" "${TMP_BAM}.bai" "$SORT_TMP" "$MMI_TMP"; }
trap cleanup EXIT

# Step 1: Build the minimap2 index (one-time, ~30 min), under a temporary
# name first: an interrupted build leaves no file at the name step 2 trusts.
if [ ! -s "$MMI" ]; then
  if [ -f "${GENOME_DIR}/reference/GRCh38.mmi" ]; then
    echo "NOTE: reference/GRCh38.mmi was built without the sr preset and is no longer used; you can delete it."
  fi
  echo "Building minimap2 index (one-time, ~30 min)..."
  # The index is shared by every sample, so its directory is writable here.
  run_in --rw "$(dirname "$MMI")" \
    --cpus "${THREADS}" --memory 16g \
    "${MINIMAP2_IMAGE}" \
    minimap2 -x sr -t "${THREADS}" -d "$(cpath "$MMI_TMP")" \
      "${REF_FASTA_C}"
  mv -f "$MMI_TMP" "$MMI"
fi

# Step 2: Align, mark duplicates and sort (1-2 hours for 30X WGS)
# minimap2 runs in its own container and pipes SAM to samtools; the -i flag on
# the samtools container keeps stdin open for the pipe. minimap2 writes the
# two reads of a pair next to each other, which fixmate needs; markdup needs
# the mate tags fixmate -m adds and a coordinate-sorted input.
# -R writes a read group: GATK steps (20, 03a, 29) reject reads without one,
# and callers take the sample name from its SM field.
# samtools sort spills to SORT_TMP in the sample directory, not to the
# container's own disk; -m is per thread, so the container gets THREADS + 4 GB.
SORT_MEM_GB=$((THREADS + 4))
echo "Aligning reads and marking duplicates (this takes 1-2 hours for 30X WGS)..."
mkdir -p "$SORT_TMP"
# shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
run_in \
  --cpus "${THREADS}" --memory 16g \
  "${MINIMAP2_IMAGE}" \
  minimap2 -t "${THREADS}" -a -x sr \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
    "$(cpath "$MMI")" \
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

# Step 3: Index and check, then rename. The old index goes first, so an index
# never sits next to a BAM it was not built from.
echo "Indexing and checking BAM..."
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

echo "=== Alignment complete ==="
echo "BAM: ${BAM}"
echo "Index: ${BAM}.bai"
ls -lh "$BAM" 2>/dev/null || true
