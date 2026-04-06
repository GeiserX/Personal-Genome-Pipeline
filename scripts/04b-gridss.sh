#!/usr/bin/env bash
# GRIDSS — Assembly-based structural variant caller
# Input: sorted BAM + GRCh38 reference (with BWA index)
# Output: VCF with BND-notation breakpoints in $GENOME_DIR/<sample>/sv_gridss/
#
# GRIDSS excels at complex rearrangements that Manta/Delly miss.
# Output is BND notation — use the SV consensus merge (step 22) for integration.
#
# HEAVY: Requires 31 GB JVM heap, 8 threads, and ~50 GB intermediate disk space.
# Expected runtime: 4-8 hours for 30X WGS.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-8}
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/sv_gridss"

echo "=== GRIDSS: ${SAMPLE} ==="
echo "BAM: ${BAM}"
echo "WARNING: GRIDSS requires ~31 GB memory and 4-8 hours for 30X WGS."

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# GRIDSS requires classic BWA index files (.amb .ann .bwt .pac .sa) alongside the reference.
# NOTE: BWA-MEM2 index files (.bwt.2bit.64 etc.) are NOT compatible — GRIDSS bundles
# classic bwa internally for its realignment step and needs the classic format.
BWA_MISSING=""
for ext in amb ann bwt pac sa; do
  if [ ! -f "${REF}.${ext}" ]; then
    BWA_MISSING="${BWA_MISSING} .${ext}"
  fi
done
if [ -n "$BWA_MISSING" ]; then
  echo "ERROR: Classic BWA index files missing:${BWA_MISSING}" >&2
  echo "GRIDSS requires classic bwa index files (NOT BWA-MEM2's .bwt.2bit.64)." >&2
  echo "Generate them (~1 hour) with:" >&2
  echo "  docker run --rm -v \"\${GENOME_DIR}:/genome\" quay.io/biocontainers/bwa:0.7.18--he4a0461_1 \\" >&2
  echo "    bwa index /genome/reference/Homo_sapiens_assembly38.fasta" >&2
  exit 1
fi

# Skip if output already exists
if [ -f "${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz" ]; then
  echo "GRIDSS output already exists, skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

mkdir -p "$OUTPUT_DIR"

# Download ENCODE blacklist for hg38 if not present
BLACKLIST="${GENOME_DIR}/reference/ENCFF356LFX.bed"
if [ ! -f "$BLACKLIST" ]; then
  echo "Downloading ENCODE blacklist for GRCh38..."
  wget -q -O "$BLACKLIST" \
    "https://raw.githubusercontent.com/PapenfussLab/gridss/master/example/ENCFF356LFX.bed" || {
    echo "WARNING: Failed to download blacklist. GRIDSS will run without it."
    BLACKLIST=""
  }
fi

# Build GRIDSS command
GRIDSS_ARGS=(
  gridss
  -r /genome/reference/Homo_sapiens_assembly38.fasta
  -o "/genome/${SAMPLE}/sv_gridss/${SAMPLE}_gridss.vcf.gz"
  -a "/genome/${SAMPLE}/sv_gridss/${SAMPLE}_assembly.bam"
  -t "${THREADS}"
  --jvmheap 28g
)

if [ -n "${BLACKLIST}" ] && [ -f "${BLACKLIST}" ]; then
  GRIDSS_ARGS+=(-b /genome/reference/ENCFF356LFX.bed)
fi

GRIDSS_ARGS+=("/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam")

# GRIDSS via Docker Hub image (1.4 GB, includes all dependencies: Java 11, R, bwa, samtools)
echo "Running GRIDSS (this takes 4-8 hours for 30X WGS)..."
docker run --rm --user root \
  --cpus "${THREADS}" --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  -e JAVA_TOOL_OPTIONS="-Xmx28g" \
  quay.io/biocontainers/gridss:2.13.2--h96c455f_6 \
  "${GRIDSS_ARGS[@]}"

echo "=== GRIDSS complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz"
echo ""
echo "NOTE: GRIDSS outputs BND-notation breakpoints. For standard SV types"
echo "  (DEL/DUP/INV/INS), use the SV consensus merge step (22-survivor-merge.sh)"
echo "  which converts and integrates calls from all SV callers."
echo ""
echo "Quality filtering: QUAL >= 1000 with assembly support (AS > 0 & RAS > 0)"
echo "  is a good threshold for high-confidence calls."
