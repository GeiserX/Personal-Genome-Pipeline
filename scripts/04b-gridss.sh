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
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
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
  echo "  docker run --rm -v \"${GENOME_DIR}:/genome\" ${BWA_IMAGE} \\" >&2
  echo "    bwa index ${REF_FASTA_C}" >&2
  exit 1
fi

# Skip only a finished VCF (complete BGZF file with a VCF header); a file cut
# short by a killed run is called again.
if have_output "${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz"; then
  echo "GRIDSS output already exists, skipping."
  echo "Delete to re-run: rm -rf ${OUTPUT_DIR}"
  exit 0
fi

# GRIDSS fails without a clear message when its JVM cannot get the memory
# (docs/lessons-learned.md). Docker's own limit is what counts, so the step is
# skipped, with the reason, when Docker has less than GRIDSS_MIN_MEM_GB (32).
GRIDSS_MIN_MEM_GB=${GRIDSS_MIN_MEM_GB:-32}
DOCKER_MEM=$("$CONTAINER_ENGINE" info --format '{{.MemTotal}}' 2>/dev/null || true)
if [[ "$DOCKER_MEM" =~ ^[0-9]+$ ]]; then
  if [ $((DOCKER_MEM / 1000000000)) -lt "$GRIDSS_MIN_MEM_GB" ]; then
    echo "SKIPPED: Docker has $((DOCKER_MEM / 1000000000)) GB of memory; GRIDSS needs ${GRIDSS_MIN_MEM_GB} GB."
    echo "  Give Docker more memory, or set GRIDSS_MIN_MEM_GB to try with less."
    exit 0
  fi
else
  echo "WARNING: could not read Docker's memory limit; GRIDSS needs ${GRIDSS_MIN_MEM_GB} GB."
fi

mkdir -p "$OUTPUT_DIR"

# Download the ENCODE blacklist for hg38 if not present, from the commit of
# the GRIDSS repository pinned in versions.env. An empty file left by an
# earlier version of this script is fetched again.
BLACKLIST="${GENOME_DIR}/reference/ENCFF356LFX.bed"
if [ ! -s "$BLACKLIST" ]; then
  echo "Downloading ENCODE blacklist for GRCh38..."
  rm -f "$BLACKLIST"
  fetch "https://raw.githubusercontent.com/PapenfussLab/gridss/${GRIDSS_BLACKLIST_COMMIT}/example/ENCFF356LFX.bed" "$BLACKLIST" || {
    echo "WARNING: Failed to download blacklist. GRIDSS will run without it."
    BLACKLIST=""
  }
fi

# Build GRIDSS command
GRIDSS_ARGS=(
  gridss
  -r "${REF_FASTA_C}"
  -o "/genome/${SAMPLE}/sv_gridss/${SAMPLE}_gridss.vcf.gz"
  -a "/genome/${SAMPLE}/sv_gridss/${SAMPLE}_assembly.bam"
  -t "${THREADS}"
  --jvmheap 28g
  --workingdir "/genome/${SAMPLE}/sv_gridss/work"
)

if [ -n "${BLACKLIST}" ] && [ -f "${BLACKLIST}" ]; then
  GRIDSS_ARGS+=(-b /genome/reference/ENCFF356LFX.bed)
fi

GRIDSS_ARGS+=("/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam")

# GRIDSS via Docker Hub image (1.4 GB, includes all dependencies: Java 11, R, bwa, samtools)
echo "Running GRIDSS (this takes 4-8 hours for 30X WGS)..."
# --rw reference/: on its first run GRIDSS writes <reference>.gridsscache,
# <reference>.img and <reference>.dict next to the FASTA (and a lock directory
# while it does). --workingdir: its ~50 GB of intermediate files go to
# sv_gridss/work/, removed once the VCF is written; -w puts its log in sv_gridss/.
run_in --rw "$(dirname "$REF_FASTA")" -w "/genome/${SAMPLE}/sv_gridss" \
  --cpus "${THREADS}" --memory 32g \
  -e JAVA_TOOL_OPTIONS="-Xmx28g" \
  "${GRIDSS_IMAGE}" \
  "${GRIDSS_ARGS[@]}"

if ! wrote_vcf "${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz"; then
  echo "ERROR: GRIDSS exited without a ${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz." >&2
  echo "  Its intermediate files are kept in ${OUTPUT_DIR}/work/ for a look." >&2
  exit 1
fi
rm -rf "${OUTPUT_DIR}/work"

echo "=== GRIDSS complete ==="
echo "VCF: ${OUTPUT_DIR}/${SAMPLE}_gridss.vcf.gz"
echo ""
echo "NOTE: GRIDSS outputs BND-notation breakpoints. For standard SV types"
echo "  (DEL/DUP/INV/INS), use the SV consensus merge step (22-survivor-merge.sh)"
echo "  which converts and integrates calls from all SV callers."
echo ""
echo "Quality filtering: QUAL >= 1000 with assembly support (AS > 0 & RAS > 0)"
echo "  is a good threshold for high-confidence calls."
