#!/usr/bin/env bash
# Manta — Structural variant calling (deletions, duplications, inversions, translocations)
# Input: sorted BAM + GRCh38 reference
# Output: diploidSV.vcf.gz (~7-9K structural variants per 30X WGS)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"
MANTA_DIR="${SAMPLE_DIR}/manta"

echo "=== Manta SV Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${MANTA_DIR}/results/variants/diploidSV.vcf.gz"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# A finished run is not repeated: configManta.py refuses a runDir that
# already holds a workflow, which made every second run-all report Manta
# as failed.
if [ -f "${MANTA_DIR}/results/variants/diploidSV.vcf.gz.tbi" ]; then
  echo "Manta already done for ${SAMPLE}: ${MANTA_DIR}/results/variants/diploidSV.vcf.gz"
  echo "To run it again, delete ${MANTA_DIR}/ first."
  exit 0
fi

# Step 1: Configure Manta (skipped when an interrupted run can be resumed)
if [ -f "${MANTA_DIR}/runWorkflow.py" ]; then
  echo "Found an unfinished Manta run in ${MANTA_DIR}/; resuming it."
else
  if [ -d "$MANTA_DIR" ]; then
    # Files written by the container belong to root, so remove them from a container
    echo "Removing leftover ${MANTA_DIR}/ (no workflow and no results)..."
    docker run --rm \
      -v "${GENOME_DIR}:/genome" \
      quay.io/biocontainers/manta:1.6.0--h9ee0642_2 \
      rm -rf "/genome/${SAMPLE}/manta"
  fi
  echo "Configuring Manta..."
  docker run --rm \
    --cpus 8 --memory 16g \
    -v "${GENOME_DIR}:/genome" \
    quay.io/biocontainers/manta:1.6.0--h9ee0642_2 \
    configManta.py \
      --bam "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" \
      --referenceFasta /genome/reference/Homo_sapiens_assembly38.fasta \
      --runDir "/genome/${SAMPLE}/manta"
fi

# Step 2: Run Manta workflow
echo "Running Manta (this takes 1-3 hours for 30X WGS)..."
docker run --rm \
  --cpus 8 --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  quay.io/biocontainers/manta:1.6.0--h9ee0642_2 \
  "/genome/${SAMPLE}/manta/runWorkflow.py" -j 8

echo "=== Manta complete ==="
echo "Diploid SVs: ${MANTA_DIR}/results/variants/diploidSV.vcf.gz"
echo "Candidates: ${MANTA_DIR}/results/variants/candidateSV.vcf.gz"
echo ""
echo "SV count: $(zgrep -c -v '^#' "${MANTA_DIR}/results/variants/diploidSV.vcf.gz" 2>/dev/null || echo 'check output manually')"
echo ""
echo "Next: run 05-annotsv.sh to classify pathogenicity (ACMG)"
