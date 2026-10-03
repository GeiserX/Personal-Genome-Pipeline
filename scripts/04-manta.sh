#!/usr/bin/env bash
# Manta — Structural variant calling (deletions, duplications, inversions, translocations)
# Input: sorted BAM + GRCh38 reference
# Output: diploidSV.vcf.gz (~7-9K structural variants per 30X WGS), with
#         inversions as SVTYPE=INV records; diploidSV.raw.vcf.gz is Manta's
#         own file, where an inversion is a pair of breakend (BND) records.
# Env: THREADS (default 8) caps the container and Manta's -j; ALIGN_DIR
#      (default aligned) is the folder under the sample that holds the BAM;
#      MANTA_CALL_REGIONS, a bgzipped BED under GENOME_DIR with its .tbi
#      beside it, limits calling to those regions (configManta.py
#      --callRegions), e.g. chr1-22, X and Y without the ALT and decoy contigs.
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
MANTA_DIR="${SAMPLE_DIR}/manta"
VARIANTS="${MANTA_DIR}/results/variants"
DIPLOID="${VARIANTS}/diploidSV.vcf.gz"
RAW="${VARIANTS}/diploidSV.raw.vcf.gz"
VARIANTS_C="/genome/${SAMPLE}/manta/results/variants"

echo "=== Manta SV Calling: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Output: ${DIPLOID}"

REGION_ARGS=()
if [ -n "${MANTA_CALL_REGIONS:-}" ]; then
  REGION_ARGS=(--callRegions "$(cpath "$MANTA_CALL_REGIONS")")
  echo "Call regions: ${MANTA_CALL_REGIONS}"
fi

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai" ${MANTA_CALL_REGIONS:+"$MANTA_CALL_REGIONS" "${MANTA_CALL_REGIONS}.tbi"}; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# A finished run is not repeated: configManta.py refuses a runDir that
# already holds a workflow, which made every second run-all report Manta
# as failed.
if [ -f "$DIPLOID" ] && [ -f "${DIPLOID}.tbi" ] && [ -f "$RAW" ] && [ -f "${RAW}.tbi" ]; then
  echo "Manta already done for ${SAMPLE}: ${DIPLOID}"
  echo "To run it again, delete ${MANTA_DIR}/ first."
  exit 0
fi

if [ -f "$RAW" ] || { [ -f "$DIPLOID" ] && [ -f "${DIPLOID}.tbi" ]; }; then
  # Manta finished earlier (a run from before the inversion conversion, or one
  # stopped during it): only the conversion below is left.
  echo "Manta results found in ${VARIANTS}/; converting the inversions only."
else
  # Step 1: Configure Manta (skipped when an interrupted run can be resumed)
  if [ -f "${MANTA_DIR}/runWorkflow.py" ]; then
    echo "Found an unfinished Manta run in ${MANTA_DIR}/; resuming it."
  else
    if [ -d "$MANTA_DIR" ]; then
      # Files written by the container belong to root, so remove them from a container
      echo "Removing leftover ${MANTA_DIR}/ (no workflow and no results)..."
      # --root: a directory left by an earlier version of this script belongs to root.
      run_in --root \
        "${MANTA_IMAGE}" \
        rm -rf "/genome/${SAMPLE}/manta"
    fi
    echo "Configuring Manta..."
    run_in \
      --cpus "$THREADS" --memory 16g \
      "${MANTA_IMAGE}" \
      configManta.py \
        --bam "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
        --referenceFasta "${REF_FASTA_C}" \
        ${REGION_ARGS[@]+"${REGION_ARGS[@]}"} \
        --runDir "/genome/${SAMPLE}/manta"
  fi

  # Step 2: Run Manta workflow
  echo "Running Manta (this takes 1-3 hours for 30X WGS)..."
  run_in \
    --cpus "$THREADS" --memory 16g \
    "${MANTA_IMAGE}" \
    "/genome/${SAMPLE}/manta/runWorkflow.py" -j "$THREADS"

  if [ ! -f "$DIPLOID" ] || [ ! -f "${DIPLOID}.tbi" ]; then
    echo "ERROR: Manta finished without ${DIPLOID}(.tbi)." >&2
    echo "  Delete ${MANTA_DIR}/ and run this step again." >&2
    exit 1
  fi
fi

# Step 3: Inversions. Manta writes an inversion as two breakend (BND) records;
# AnnotSV, duphold and the SV consensus (steps 05, 15, 22) read SVTYPE, so
# libexec/convertInversion.py, which ships with Manta, turns each pair into
# one SVTYPE=INV record. It needs a samtools and the reference: both come from
# the Manta image. Manta's own file stays as diploidSV.raw.vcf.gz, and the
# converted file takes the name the later steps read.
# Each file moves on its own, so a run stopped between the two moves finishes
# them on the next run.
if [ ! -f "$RAW" ]; then
  mv -f "$DIPLOID" "$RAW"
fi
if [ ! -f "${RAW}.tbi" ] && [ ! -f "$DIPLOID" ] && [ -f "${DIPLOID}.tbi" ]; then
  mv -f "${DIPLOID}.tbi" "${RAW}.tbi"
fi
if [ ! -f "${RAW}.tbi" ]; then
  echo "ERROR: ${RAW} has no .tbi beside it." >&2
  echo "  Delete ${MANTA_DIR}/ and run this step again." >&2
  exit 1
fi
# The breakends convertInversion.py converts: ALT [p[t or t]p] with the mate
# on the same chromosome, two records per inversion.
N_BND=$(gzip -cd "$RAW" | awk -F'\t' '!/^#/ && ($5 ~ /^\[/ || $5 ~ /\]$/) {
  m = $5; sub(/^[^][]*[][]/, "", m); sub(/:.*/, "", m); if (m == $1) n++ } END { print n + 0 }')
if [ "$N_BND" -eq 0 ]; then
  # Nothing to convert, so no container: the calls stay as Manta wrote them.
  cp -f "$RAW" "$DIPLOID"
  cp -f "${RAW}.tbi" "${DIPLOID}.tbi"
  echo "Inversion conversion: 0 inversion breakend records in Manta's calls, nothing to convert; ${DIPLOID} is a copy of ${RAW}"
else
  echo "Converting inversions (libexec/convertInversion.py)..."
  # shellcheck disable=SC2016  # the $ expressions expand inside the container
  atomic_out "$DIPLOID" run_in "${MANTA_IMAGE}" bash -c '
    set -euo pipefail
    libexec="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"
    "${libexec}/convertInversion.py" "${libexec}/samtools" "$1" "$2" | "${libexec}/bgzip" -c
  ' _ "$REF_FASTA_C" "${VARIANTS_C}/diploidSV.raw.vcf.gz"
  # shellcheck disable=SC2016  # expands inside the container
  run_in "${MANTA_IMAGE}" bash -c \
    '"$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec/tabix" -f -p vcf "$1"' \
    _ "${VARIANTS_C}/diploidSV.vcf.gz"
  if ! have_output "$DIPLOID" || [ ! -s "${DIPLOID}.tbi" ]; then
    echo "ERROR: the inversion conversion did not write ${DIPLOID}(.tbi)." >&2
    exit 1
  fi
  N_INV=$(gzip -cd "$DIPLOID" | awk -F'\t' '!/^#/ && $8 ~ /(^|;)SVTYPE=INV(;|$)/ { n++ } END { print n + 0 }')
  echo "Inversion conversion: ${N_BND} inversion breakend records in, ${N_INV} SVTYPE=INV records out"
fi

echo "=== Manta complete ==="
echo "Diploid SVs: ${DIPLOID} (Manta's own file: ${RAW})"
echo "Candidates: ${VARIANTS}/candidateSV.vcf.gz"
echo ""
echo "SV count: $(gzip -cd "$DIPLOID" | grep -vc '^#' || true)"
echo ""
echo "Next: run 05-annotsv.sh to classify pathogenicity (ACMG)"
