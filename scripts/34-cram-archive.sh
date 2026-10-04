#!/usr/bin/env bash
# Step 34: keep the alignments as CRAM (samtools)
# Input:  ${SAMPLE}/aligned/${SAMPLE}_sorted.bam (+ .bai) and the reference
# Output: ${SAMPLE}/aligned/${SAMPLE}_sorted.cram (+ .crai), about half the
#         size of the BAM, checked against it before anything else happens
#
# Usage: ./scripts/34-cram-archive.sh <sample_name>               write and check the CRAM
#        ./scripts/34-cram-archive.sh <sample_name> --delete-bam  the same, then delete the BAM
#        ./scripts/34-cram-archive.sh <sample_name> --restore     write the BAM back from the CRAM
#
# A CRAM stores each read as its difference from the reference, so it can only
# be read with the same reference FASTA it was written with: keep that file
# (REF_FASTA) as long as you keep the CRAM. The check is that the CRAM passes
# `samtools quickcheck` and that `samtools flagstat` of the CRAM and of the BAM
# agree line for line (the same number of reads, mapped, paired, duplicates).
# The BAM is deleted only with --delete-bam, only after that check passed in
# this run, and never when it failed. --restore writes the BAM back the same
# way, checked against the CRAM, for the bash steps, which read the BAM (the
# Nextflow pipeline reads a CRAM row itself; see docs/34-cram-archive.md).
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [--delete-bam | --restore]}
MODE=${2:-archive}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
case "$MODE" in
  archive|--delete-bam|--restore) ;;
  *) echo "ERROR: unknown option '${MODE}'. Usage: $0 <sample_name> [--delete-bam | --restore]" >&2; exit 2 ;;
esac
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALN="${SAMPLE_DIR}/aligned"
BAM="${ALN}/${SAMPLE}_sorted.bam"
CRAM="${ALN}/${SAMPLE}_sorted.cram"
C="/genome/${SAMPLE}/aligned/${SAMPLE}_sorted"   # the same files inside a container

echo "=== Step 34: CRAM archive: ${SAMPLE} (${MODE#--}) ==="
echo "Reference: ${REF_FASTA}"

for f in "$REF_FASTA" "${REF_FASTA}.fai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done

# flagstat FILE_IN_CONTAINER OUT: samtools flagstat of a BAM or CRAM, to OUT.
# flagstat takes no --reference; a CRAM gets the reference as an input option.
flagstat() {
  local -a ref=()
  case "$1" in *.cram) ref=(--input-fmt-option "reference=${REF_FASTA_C}") ;; esac
  atomic_out "$2" run_in --cpus "$THREADS" --memory 2g "$SAMTOOLS_IMAGE" \
    samtools flagstat -@ "$THREADS" ${ref[@]+"${ref[@]}"} "$1"
}

# same_reads A B: the two flagstat files agree line for line. Prints the reads.
same_reads() {
  if [ ! -s "$1" ] || [ ! -s "$2" ]; then
    echo "  flagstat missing: ${1} or ${2}" >&2
    return 1
  fi
  echo "  reads (BAM):  $(awk 'NR == 1 { print $1 + $3 }' "$1")"
  echo "  reads (CRAM): $(awk 'NR == 1 { print $1 + $3 }' "$2")"
  if ! diff "$1" "$2" >&2; then
    echo "  flagstat differs between the two files (lines above: < BAM, > CRAM)" >&2
    return 1
  fi
}

if [ "$MODE" = --restore ]; then
  for f in "$CRAM" "${CRAM}.crai"; do
    [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
  done
  if [ -e "$BAM" ]; then
    echo "ERROR: ${BAM} exists already; there is nothing to restore. Remove it first to write it again." >&2
    exit 1
  fi
  echo "Writing ${BAM} from ${CRAM}..."
  rm -f "${ALN}/${SAMPLE}_sorted.part.bam" "${ALN}/${SAMPLE}_sorted.part.bam.bai"
  run_in --cpus "$THREADS" --memory 4g "$SAMTOOLS_IMAGE" \
    samtools view -@ "$THREADS" -b --reference "$REF_FASTA_C" \
      -o "${C}.part.bam" "${C}.cram"
  run_in --cpus "$THREADS" --memory 2g "$SAMTOOLS_IMAGE" \
    samtools index -@ "$THREADS" "${C}.part.bam"
  flagstat "${C}.cram" "${ALN}/${SAMPLE}_sorted.cram.flagstat"
  flagstat "${C}.part.bam" "${ALN}/${SAMPLE}_sorted.part.bam.flagstat"
  if ! run_in "$SAMTOOLS_IMAGE" samtools quickcheck -v "${C}.part.bam" \
     || ! same_reads "${ALN}/${SAMPLE}_sorted.part.bam.flagstat" "${ALN}/${SAMPLE}_sorted.cram.flagstat"; then
    rm -f "${ALN}/${SAMPLE}_sorted.part.bam" "${ALN}/${SAMPLE}_sorted.part.bam.bai" "${ALN}/${SAMPLE}_sorted.part.bam.flagstat"
    echo "ERROR: the BAM written from the CRAM does not match it; nothing was kept." >&2
    exit 1
  fi
  mv -f "${ALN}/${SAMPLE}_sorted.part.bam.bai" "${BAM}.bai"
  mv -f "${ALN}/${SAMPLE}_sorted.part.bam" "$BAM"
  rm -f "${ALN}/${SAMPLE}_sorted.part.bam.flagstat"
  # The index is newer than the BAM once both are in place.
  touch "${BAM}.bai"
  echo "=== Step 34 complete: ${BAM} restored and checked against the CRAM ==="
  exit 0
fi

for f in "$BAM" "${BAM}.bai"; do
  [ -f "$f" ] || { echo "ERROR: File not found: ${f}" >&2; exit 1; }
done

# The CRAM is written under a .part name and takes its real name only after
# the check, so a CRAM that exists was checked. An older CRAM is replaced, not
# trusted: the BAM may have been aligned again since.
PART="${ALN}/${SAMPLE}_sorted.part.cram"
discard() {
  rm -f "$PART" "${PART}.crai" "${ALN}/${SAMPLE}_sorted.part.cram.flagstat"
  echo "ERROR: $1 The CRAM was removed; the BAM is kept." >&2
  exit 1
}
rm -f "$PART" "${PART}.crai" "${ALN}/${SAMPLE}_sorted.part.cram.flagstat"
echo "Writing ${CRAM}..."
run_in --cpus "$THREADS" --memory 4g "$SAMTOOLS_IMAGE" \
  samtools view -@ "$THREADS" -C --reference "$REF_FASTA_C" \
    -o "${C}.part.cram" "${C}.bam" || discard "samtools could not write the CRAM."

echo "Checking the CRAM against the BAM..."
run_in "$SAMTOOLS_IMAGE" samtools quickcheck -v "${C}.part.cram" \
  || discard "The CRAM fails samtools quickcheck (truncated, or not a CRAM)."
run_in --cpus "$THREADS" --memory 2g "$SAMTOOLS_IMAGE" \
  samtools index -@ "$THREADS" "${C}.part.cram" || discard "samtools could not index the CRAM."
flagstat "${C}.bam" "${ALN}/${SAMPLE}_sorted.bam.flagstat" || discard "samtools flagstat failed on the BAM."
flagstat "${C}.part.cram" "${ALN}/${SAMPLE}_sorted.part.cram.flagstat" || discard "samtools flagstat failed on the CRAM."
same_reads "${ALN}/${SAMPLE}_sorted.bam.flagstat" "${ALN}/${SAMPLE}_sorted.part.cram.flagstat" \
  || discard "The CRAM does not hold the same reads as the BAM."
mv -f "${PART}.crai" "${CRAM}.crai"
mv -f "$PART" "$CRAM"
mv -f "${ALN}/${SAMPLE}_sorted.part.cram.flagstat" "${ALN}/${SAMPLE}_sorted.cram.flagstat"
touch "${CRAM}.crai"

BAM_KB=$(du -k "$BAM" | awk '{ print $1 }')
CRAM_KB=$(du -k "$CRAM" | awk '{ print $1 }')
echo "  BAM:  $((BAM_KB / 1024)) MB"
echo "  CRAM: $((CRAM_KB / 1024)) MB ($(awk -v c="$CRAM_KB" -v b="$BAM_KB" 'BEGIN { printf "%d", b ? 100 * c / b : 0 }')% of the BAM)"
echo "The CRAM can only be read with ${REF_FASTA}: keep that file as long as you keep the CRAM."

if [ "$MODE" = --delete-bam ]; then
  rm -f "$BAM" "${BAM}.bai"
  echo "Deleted ${BAM} and its index (the CRAM was checked against it in this run)."
  echo "The bash steps read the BAM: write it back with $0 ${SAMPLE} --restore before running one of them."
else
  echo "The BAM is kept. Delete it with: $0 ${SAMPLE} --delete-bam"
fi
echo "=== Step 34 complete: ${CRAM} ==="
