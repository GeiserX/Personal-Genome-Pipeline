#!/usr/bin/env bash
# HLA Typing — T1K (Class I + II, 4-digit resolution)
# Types HLA-A, B, C (Class I) and DRB1, DQB1, DPB1 (Class II)
# Uses IPD-IMGT/HLA database aligned against GRCh38 reference
#
# The database is a named IPD-IMGT/HLA release (HLA_DB_RELEASE, set in
# scripts/lib/common.sh), not whatever is current on the day of the run. Its
# T1K index lives in a directory named after the T1K version and that release,
# so a change of either builds a new index. The release is written next to the
# genotype file, in database_release.txt.
#
# The coordinate file (T1K's -c) takes each gene's GRCh38 position from a
# gene annotation (GENCODE's basic GTF), as T1K's README says. A FASTA or its
# .fai there gives every gene "-1 -1" coordinates, and then T1K extracts no
# reads from the BAM; the step stops when a typed gene has no coordinates.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${GENOME_DIR}/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/hla_t1k"
# The genes the summary below reports; each must have coordinates.
TYPED_GENES="HLA-A HLA-B HLA-C HLA-DRB1 HLA-DQB1 HLA-DPB1"

echo "=== T1K HLA Typing: ${SAMPLE} ==="

for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"

# The HLA database and the gene coordinates are installed by setup.sh. Without
# them there is nothing to type against: the step says so and stops without
# failing, as other steps do when their data is not installed.
HLA_DAT=$(data_file hla_dat) || {
  echo "SKIPPED: IPD-IMGT/HLA ${HLA_DB_RELEASE} is not installed (${HLA_DAT})."
  echo "  Install it with: ./scripts/setup.sh ${GENOME_DIR}"
  exit 0
}
GENES_GTF=$(data_file gencode_genes) || {
  echo "SKIPPED: the GENCODE ${GENCODE_RELEASE} gene coordinates are not installed (${GENES_GTF})."
  echo "  Install them with: ./scripts/setup.sh ${GENOME_DIR}"
  exit 0
}

T1K_VERSION=${T1K_IMAGE##*:}
T1K_VERSION=${T1K_VERSION%%--*}
IDX_ROOT="${GENOME_DIR}/t1k_idx"
IDX_DIR="${IDX_ROOT}/t1k-${T1K_VERSION}_imgt-${HLA_DB_RELEASE}_gencode-${GENCODE_RELEASE}"

# Step 1: the HLA index with coordinates (one-time per T1K version and database
# release, ~5 min). The index is shared by every sample: one run builds it
# while the others wait, into a .part directory renamed when complete.
mkdir -p "$IDX_ROOT"
LOCK="${IDX_ROOT}/.build.lock"
lock_acquire "$LOCK"
trap 'lock_release "$LOCK"' EXIT
if [ ! -d "$IDX_DIR" ]; then
  echo "Building the T1K ${T1K_VERSION} index for IPD-IMGT/HLA ${HLA_DB_RELEASE}..."
  rm -rf "${IDX_DIR}.part"
  run_in --rw "$IDX_ROOT" --cpus 2 --memory 4g \
    "${T1K_IMAGE}" \
    t1k-build.pl \
      -d "$(cpath "$HLA_DAT")" \
      -g "$(cpath "$GENES_GTF")" \
      --prefix hla \
      -o "$(cpath "${IDX_DIR}.part")"
  COORD=$(find "${IDX_DIR}.part" -maxdepth 1 -name '*dna_coord.fa' | head -n 1)
  if [ -z "$COORD" ] || [ -z "$(find "${IDX_DIR}.part" -maxdepth 1 -name '*dna_seq.fa' | head -n 1)" ]; then
    echo "ERROR: t1k-build.pl wrote no *dna_seq.fa and *dna_coord.fa in ${IDX_DIR}.part" >&2
    exit 1
  fi
  NO_COORD=""
  for g in $TYPED_GENES; do
    if grep -q "^>${g}\*[^ ]* [^ ]* -1 -1 " "$COORD" || ! grep -q "^>${g}\*" "$COORD"; then
      NO_COORD="${NO_COORD} ${g}"
    fi
  done
  if [ -n "$NO_COORD" ]; then
    echo "ERROR: genes without GRCh38 coordinates in ${COORD}:${NO_COORD}" >&2
    echo "  T1K would extract no reads for them. Check ${GENES_GTF}." >&2
    exit 1
  fi
  echo "  $(grep -c ' -1 -1 ' "$COORD" || true) alleles of genes outside the GRCh38 primary assembly have no coordinates (not typed)."
  mv "${IDX_DIR}.part" "$IDX_DIR"
fi
lock_release "$LOCK"
trap - EXIT

# File names as the Nextflow module finds them: T1K's prefix has changed
# between versions.
SEQ_FA=$(find "$IDX_DIR" -maxdepth 1 -name '*dna_seq.fa' | head -n 1)
COORD_FA=$(find "$IDX_DIR" -maxdepth 1 -name '*dna_coord.fa' | head -n 1)
if [ -z "$SEQ_FA" ] || [ -z "$COORD_FA" ]; then
  echo "ERROR: ${IDX_DIR} has no *dna_seq.fa or *dna_coord.fa; remove it to rebuild." >&2
  exit 1
fi

# Step 2: Run HLA typing
echo "Running T1K genotyping..."
run_in \
  --cpus "${THREADS}" --memory 8g \
  "${T1K_IMAGE}" \
  run-t1k \
    -b "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
    -f "$(cpath "$SEQ_FA")" \
    -c "$(cpath "$COORD_FA")" \
    --preset hla-wgs \
    -t "${THREADS}" \
    --od "/genome/${SAMPLE}/hla_t1k/" \
    -o "${SAMPLE}_hla"

# The database the calls come from, read from hla.dat itself, beside the calls.
RELEASE_LINE=$(grep -m 1 'IPD-IMGT/HLA Release' "$HLA_DAT" | sed 's/^CC *//') || RELEASE_LINE=""
{
  echo "database: ${RELEASE_LINE:-IPD-IMGT/HLA ${HLA_DB_RELEASE} (no release line in hla.dat)}"
  echo "t1k: ${T1K_VERSION}"
  echo "coordinates: GENCODE ${GENCODE_RELEASE} basic annotation"
} > "${OUTPUT_DIR}/database_release.txt"

echo "=== T1K complete ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_hla_genotype.tsv"
echo "Database: ${OUTPUT_DIR}/database_release.txt ($(head -n 1 "${OUTPUT_DIR}/database_release.txt"))"
