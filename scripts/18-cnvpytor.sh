#!/usr/bin/env bash
# CNVpytor — Depth-based CNV detection (maintained Python successor to CNVnator; orthogonal to Manta/Delly)
# Input:  Sorted BAM (GRCh38) + pinned CNVpytor GC/mask resource files (see docs/00-reference-setup.md)
# Output: CNV calls (TSV) + normalized VCF (deletions / duplications / LOH)
# Runtime: ~1-3 hours per 30X genome
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/aligned/${SAMPLE}_sorted.bam"
REF_FAI="reference/Homo_sapiens_assembly38.fasta.fai"   # relative to the /genome mount
CNVPYTOR_DATA="${GENOME_DIR}/reference/cnvpytor"          # pinned GC/mask .pytor resource files
OUTPUT_DIR="${SAMPLE_DIR}/cnvpytor"
BIN_SIZE=1000

# The pinned biocontainer ships WITHOUT the reference GC/mask resources and its
# built-in `-download` is broken in 1.3.2, so we bind-mount pre-fetched pinned
# files onto the container's package data dir. This path is stable for the
# pinned CNVPYTOR_IMAGE; if the image is bumped, re-verify with:
#   docker run --rm <img> python -c 'import cnvpytor,os;print(os.path.dirname(cnvpytor.__file__)+"/data")'
CNVPYTOR_IMG_DATA="/usr/local/lib/python3.12/site-packages/cnvpytor/data"

# shellcheck source=/dev/null
source "$(dirname "$0")/../versions.env" 2>/dev/null || {
  CNVPYTOR_IMAGE="quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0"
  BCFTOOLS_IMAGE="staphb/bcftools:1.21"
}

echo "=== CNVpytor: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Bin size:  ${BIN_SIZE} bp"
echo "Output:    ${OUTPUT_DIR}"

# Validate inputs
for f in "$BAM" "${BAM}.bai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done
# hg38 GC/mask resources must be present and non-empty (see docs/00-reference-setup.md)
for f in gc_hg38.pytor mask_hg38.pytor; do
  if [ ! -s "${CNVPYTOR_DATA}/${f}" ]; then
    echo "ERROR: CNVpytor resource '${f}' missing or empty in ${CNVPYTOR_DATA}." >&2
    echo "       Download the pinned resource files first — see docs/00-reference-setup.md." >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
PYTOR="/genome/${SAMPLE}/cnvpytor/${SAMPLE}.pytor"

# cnvpytor invocation with the genome data + pinned resource mounts
cnvpytor_run() {
  docker run --rm --user root \
    --cpus 4 --memory 8g \
    -v "${GENOME_DIR}:/genome" \
    -v "${CNVPYTOR_DATA}:${CNVPYTOR_IMG_DATA}" \
    "${CNVPYTOR_IMAGE}" "$@"
}

echo "[1/6] Importing read depth from BAM..."
cnvpytor_run cnvpytor -root "$PYTOR" -rd "/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam" -j 4

echo "[2/6] Read-depth histogram + GC correction..."
cnvpytor_run cnvpytor -root "$PYTOR" -his "$BIN_SIZE"

echo "[3/6] Partitioning (mean-shift segmentation)..."
cnvpytor_run cnvpytor -root "$PYTOR" -partition "$BIN_SIZE"

echo "[4/6] Calling CNVs..."
cnvpytor_run cnvpytor -root "$PYTOR" -call "$BIN_SIZE" > "${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"

echo "[5/6] Exporting VCF..."
# `-view` reads its commands from stdin when stdin is not a TTY, so docker needs -i.
docker run --rm --user root -i \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  -v "${CNVPYTOR_DATA}:${CNVPYTOR_IMG_DATA}" \
  "${CNVPYTOR_IMAGE}" \
  cnvpytor -root "$PYTOR" -view "$BIN_SIZE" > /dev/null <<VIEW
set print_filename /genome/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.raw.vcf
print calls
VIEW

echo "[6/6] Normalizing VCF (full contig headers, sort, compress, index)..."
# CNVpytor's VCF only carries ##contig lines for processed chromosomes; reheader
# from the reference .fai so headers match the other SV callers for consensus
# merging (step 22). Emit a valid header-only VCF when there are no calls.
docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -c "
    set -euo pipefail
    RAW=/genome/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.raw.vcf
    OUT=/genome/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.vcf.gz
    if [ -s \"\$RAW\" ] && grep -qv '^#' \"\$RAW\"; then
      bcftools reheader --fai /genome/${REF_FAI} \"\$RAW\" | bcftools sort -Oz -o \"\$OUT\" -
    else
      { [ -s \"\$RAW\" ] && bcftools view -h \"\$RAW\" \
          || printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'; } \
        | bcftools view -Oz -o \"\$OUT\" -
    fi
    bcftools index -t \"\$OUT\"
  "
rm -f "${OUTPUT_DIR}/${SAMPLE}_cnvs.raw.vcf"

CNV_COUNT=$(wc -l < "${OUTPUT_DIR}/${SAMPLE}_cnvs.txt")
echo "=== CNVpytor complete ==="
echo "Total CNVs called: ${CNV_COUNT}"
echo "Calls (TSV): ${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"
echo "Calls (VCF): ${OUTPUT_DIR}/${SAMPLE}_cnvs.vcf.gz"
echo ""
echo "Filter significant CNVs (e-val1 < 0.05 in col 5, size > 1kb in col 3):"
echo "  awk -F'\t' '\$5 < 0.05 && \$3 > 1000' ${OUTPUT_DIR}/${SAMPLE}_cnvs.txt"
