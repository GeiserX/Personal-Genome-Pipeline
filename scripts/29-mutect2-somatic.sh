#!/usr/bin/env bash
# [EXPERIMENTAL] Somatic variant calling — Mutect2 tumor-only mode
# Finds mutations acquired during life (not inherited), such as clonal hematopoiesis
# Input: Sorted BAM + GRCh38 reference
# Output: Filtered somatic VCF in $GENOME_DIR/<sample>/somatic/
# Runtime: minutes on the default CHIP genes, ~2-6 hours with INTERVALS=genome
#
# By default Mutect2 calls only the clonal hematopoiesis (CHIP) driver genes in
# assets/chip_genes_grch38.bed: gene bodies from GENCODE 50's basic annotation,
# 100 bp either side. INTERVALS=genome calls the whole genome, and any other
# INTERVALS value (a region such as chr22, or a BED under /genome) is passed
# to Mutect2 as it is. The read orientation model (LearnReadOrientationModel)
# always runs; the contamination estimate (GetPileupSummaries,
# CalculateContamination) runs when the common-sites VCF is in somatic/.
#
# WARNING: Tumor-only mode (no matched normal) has a HIGH false positive rate.
# Many germline variants will be called as somatic. Use gnomAD and PoN resources
# to reduce false positives, and treat results as exploratory only.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
THREADS=${THREADS:-4}   # common.sh defaults to 8
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
INTERVALS=${INTERVALS:-chip}

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
ALIGN_DIR=${ALIGN_DIR:-aligned}
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/somatic"


# Optional resources (improve filtering if present)
GNOMAD_VCF="${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz"
PON_VCF="${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz"
COMMON_VCF="${GENOME_DIR}/somatic/small_exac_common_3.hg38.vcf.gz"
# has_vcf FILE: "yes" when FILE and its .tbi are present, else "no".
has_vcf() { if [ -f "$1" ] && [ -f "${1}.tbi" ]; then echo yes; else echo no; fi; }
HAVE_GNOMAD=$(has_vcf "$GNOMAD_VCF")
HAVE_PON=$(has_vcf "$PON_VCF")
HAVE_COMMON=$(has_vcf "$COMMON_VCF")

echo "=== [EXPERIMENTAL] Somatic Variant Calling (Mutect2 Tumor-Only): ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}"
case "$INTERVALS" in
  chip) echo "Intervals: CHIP driver genes (assets/chip_genes_grch38.bed; INTERVALS=genome for the whole genome)" ;;
  genome) echo "Intervals: whole genome" ;;
  *) echo "Intervals: ${INTERVALS}" ;;
esac
echo "Output: ${OUTPUT_DIR}/"
echo ""

# Check for idempotent skip. A finished result is reused only when it was
# called the way this run would call it: RUN_FILE holds the INTERVALS value and
# which optional resources (gnomAD, Panel of Normals, common sites) were there,
# and is written last. A CHIP result never stands in for a whole-genome
# request, a result filtered without the contamination estimate is redone once
# the common-sites VCF is installed, and a result without the record (from an
# older version of this step) is called again.
FINAL_OUTPUT="${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz"
RUN_FILE="${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.run"
RUN_KEY="INTERVALS=${INTERVALS} germline_resource=${HAVE_GNOMAD} panel_of_normals=${HAVE_PON} common_sites=${HAVE_COMMON}"
if have_output "$FINAL_OUTPUT" && [ -f "$RUN_FILE" ] && [ "$(cat "$RUN_FILE")" = "$RUN_KEY" ]; then
  echo "Output already exists: ${FINAL_OUTPUT} (${RUN_KEY})"
  echo "Skipping. Delete the file to re-run."
  exit 0
fi
if [ -e "$FINAL_OUTPUT" ]; then
  echo "Calling again: ${FINAL_OUTPUT} is unfinished or was not called with ${RUN_KEY}."
fi
rm -f "$RUN_FILE"

# Validate required inputs
for f in "$BAM" "${BAM}.bai" "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# GATK needs .dict file
if [ ! -f "$REF_DICT" ]; then
  echo "ERROR: Sequence dictionary not found: ${REF_DICT}" >&2
  echo "Generate it with: docker run --rm -v \"${GENOME_DIR}:/genome\" ${GATK_IMAGE} gatk CreateSequenceDictionary -R ${REF_FASTA_C}" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
S="/genome/${SAMPLE}/somatic/${SAMPLE}"

# The CHIP BED is copied into the sample's directory: containers see only GENOME_DIR.
INTERVAL_ARGS=()
case "$INTERVALS" in
  chip)
    cp "${PGP_ROOT}/assets/chip_genes_grch38.bed" "${OUTPUT_DIR}/chip_genes_grch38.bed"
    INTERVAL_ARGS=(--intervals "/genome/${SAMPLE}/somatic/chip_genes_grch38.bed") ;;
  genome) ;;
  *) INTERVAL_ARGS=(--intervals "$INTERVALS") ;;
esac

# Build Mutect2 command. --f1r2-tar-gz collects the read orientation counts
# that LearnReadOrientationModel turns into priors for FilterMutectCalls.
MUTECT2_CMD=(
  gatk Mutect2
  -R "${REF_FASTA_C}"
  -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  -O "${S}_somatic_unfiltered.vcf.gz"
  --f1r2-tar-gz "${S}_f1r2.tar.gz"
  --native-pair-hmm-threads "$THREADS"
  --max-mnp-distance 0
)

# Add gnomAD germline resource if available (reduces germline false positives)
if [ "$HAVE_GNOMAD" = yes ]; then
  echo "Using gnomAD germline resource: ${GNOMAD_VCF}"
  MUTECT2_CMD+=(--germline-resource "/genome/somatic/af-only-gnomad.hg38.vcf.gz")
else
  echo "WARNING: gnomAD AF-only VCF not found at ${GNOMAD_VCF}"
  echo "  Without gnomAD, many common germline variants will appear as somatic calls."
  echo "  See docs/29-mutect2-somatic.md for download instructions."
  echo ""
fi

# Add Panel of Normals if available (reduces recurrent technical artifacts)
if [ "$HAVE_PON" = yes ]; then
  echo "Using Panel of Normals: ${PON_VCF}"
  MUTECT2_CMD+=(-pon "/genome/somatic/1000g_pon.hg38.vcf.gz")
else
  echo "INFO: Panel of Normals not found at ${PON_VCF} (optional, reduces artifacts)."
  echo ""
fi

MUTECT2_CMD+=(${INTERVAL_ARGS[@]+"${INTERVAL_ARGS[@]}"})

echo "=== [1/4] Running Mutect2 in tumor-only mode ==="
echo "  The slowest step: minutes on the CHIP genes, ~2-6 hours with INTERVALS=genome."
echo ""
run_in --cpus "$THREADS" --memory 8g \
  "$GATK_IMAGE" \
  "${MUTECT2_CMD[@]}"

echo ""
echo "=== [2/4] Read orientation model (LearnReadOrientationModel) ==="
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  gatk LearnReadOrientationModel \
    -I "${S}_f1r2.tar.gz" \
    -O "${S}_read-orientation-model.tar.gz"

FILTER_CMD=(
  gatk FilterMutectCalls
  -R "${REF_FASTA_C}"
  -V "${S}_somatic_unfiltered.vcf.gz"
  --ob-priors "${S}_read-orientation-model.tar.gz"
  -O "${S}_somatic_filtered.vcf.gz"
)

echo ""
echo "=== [3/4] Contamination (GetPileupSummaries, CalculateContamination) ==="
if [ "$HAVE_COMMON" = yes ]; then
  # Pileups at common sites only, inside the same intervals as the calls.
  PILEUP_L=(-L "$(cpath "$COMMON_VCF")")
  if [ "${#INTERVAL_ARGS[@]}" -gt 0 ]; then
    PILEUP_L+=(-L "${INTERVAL_ARGS[1]}" --interval-set-rule INTERSECTION)
  fi
  run_in --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk GetPileupSummaries \
      -I "/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam" \
      -V "$(cpath "$COMMON_VCF")" \
      "${PILEUP_L[@]}" \
      -O "${S}_pileups.table"
  run_in --cpus 2 --memory 4g \
    "$GATK_IMAGE" \
    gatk CalculateContamination \
      -I "${S}_pileups.table" \
      --tumor-segmentation "${S}_segments.table" \
      -O "${S}_contamination.table"
  FILTER_CMD+=(--contamination-table "${S}_contamination.table" --tumor-segmentation "${S}_segments.table")
else
  echo "  Skipped: no common-sites VCF at ${COMMON_VCF} (see docs/29-mutect2-somatic.md)."
fi

echo ""
echo "=== [4/4] Filtering somatic calls (FilterMutectCalls) ==="
run_in --cpus 2 --memory 4g \
  "$GATK_IMAGE" \
  "${FILTER_CMD[@]}"
printf '%s\n' "$RUN_KEY" > "$RUN_FILE"

echo ""
echo "=== Somatic variant statistics ==="

# Count PASS variants
PASS_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools view -f PASS "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz" \
  2>/dev/null | grep -c "^[^#]" || true)

TOTAL_COUNT=$(run_in \
  "$BCFTOOLS_IMAGE" \
  bcftools view "/genome/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz" \
  2>/dev/null | grep -c "^[^#]" || true)

echo "  Total calls: ${TOTAL_COUNT}"
echo "  PASS calls:  ${PASS_COUNT}"
echo ""

echo "=== [EXPERIMENTAL] Somatic variant calling complete ==="
echo ""
echo "Output files:"
echo "  Unfiltered: ${OUTPUT_DIR}/${SAMPLE}_somatic_unfiltered.vcf.gz"
echo "  Filtered:   ${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz"
echo "  Stats:      ${OUTPUT_DIR}/${SAMPLE}_somatic_unfiltered.vcf.gz.stats"
echo ""
echo "IMPORTANT: Tumor-only mode produces many false positives."
echo "  - Most PASS variants in a healthy individual are germline, not somatic."
echo "  - True somatic variants (e.g., clonal hematopoiesis) typically have low AF (<0.1)."
echo "  - Cross-reference with ClinVar and gnomAD before interpreting any variant."
echo ""
echo "View low-AF PASS variants (potential somatic):"
echo "  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' -i 'FILTER=\"PASS\"' \\"
echo "    ${OUTPUT_DIR}/${SAMPLE}_somatic_filtered.vcf.gz | awk '\$5 < 0.1'"
