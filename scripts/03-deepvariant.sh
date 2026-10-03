#!/usr/bin/env bash
# DeepVariant — Variant calling (BAM to VCF and gVCF)
# Usage: ./scripts/03-deepvariant.sh <sample_name> [male|female]
# Input: sorted BAM + GRCh38 reference
# Output: <sample>.vcf.gz with SNPs and small indels (~5.5M variants per 30X
#   WGS) and <sample>.g.vcf.gz, which also records the stretches where the
#   sample matches the reference (steps 07, 14 and 25 read hom-ref genotypes
#   from it), both with a .tbi index, in $GENOME_DIR/<sample>/<VCF_OUT_DIR>/.
#
# Sex: with `male`, chrX and chrY are called haploid outside the
#   pseudoautosomal regions (assets/par_grch38.bed), so a male sample gets no
#   impossible heterozygous call there. `female` and no sex call every contig
#   diploid (the alternative callers 03a, 03b and 03d always do).
# Optional environment:
#   INTERVALS="chr20:10000001-10500000 chr22:1-50818468" calls only those
#     regions (space-separated region literals, passed to --regions). Unset
#     means the whole genome.
#   ALIGN_DIR    directory of the BAM inside the sample directory (default aligned)
#   VCF_OUT_DIR  output directory inside the sample directory (default vcf).
#     Use another name for a second BAM, or its VCF replaces the primary one.
#   MODEL_TYPE   WGS (default), WES, PACBIO or ONT_R104
#   THREADS      CPUs and DeepVariant shards (default 8)
#   DV_MEM       container memory (default 32g)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [male|female]}
SEX=${2:-}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
case "$SEX" in
  ''|male|female) ;;
  *) echo "ERROR: sex must be 'male' or 'female' (or left out), got '${SEX}'." >&2; exit 2 ;;
esac
ALIGN_DIR=${ALIGN_DIR:-aligned}
VCF_OUT_DIR=${VCF_OUT_DIR:-vcf}
if [[ ! "$VCF_OUT_DIR" =~ ^[A-Za-z0-9._-]+$ ]] || [ "$VCF_OUT_DIR" = . ] || [ "$VCF_OUT_DIR" = .. ]; then
  echo "ERROR: VCF_OUT_DIR must be a directory name inside the sample directory, got '${VCF_OUT_DIR}'." >&2
  exit 2
fi
INTERVALS=${INTERVALS:-}
DV_MEM=${DV_MEM:-32g}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
BAM="${SAMPLE_DIR}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
REF="$REF_FASTA"
OUTPUT_DIR="${SAMPLE_DIR}/${VCF_OUT_DIR}"
VCF="${OUTPUT_DIR}/${SAMPLE}.vcf.gz"
GVCF="${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz"
PAR_BED="${PGP_ROOT}/assets/par_grch38.bed"
# DeepVariant writes under these names; they are renamed once both files and
# their indexes are complete, so a killed run leaves no VCF that looks finished.
PART="${OUTPUT_DIR}/${SAMPLE}.part"
TMP_DIR="${OUTPUT_DIR}/deepvariant_tmp"

# Select DeepVariant model type: WGS (default), WES, or PACBIO/ONT_R104
# WES uses a model trained on exome depth profiles and capture boundaries.
MODEL_TYPE=${MODEL_TYPE:-WGS}
case "$MODEL_TYPE" in
  WGS|WES|PACBIO|ONT_R104) ;;
  *) echo "ERROR: MODEL_TYPE must be WGS, WES, PACBIO, or ONT_R104, got '${MODEL_TYPE}'" >&2; exit 1 ;;
esac

echo "=== DeepVariant: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Model type: ${MODEL_TYPE}"
echo "Reference: ${REF}"
echo "Threads: ${THREADS}, memory: ${DV_MEM}"
if [ -n "$INTERVALS" ]; then
  echo "Regions: ${INTERVALS}"
fi
case "$SEX" in
  male) echo "Ploidy: male, chrX and chrY haploid outside the PARs (${PAR_BED#"${PGP_ROOT}"/})" ;;
  female) echo "Ploidy: female, every contig diploid" ;;
  *) echo "Ploidy: no sex given, chrX and chrY are called diploid. Pass male or female as the second argument." ;;
esac
echo "Output: ${VCF} and ${GVCF}"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$REF" "$PAR_BED"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
cleanup() { rm -rf "$TMP_DIR" "${PART}".*; }
cleanup
trap cleanup EXIT
mkdir -p "$TMP_DIR"

OUT_C=$(cpath "$OUTPUT_DIR")
DV_ARGS=(
  --model_type="${MODEL_TYPE}"
  --ref="${REF_FASTA_C}"
  --reads="/genome/${SAMPLE}/${ALIGN_DIR}/${SAMPLE}_sorted.bam"
  --output_vcf="${OUT_C}/${SAMPLE}.part.vcf.gz"
  --output_gvcf="${OUT_C}/${SAMPLE}.part.g.vcf.gz"
  --intermediate_results_dir="${OUT_C}/deepvariant_tmp"
  --sample_name="${SAMPLE}"
  --num_shards="${THREADS}"
)
if [ "$SEX" = male ]; then
  DV_ARGS+=("--haploid_contigs=chrX,chrY" --par_regions_bed=/pgp/par_grch38.bed)
fi
if [ -n "$INTERVALS" ]; then
  DV_ARGS+=(--regions "$INTERVALS")
fi

run_in \
  --cpus "${THREADS}" --memory "${DV_MEM}" \
  -v "${PAR_BED}:/pgp/par_grch38.bed:ro" \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant "${DV_ARGS[@]}"

for f in "${PART}.vcf.gz" "${PART}.g.vcf.gz"; do
  if ! have_output "$f" || [ ! -s "${f}.tbi" ]; then
    echo "ERROR: DeepVariant exited 0 but left no complete ${f} with its .tbi index." >&2
    exit 1
  fi
done
# The old indexes go first, so an index never sits next to a file it was not
# built from. The gVCF goes before the VCF: run-all.sh takes a VCF with its
# index as a finished call, so the VCF pair is the last thing to appear.
rm -f "${VCF}.tbi" "${GVCF}.tbi"
mv -f "${PART}.g.vcf.gz" "$GVCF"
mv -f "${PART}.g.vcf.gz.tbi" "${GVCF}.tbi"
mv -f "${PART}.vcf.gz" "$VCF"
mv -f "${PART}.vcf.gz.tbi" "${VCF}.tbi"
if [ -f "${PART}.visual_report.html" ]; then
  mv -f "${PART}.visual_report.html" "${OUTPUT_DIR}/${SAMPLE}.visual_report.html"
fi

echo "=== DeepVariant complete ==="
echo "VCF:  ${VCF}"
echo "gVCF: ${GVCF}"
echo ""
echo "Quick stats:"
echo "  Total variants: $(run_in "${BCFTOOLS_IMAGE}" bcftools stats "${OUT_C}/${SAMPLE}.vcf.gz" | grep '^SN' | grep 'number of records' | awk '{print $NF}' 2>/dev/null || echo 'run bcftools stats manually')"
