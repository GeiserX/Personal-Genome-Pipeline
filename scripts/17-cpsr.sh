#!/usr/bin/env bash
# CPSR — Cancer Predisposition Sequencing Reporter
# Input: Germline VCF + PCGR 2.x data bundle + VEP cache
# Output: HTML report + classified variant TSV
# Requires: ~7 GB ref data bundle + VEP cache (download once, reuse for all samples)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
VCF_DIR=${VCF_DIR:-vcf}
VCF="${SAMPLE_DIR}/${VCF_DIR}/${SAMPLE}.vcf.gz"
VEP_DIR="${GENOME_DIR}/vep_cache"
REFDATA_DIR="${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}"
OUTPUT_DIR="${SAMPLE_DIR}/cpsr"

echo "=== CPSR Cancer Predisposition: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "VEP cache: ${VEP_DIR}"
echo "Ref data bundle: ${REFDATA_DIR}"
echo "Output: ${OUTPUT_DIR}"

# Validate inputs
if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}" >&2
  exit 1
fi

# CPSR runs the VEP release inside the PCGR image (PCGR_VEP_CACHE_RELEASE in
# versions.env), not the one step 13 uses.
if [ ! -f "${VEP_DIR}/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  CACHE_URL=$(vep_cache_url "$PCGR_VEP_CACHE_RELEASE")
  echo "ERROR: VEP release-${PCGR_VEP_CACHE_RELEASE} cache not found at ${VEP_DIR}/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/" >&2
  echo "${PCGR_IMAGE} requires the VEP ${PCGR_VEP_CACHE_RELEASE} cache (different from step 13's release-${VEP_CACHE_RELEASE})." >&2
  echo "Download and extract it:" >&2
  echo "  mkdir -p ${VEP_DIR}" >&2
  echo "  curl -fL -C - -o ${VEP_DIR}/$(basename "$CACHE_URL") ${CACHE_URL}" >&2
  echo "  tar xzf ${VEP_DIR}/$(basename "$CACHE_URL") -C ${VEP_DIR}" >&2
  exit 1
fi

if [ ! -d "${REFDATA_DIR}/data" ]; then
  echo "ERROR: PCGR 2.x ref data bundle not found at ${REFDATA_DIR}/data/" >&2
  echo "Download and extract it first:" >&2
  echo "  cd ${GENOME_DIR}/pcgr_data" >&2
  echo "  curl -fL -C - -O https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz" >&2
  echo "  tar xzf pcgr_ref_data.${PCGR_DATA_BUNDLE}.grch38.tgz" >&2
  echo "  mkdir -p ${PCGR_DATA_BUNDLE} && mv data/ ${PCGR_DATA_BUNDLE}/" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# CPSR stops on a --sample_id outside 3 to 40 characters (PCGR 2.3.2): it gets
# bin/cpsr_sample_id's id, and its files are renamed back to ${SAMPLE}.cpsr.*
# after the run. CPSR's own HTML title then shows that id.
CPSR_ID=$(sh "${PGP_ROOT}/bin/cpsr_sample_id" "$SAMPLE")
[ "$CPSR_ID" = "$SAMPLE" ] \
  || echo "CPSR takes a sample id of 3 to 40 characters: it runs as '${CPSR_ID}', and its files are renamed to ${SAMPLE}.cpsr.*"

# CPSR 2.3 gives every variant in the panel genes its own class
# (CPSR_CLASSIFICATION), so 2.2's --classify_all is gone (2.3 refuses it).
# The final class (CLASSIFICATION, source in ASSERTION_AUTHORITY) is
# ClinVar's unless ClinVar has no record or a conflicted one: the default
# --clinvar_trust_level 0. The reports (bin/collect_summary.py) count the
# final class and list every row where CPSR's own class and ClinVar's
# disagree on P/LP.

# ACMG secondary findings (genes outside CPSR's cancer panels, such as cardiac
# and metabolic ones) are reported unless CPSR_SECONDARY_FINDINGS=false: some
# people do not want to learn about them.
CPSR_EXTRA=()
case "${CPSR_SECONDARY_FINDINGS:-true}" in
  true) CPSR_EXTRA+=(--secondary_findings) ;;
  false) echo "Secondary findings: off (CPSR_SECONDARY_FINDINGS=false)" ;;
  *) echo "ERROR: CPSR_SECONDARY_FINDINGS must be true or false, got '${CPSR_SECONDARY_FINDINGS}'" >&2; exit 1 ;;
esac

# --root: the PCGR image has not been shown to run as an unprivileged user.
# The VEP cache and the PCGR bundle stay writable as before: no CI run shows
# that PCGR and its VEP never write into them.
run_in --root --cpus 4 --memory 8g \
  -v "${VEP_DIR}:/mnt/.vep" \
  -v "${REFDATA_DIR}:/mnt/bundle" \
  -v "${SAMPLE_DIR}/${VCF_DIR}:/mnt/inputs" \
  -v "${SAMPLE_DIR}/cpsr:/mnt/outputs" \
  "${PCGR_IMAGE}" \
  cpsr \
    --input_vcf "/mnt/inputs/${SAMPLE}.vcf.gz" \
    --vep_dir /mnt/.vep \
    --refdata_dir /mnt/bundle \
    --output_dir /mnt/outputs \
    --genome_assembly grch38 \
    --sample_id "${CPSR_ID}" \
    --panel_id 0 \
    ${CPSR_EXTRA[@]+"${CPSR_EXTRA[@]}"} \
    --force_overwrite
sh "${PGP_ROOT}/bin/cpsr_sample_id" --rename "$OUTPUT_DIR" "$SAMPLE"

# The cDNA read check of step 06 on CPSR's P/LP rows (final class or CPSR's
# own): spliced cDNA or expression-vector DNA in the sample gives alt reads
# clipped exactly at an exon edge, and CPSR classes a novel variant such a
# construct carries like any other. bin/clip_check.py --cpsr reads the rows
# from the classification table and writes ${SAMPLE}_cpsr_cdna_check.tsv,
# which the reports put beside each row. Same BAM and exons as step 06.
CLASS_TSV="${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.classification.tsv.gz"
CDNA_TSV="${OUTPUT_DIR}/${SAMPLE}_cpsr_cdna_check.tsv"
BAM="${SAMPLE_DIR}/${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam"
EXONS_BED=$(data_file gencode_exons) || true
if [ ! -s "$CLASS_TSV" ]; then
  rm -f "$CDNA_TSV"
  echo "cDNA read check: skipped, CPSR wrote no classification table"
else
  if [ ! -f "$BAM" ] || [ ! -f "${BAM}.bai" ]; then
    CDNA_ARGS=(--skip "no BAM (${ALIGN_DIR:-aligned}/${SAMPLE}_sorted.bam with its .bai)")
  elif [ ! -s "$EXONS_BED" ]; then
    CDNA_ARGS=(--skip "no GENCODE exons (scripts/setup.sh installs them)")
  else
    CDNA_ARGS=(--bam "$(cpath "$BAM")" --exons "$(cpath "$EXONS_BED")")
  fi
  # pysam is in PYPGX_IMAGE.
  run_in --cpus 1 --memory 2g -v "${PGP_ROOT}/bin:/pgp-bin:ro" "$PYPGX_IMAGE" \
    python3 /pgp-bin/clip_check.py --cpsr "$(cpath "$CLASS_TSV")" --out "$(cpath "${CDNA_TSV}.part")" "${CDNA_ARGS[@]}"
  mv "${CDNA_TSV}.part" "$CDNA_TSV"
fi

echo "=== CPSR complete ==="
echo "HTML report: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.html"
echo "Variant table: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.classification.tsv.gz"
[ ! -f "$CDNA_TSV" ] || echo "cDNA read check: ${CDNA_TSV}"
