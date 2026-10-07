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

# CPSR 2.3 gives every variant in the panel genes its own class
# (CPSR_CLASSIFICATION), so 2.2's --classify_all is gone (2.3 refuses it).
# The final class (CLASSIFICATION, source in ASSERTION_AUTHORITY) is
# ClinVar's unless ClinVar has no record or a conflicted one: the default
# --clinvar_trust_level 0.

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
    --sample_id "${SAMPLE}" \
    --panel_id 0 \
    ${CPSR_EXTRA[@]+"${CPSR_EXTRA[@]}"} \
    --force_overwrite

echo "=== CPSR complete ==="
echo "HTML report: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.html"
echo "Variant table: ${OUTPUT_DIR}/${SAMPLE}.cpsr.grch38.classification.tsv.gz"
