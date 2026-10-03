#!/usr/bin/env bash
# PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
# Input: VCF.gz + GRCh38 reference
# Output: HTML + JSON reports with metabolizer status for 23 pharmacogenes
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
REF="$REF_FASTA"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/vcf"

echo "=== PharmCAT: ${SAMPLE} ==="
echo "Input VCF: ${VCF}"
echo "Outputs: ${OUTPUT_DIR}/${SAMPLE}.report.html and ${OUTPUT_DIR}/${SAMPLE}.report.json"

# Validate inputs (PharmCAT requires both VCF and its tabix index)
for f in "$VCF" "${VCF}.tbi" "$REF"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    if [ "$f" = "${VCF}.tbi" ]; then
      echo "  PharmCAT requires a tabix index. Generate it with:" >&2
      echo "  docker run --rm -v \"${GENOME_DIR}:/genome\" ${BCFTOOLS_IMAGE} bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" >&2
    fi
    exit 1
  fi
done

# Step 1: Preprocess VCF (normalize, filter to PGx positions)
run_in \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "$(dirname "$REF_FASTA"):/ref:ro" \
  "${PHARMCAT_IMAGE}" \
  python3 /pharmcat/pharmcat_vcf_preprocessor \
    -vcf "/data/${SAMPLE}.vcf.gz" \
    -refFna "/ref/$(basename "$REF_FASTA")" \
    -o /data/ \
    -bf "$SAMPLE"

# Step 2: PharmCAT up to 3.4.0 bundles vcf-parser 0.3.1, which stops with
# "Error parsing metadata: character to be escaped is missing" on a backslash
# in a ## header line. Such lines are valid VCF: bcftools writes one for a soft
# filter with a quoted string (-s LowDP -e 'GT!="0/0"'). Rewrite PharmCAT's own
# copy only, on ## lines: \" becomes ' and any other \ becomes /. Remove this
# once a PharmCAT release bundles vcf-parser newer than 0.3.1 (the PHARMCAT
# module in modules/local/pharmcat/main.nf does the same).
HEADER_FIX='/^##/ { gsub(/\\"/, "\047"); gsub(/\\/, "/") } { print }'
# shellcheck disable=SC2016  # $1 to $3 belong to the inner sh
run_in \
  --cpus 1 --memory 1g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${PHARMCAT_IMAGE}" \
  sh -c 'gzip -dc "$1" | awk "$3" > "$2"' sh \
    "/data/${SAMPLE}.preprocessed.vcf.bgz" "/data/${SAMPLE}.pharmcat_input.vcf" "$HEADER_FIX"

# Step 3: Run PharmCAT on preprocessed VCF
run_in \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  "${PHARMCAT_IMAGE}" \
  java -jar /pharmcat/pharmcat.jar \
    -vcf "/data/${SAMPLE}.pharmcat_input.vcf" \
    -o /data/ \
    -bf "$SAMPLE" \
    -reporterJson \
    -reporterHtml

rm -f "${OUTPUT_DIR}/${SAMPLE}.pharmcat_input.vcf"

echo "=== PharmCAT complete ==="
echo "Reports: ${OUTPUT_DIR}/${SAMPLE}.report.html and ${OUTPUT_DIR}/${SAMPLE}.report.json"
echo ""
echo "Key genes covered: CYP2C19, CYP2D6, CYP2B6, CYP3A5, UGT1A1, DPYD, NAT2, TPMT"
echo "NOTE: CYP2D6 may return 'Not called' — use Cyrius (BAM-based) for CYP2D6."
