#!/usr/bin/env bash
# PharmCAT — Clinical pharmacogenomics (star alleles + drug recommendations)
# Input: VCF.gz + GRCh38 reference; the gVCF from step 03 when there is one;
#        the outside calls of step 36 (HLA-A, HLA-B, an agreed CYP2D6) when
#        that step wrote any
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

# Step 0: PharmCAT reads a PGx position missing from its input as "not
# covered", and a variants-only VCF lists only where the sample differs from
# the reference. Step 03's gVCF also records where it matches, so when there
# is one its reference blocks are expanded over PharmCAT's gene regions into a
# plain VCF: every covered position becomes a 0/0 call, and an uncovered one
# (./.) stays missing. PharmCAT refuses a gVCF, by content and by a name with
# .g.vcf in it, so the expanded file has neither.
GVCF="${OUTPUT_DIR}/${SAMPLE}.g.vcf.gz"
PGX_INPUT="${SAMPLE}.vcf.gz"
rm -f "${OUTPUT_DIR}/${SAMPLE}.pgx_regions.vcf.gz" "${OUTPUT_DIR}/${SAMPLE}.pgx_regions.vcf.gz.tbi"
if [ -f "$GVCF" ] && [ -f "${GVCF}.tbi" ]; then
  echo "Input: ${GVCF} (reference blocks expanded over PharmCAT's gene regions)"
  PGX_INPUT="${SAMPLE}.pgx_regions.vcf.gz"
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  run_in \
    --cpus 1 --memory 2g \
    -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
    "${PHARMCAT_IMAGE}" \
    bash -euo pipefail -c '
      bcftools convert --gvcf2vcf -f "$1" -R /pharmcat/pharmcat_regions.bed -Ou "$2" \
        | bcftools view --trim-alt-alleles -i "GT!=\"mis\"" -Oz -o "$3.part" --write-index=tbi
      mv -f "$3.part" "$3"
      mv -f "$3.part.tbi" "$3.tbi"' \
    _ "${REF_FASTA_C}" "/data/${SAMPLE}.g.vcf.gz" "/data/${PGX_INPUT}"
else
  echo "Input: ${VCF} (no gVCF: PGx positions where you match the reference read as missing)"
fi

# Step 1: Preprocess VCF (normalize, filter to PGx positions)
run_in \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  -v "$(dirname "$REF_FASTA"):/ref:ro" \
  "${PHARMCAT_IMAGE}" \
  python3 /pharmcat/pharmcat_vcf_preprocessor \
    -vcf "/data/${PGX_INPUT}" \
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

# Step 3: Run PharmCAT on preprocessed VCF. PharmCAT types neither HLA nor
# CYP2D6 from a VCF; step 36 writes those calls from the BAM-based callers
# (T1K, and CYP2D6 only when pypgx and Cyrius agree), and PharmCAT reads them
# with -po. An empty file means nothing was agreed, and is not passed.
OUTSIDE="${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_outside_calls.tsv"
OUTSIDE_ARGS=() PO_ARGS=()
if [ -s "$OUTSIDE" ]; then
  echo "Outside calls (step 36): ${OUTSIDE}"
  sed 's/^/  /' "$OUTSIDE"
  OUTSIDE_ARGS=(-v "${OUTSIDE}:/outside_calls.tsv:ro")
  PO_ARGS=(-po /outside_calls.tsv)
else
  echo "No outside calls (step 36 not run, or it passed none): HLA and CYP2D6 stay uncalled."
fi
run_in \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}/${SAMPLE}/vcf:/data" \
  ${OUTSIDE_ARGS[@]+"${OUTSIDE_ARGS[@]}"} \
  "${PHARMCAT_IMAGE}" \
  java -jar /pharmcat/pharmcat.jar \
    -vcf "/data/${SAMPLE}.pharmcat_input.vcf" \
    ${PO_ARGS[@]+"${PO_ARGS[@]}"} \
    -o /data/ \
    -bf "$SAMPLE" \
    -reporterJson \
    -reporterHtml

rm -f "${OUTPUT_DIR}/${SAMPLE}.pharmcat_input.vcf" "${OUTPUT_DIR}/${SAMPLE}.pgx_regions.vcf.gz" \
  "${OUTPUT_DIR}/${SAMPLE}.pgx_regions.vcf.gz.tbi"

echo "=== PharmCAT complete ==="
echo "Reports: ${OUTPUT_DIR}/${SAMPLE}.report.html and ${OUTPUT_DIR}/${SAMPLE}.report.json"
echo ""
echo "Key genes covered: CYP2C19, CYP2D6, CYP2B6, CYP3A5, UGT1A1, DPYD, NAT2, TPMT"
echo "NOTE: PharmCAT calls no CYP2D6 or HLA from a VCF. Step 36 gives it T1K's HLA types and a"
echo "  CYP2D6 call that pypgx (step 32) and Cyrius (step 21) agree on; run it, then this step again."
