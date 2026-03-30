#!/usr/bin/env bash
# 23-clinical-filter.sh — Extract clinically interesting variants from VEP-annotated VCF
# Usage: ./scripts/23-clinical-filter.sh <sample_name>
#
# Produces a small VCF of variants that are:
#   - Rare (gnomAD AF < 1%) AND
#   - Functionally impactful (HIGH/MODERATE impact, or ClinVar pathogenic/likely pathogenic)
#
# Requires: VEP-annotated VCF from step 13
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

VEP_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"
VEP_VCF_GZ="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/clinical"
mkdir -p "$OUTDIR"

# Find VEP output (may be compressed or not)
INPUT=""
if [ -f "$VEP_VCF_GZ" ]; then
  INPUT="$VEP_VCF_GZ"
elif [ -f "$VEP_VCF" ]; then
  INPUT="$VEP_VCF"
else
  echo "ERROR: VEP-annotated VCF not found. Run step 13 first."
  echo "  Expected: ${VEP_VCF} or ${VEP_VCF_GZ}"
  exit 1
fi

echo "============================================"
echo "  Step 23: Clinical Variant Filter"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${INPUT}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Determine container input path
CONTAINER_INPUT="/genome/${SAMPLE}/vep/$(basename "$INPUT")"

# If input is uncompressed, compress and index it first
if [[ "$INPUT" == *.vcf ]] && [ ! -f "$VEP_VCF_GZ" ]; then
  echo "[1/4] Compressing VEP VCF (required for bcftools filtering)..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bash -c "bcftools view /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf -Oz \
      -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  CONTAINER_INPUT="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  echo "  Done."
elif [[ "$INPUT" == *.vcf.gz ]] && [ ! -f "${INPUT}.tbi" ]; then
  echo "[1/4] Indexing compressed VEP VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bcftools index -t "$CONTAINER_INPUT"
  echo "  Done."
else
  echo "[1/4] VEP VCF already compressed and indexed."
fi

# Step 2: Extract HIGH impact variants (loss of function)
echo "[2/4] Extracting HIGH impact variants (stop-gain, frameshift, splice)..."
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools view -f PASS /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz | \
    grep -E '^#|HIGH' | \
    bcftools view -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz"

HIGH_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Found: ${HIGH_COUNT} HIGH impact variants"

# Step 3: Extract rare MODERATE impact variants (missense, in-frame indel)
echo "[3/4] Extracting rare MODERATE impact variants (gnomAD AF < 1%)..."
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools view -f PASS /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz | \
    grep -E '^#|MODERATE' | \
    grep -v 'gnomADe_AF=0\\.[1-9]' | \
    grep -v 'gnomADe_AF=0\\.0[1-9]' | \
    bcftools view -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"

MODERATE_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Found: ${MODERATE_COUNT} rare MODERATE impact variants"

# Step 4: Merge into single clinical VCF
echo "[4/4] Merging into combined clinical VCF..."
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools concat -a \
    /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz \
    /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz | \
    bcftools sort -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz"

TOTAL_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" 2>/dev/null | wc -l || echo 0)

# Generate summary TSV
echo ""
echo "Generating human-readable summary..."
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "echo -e 'CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tCSQ_EXCERPT' > /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv && \
    bcftools view -H /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz | \
    awk -F'\t' '{
      gt=\".\";
      split(\$9, fmt, \":\");
      split(\$10, vals, \":\");
      for(i in fmt) if(fmt[i]==\"GT\") gt=vals[i];
      impact=\".\";
      if(\$8 ~ /HIGH/) impact=\"HIGH\";
      else if(\$8 ~ /MODERATE/) impact=\"MODERATE\";
      csq=\".\";
      match(\$8, /CSQ=[^;]+/);
      if(RSTART>0) csq=substr(\$8, RSTART, RLENGTH>200?200:RLENGTH);
      print \$1\"\t\"\$2\"\t\"\$4\"\t\"\$5\"\t\"gt\"\t\"impact\"\t\"csq;
    }' >> /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv"

echo ""
echo "============================================"
echo "  Clinical filter complete: ${SAMPLE}"
echo "  Total clinically interesting variants: ${TOTAL_COUNT}"
echo "    HIGH impact (LoF):        ${HIGH_COUNT}"
echo "    Rare MODERATE impact:     ${MODERATE_COUNT}"
echo ""
echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_clinical.vcf.gz       (combined VCF)"
echo "    ${OUTDIR}/${SAMPLE}_clinical_summary.tsv  (human-readable table)"
echo "    ${OUTDIR}/${SAMPLE}_high_impact.vcf.gz    (HIGH only)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate.vcf.gz  (rare MODERATE only)"
echo "============================================"
echo ""
echo "Next: Review ${SAMPLE}_clinical_summary.tsv or load the VCF in IGV/gene.iobio"
