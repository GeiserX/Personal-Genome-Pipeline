#!/usr/bin/env bash
# 23-clinical-filter.sh — Extract clinically interesting variants from VEP-annotated VCF
# Usage: ./scripts/23-clinical-filter.sh <sample_name>
#
# Produces a small VCF of variants that are:
#   - Rare (gnomAD AF < 1%) AND functionally impactful (HIGH/MODERATE VEP impact)
#   - OR ClinVar pathogenic/likely pathogenic (if VEP includes ClinVar annotations)
#
# Uses bcftools +split-vep to parse VEP CSQ fields structurally.
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

# Container path for docker volume mount
CONTAINER_INPUT="/genome/${SAMPLE}/vep/$(basename "$INPUT")"

echo "============================================"
echo "  Step 23: Clinical Variant Filter"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${INPUT}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Detect available CSQ subfields (reads VCF header only, no index needed)
echo "Detecting VEP annotation fields..."
VEP_FIELDS=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools +split-vep -l "$CONTAINER_INPUT" 2>/dev/null || echo "")

if [ -z "$VEP_FIELDS" ]; then
  echo "ERROR: No CSQ/BCSQ annotation found in VEP VCF."
  echo "  Was VEP step 13 run correctly? The VCF must contain a CSQ INFO field."
  exit 1
fi

HAS_GNOMAD=0
HAS_CLINVAR=0
echo "$VEP_FIELDS" | grep -q 'gnomADe_AF' && HAS_GNOMAD=1
echo "$VEP_FIELDS" | grep -q 'CLIN_SIG' && HAS_CLINVAR=1

echo "  gnomAD frequencies: $([ "$HAS_GNOMAD" -eq 1 ] && echo 'available' || echo 'not in VEP output')"
echo "  ClinVar annotations: $([ "$HAS_CLINVAR" -eq 1 ] && echo 'available' || echo 'not in VEP output')"
echo ""

TOTAL_STEPS=4
[ "$HAS_CLINVAR" -eq 1 ] && TOTAL_STEPS=5

# Step 1: Compress and index if needed
if [[ "$INPUT" == *.vcf ]] && [ ! -f "$VEP_VCF_GZ" ]; then
  echo "[1/${TOTAL_STEPS}] Compressing VEP VCF (required for bcftools filtering)..."
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
  echo "[1/${TOTAL_STEPS}] Indexing compressed VEP VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bcftools index -t "$CONTAINER_INPUT"
  echo "  Done."
else
  echo "[1/${TOTAL_STEPS}] VEP VCF already compressed and indexed."
fi

# Step 2: Extract HIGH impact variants (loss of function)
# bcftools +split-vep parses the IMPACT subfield from VEP's pipe-delimited CSQ annotation,
# selecting only the worst consequence per variant (-s worst)
echo "[2/${TOTAL_STEPS}] Extracting HIGH impact variants (stop-gain, frameshift, splice)..."
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
    bcftools +split-vep - -c IMPACT -s worst -i 'IMPACT=\"HIGH\"' \
      -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz"

HIGH_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Found: ${HIGH_COUNT} HIGH impact variants"

# Step 3: Extract rare MODERATE impact variants (missense, in-frame indel)
# Uses gnomAD allele frequency from the CSQ field if available
if [ "$HAS_GNOMAD" -eq 1 ]; then
  echo "[3/${TOTAL_STEPS}] Extracting rare MODERATE impact variants (gnomAD AF < 1%)..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bash -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \
        -i 'IMPACT=\"MODERATE\" && (gnomADe_AF<0.01 || gnomADe_AF=\".\")' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"
else
  echo "[3/${TOTAL_STEPS}] Extracting MODERATE impact variants (no gnomAD AF available)..."
  echo "  WARNING: VEP output lacks gnomAD frequencies — including all MODERATE variants."
  echo "  Tip: Re-run VEP (step 13) with --af_gnomade for population frequency filtering."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bash -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT -s worst -i 'IMPACT=\"MODERATE\"' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"
fi

MODERATE_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Found: ${MODERATE_COUNT} MODERATE impact variants"

# Step 4: Extract ClinVar pathogenic/likely pathogenic
# Matches CLIN_SIG values containing "pathogenic" (covers both pathogenic and likely_pathogenic)
# Does NOT use -s worst: a variant is included if ANY transcript has a pathogenic ClinVar entry
CLINVAR_COUNT=0
CLINVAR_FILE=""
if [ "$HAS_CLINVAR" -eq 1 ]; then
  echo "[4/${TOTAL_STEPS}] Extracting ClinVar pathogenic/likely pathogenic variants..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bash -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\"' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz"

  CLINVAR_COUNT=$(docker run --rm \
    -v "${GENOME_DIR}:/genome" \
    staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "  Found: ${CLINVAR_COUNT} ClinVar pathogenic/likely pathogenic variants"
  CLINVAR_FILE="/genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz"
else
  echo "[4/${TOTAL_STEPS}] Skipping ClinVar filter (CLIN_SIG not in VEP annotations)."
  echo "  Tip: Re-run VEP with --everything or --check_existing to include ClinVar."
fi

# Final step: Merge into combined clinical VCF
echo "[${TOTAL_STEPS}/${TOTAL_STEPS}] Merging into combined clinical VCF..."
MERGE_FILES="/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"
[ -n "$CLINVAR_FILE" ] && MERGE_FILES="${MERGE_FILES} ${CLINVAR_FILE}"

docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "bcftools concat -a \
    ${MERGE_FILES} | \
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
      else if(\$8 ~ /pathogenic/) impact=\"CLINVAR\";
      csq=\".\";
      match(\$8, /CSQ=[^;]+/);
      if(RSTART>0) csq=substr(\$8, RSTART, RLENGTH>200?200:RLENGTH);
      print \$1\"\t\"\$2\"\t\"\$4\"\t\"\$5\"\t\"gt\"\t\"impact\"\t\"csq;
    }' >> /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv"

echo ""
echo "============================================"
echo "  Clinical filter complete: ${SAMPLE}"
echo "  Total clinically interesting variants: ${TOTAL_COUNT}"
echo "    HIGH impact (LoF):           ${HIGH_COUNT}"
echo "    Rare MODERATE impact:        ${MODERATE_COUNT}"
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    ClinVar pathogenic/LP:       ${CLINVAR_COUNT}"
fi
echo ""
echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_clinical.vcf.gz              (combined)"
echo "    ${OUTDIR}/${SAMPLE}_clinical_summary.tsv         (human-readable table)"
echo "    ${OUTDIR}/${SAMPLE}_high_impact.vcf.gz           (HIGH only)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate.vcf.gz         (rare MODERATE only)"
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    ${OUTDIR}/${SAMPLE}_clinvar_pathogenic.vcf.gz    (ClinVar P/LP only)"
fi
echo "============================================"
echo ""
echo "Next: Review ${SAMPLE}_clinical_summary.tsv or load the VCF in IGV/gene.iobio"
