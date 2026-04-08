#!/usr/bin/env bash
# 23-clinical-filter.sh — Extract clinically interesting variants from annotated VCF
# Usage: ./scripts/23-clinical-filter.sh <sample_name>
#
# Produces a small VCF of variants that are:
#   - Rare (gnomAD AF < 1%) AND functionally impactful (HIGH/MODERATE VEP impact)
#   - OR ClinVar pathogenic/likely pathogenic (if VEP includes ClinVar annotations)
#   - OR high CADD score (>= 20) for non-coding variants (if step 30 was run)
#   - OR high SpliceAI delta score (>= 0.2) for cryptic splice variants
#   - OR high REVEL/AlphaMissense for missense variants
#
# Uses bcftools +split-vep to parse VEP CSQ fields and bcftools view -i for
# INFO-level annotations from vcfanno (step 30).
# Requires: VEP-annotated VCF from step 13 (step 30 vcfanno enrichment recommended)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

# Validate sample name to prevent shell injection in bash -c strings
if [[ "$SAMPLE" =~ [^a-zA-Z0-9._-] ]]; then
  echo "ERROR: Sample name contains invalid characters. Use only a-z, A-Z, 0-9, ., _, -" >&2
  exit 1
fi

# Prefer vcfanno-enriched VCF (step 30), fall back to VEP VCF (step 13)
ANNOTATED_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz"
VEP_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"
VEP_VCF_GZ="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/clinical"
CONSTRAINT_TSV="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"
mkdir -p "$OUTDIR"

# Find best available input VCF
INPUT=""
if [ -f "$ANNOTATED_VCF" ]; then
  INPUT="$ANNOTATED_VCF"
elif [ -f "$VEP_VCF_GZ" ]; then
  INPUT="$VEP_VCF_GZ"
elif [ -f "$VEP_VCF" ]; then
  INPUT="$VEP_VCF"
else
  echo "ERROR: No annotated VCF found. Run step 13 (VEP) first."
  echo "  Expected: ${ANNOTATED_VCF} or ${VEP_VCF_GZ} or ${VEP_VCF}"
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
  "${BCFTOOLS_IMAGE}" \
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

# Detect vcfanno INFO fields (from step 30)
VCF_HEADER=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -h "$CONTAINER_INPUT" 2>/dev/null || echo "")

HAS_CADD=0
HAS_CADD_INDEL=0
HAS_SPLICEAI=0
HAS_SPLICEAI_INDEL=0
HAS_REVEL=0
HAS_ALPHAMISSENSE=0
echo "$VCF_HEADER" | grep -q 'ID=CADD_PHRED,' && HAS_CADD=1
echo "$VCF_HEADER" | grep -q 'ID=CADD_PHRED_indel,' && HAS_CADD_INDEL=1
echo "$VCF_HEADER" | grep -q 'ID=SpliceAI,' && HAS_SPLICEAI=1
echo "$VCF_HEADER" | grep -q 'ID=SpliceAI_indel,' && HAS_SPLICEAI_INDEL=1
echo "$VCF_HEADER" | grep -q 'ID=REVEL' && HAS_REVEL=1
echo "$VCF_HEADER" | grep -q 'ID=AM_pathogenicity' && HAS_ALPHAMISSENSE=1

HAS_CONSTRAINT=0
[ -f "$CONSTRAINT_TSV" ] && HAS_CONSTRAINT=1

echo "  gnomAD frequencies: $([ "$HAS_GNOMAD" -eq 1 ] && echo 'available' || echo 'not in VEP output')"
echo "  ClinVar annotations: $([ "$HAS_CLINVAR" -eq 1 ] && echo 'available' || echo 'not in VEP output')"
echo "  CADD scores: $([ "$HAS_CADD" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')$([ "$HAS_CADD_INDEL" -eq 1 ] && echo ' (+indels)')"
echo "  SpliceAI scores: $([ "$HAS_SPLICEAI" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')$([ "$HAS_SPLICEAI_INDEL" -eq 1 ] && echo ' (+indels)')"
echo "  REVEL scores: $([ "$HAS_REVEL" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')"
echo "  AlphaMissense: $([ "$HAS_ALPHAMISSENSE" -eq 1 ] && echo 'available' || echo 'not annotated (run step 30)')"
echo "  gnomAD constraint: $([ "$HAS_CONSTRAINT" -eq 1 ] && echo 'available' || echo 'not downloaded')"
echo ""

# Calculate total steps dynamically
TOTAL_STEPS=4
[ "$HAS_CLINVAR" -eq 1 ] && TOTAL_STEPS=$((TOTAL_STEPS + 1))
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && TOTAL_STEPS=$((TOTAL_STEPS + 1))
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && TOTAL_STEPS=$((TOTAL_STEPS + 1))
REVEL_OR_AM=0
{ [ "$HAS_REVEL" -eq 1 ] || [ "$HAS_ALPHAMISSENSE" -eq 1 ]; } && REVEL_OR_AM=1 && TOTAL_STEPS=$((TOTAL_STEPS + 1))

# Step 1: Compress and index if needed
if [[ "$INPUT" == *.vcf ]] && [ ! -f "$VEP_VCF_GZ" ]; then
  echo "[1/${TOTAL_STEPS}] Compressing VEP VCF (required for bcftools filtering)..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
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
    "${BCFTOOLS_IMAGE}" \
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
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
    bcftools +split-vep - -c IMPACT -s worst -i 'IMPACT=\"HIGH\"' \
      -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz"

HIGH_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Found: ${HIGH_COUNT} HIGH impact variants"

# Step 3: Extract rare MODERATE impact variants (missense, in-frame indel)
# Uses gnomAD allele frequency from the CSQ field if available
if [ "$HAS_GNOMAD" -eq 1 ]; then
  echo "[3/${TOTAL_STEPS}] Extracting rare MODERATE impact variants (gnomAD AF < 1%)..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
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
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT -s worst -i 'IMPACT=\"MODERATE\"' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"
fi

MODERATE_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
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
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz"

  CLINVAR_COUNT=$(docker run --rm \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "  Found: ${CLINVAR_COUNT} ClinVar pathogenic/likely pathogenic variants"
  CLINVAR_FILE="/genome/${SAMPLE}/clinical/${SAMPLE}_clinvar_pathogenic.vcf.gz"
else
  echo "[4/${TOTAL_STEPS}] Skipping ClinVar filter (CLIN_SIG not in VEP annotations)."
  echo "  Tip: Re-run VEP with --everything or --check_existing to include ClinVar."
fi

# Step N: CADD high-score non-coding variants (if annotated by step 30)
STEP_NUM=5
[ "$HAS_CLINVAR" -eq 1 ] && STEP_NUM=$((STEP_NUM + 1))
CADD_COUNT=0
CADD_FILE=""
if [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; then
  # Build CADD filter — only reference tags that exist in the header
  CADD_EXPR=''
  [ "$HAS_CADD" -eq 1 ] && CADD_EXPR='INFO/CADD_PHRED>=20'
  if [ "$HAS_CADD_INDEL" -eq 1 ]; then
    [ -n "$CADD_EXPR" ] && CADD_EXPR="${CADD_EXPR} || INFO/CADD_PHRED_indel>=20" || CADD_EXPR='INFO/CADD_PHRED_indel>=20'
  fi

  echo "[${STEP_NUM}/${TOTAL_STEPS}] Extracting high-CADD variants (PHRED >= 20, non-HIGH/MODERATE)..."
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT -s worst \
        -i 'IMPACT!=\"HIGH\" && IMPACT!=\"MODERATE\" && (${CADD_EXPR})' \
        -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_cadd_high.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_cadd_high.vcf.gz"

  CADD_COUNT=$(docker run --rm \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_cadd_high.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "  Found: ${CADD_COUNT} high-CADD non-coding variants (PHRED >= 20)"
  CADD_FILE="/genome/${SAMPLE}/clinical/${SAMPLE}_cadd_high.vcf.gz"
  STEP_NUM=$((STEP_NUM + 1))
fi

# Step N: SpliceAI cryptic splice variants (if annotated by step 30)
SPLICEAI_COUNT=0
SPLICEAI_FILE=""
if [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
  # Build SpliceAI pre-filter — only reference tags that exist in the header
  SPLICEAI_PREFILTER_PARTS=()
  [ "$HAS_SPLICEAI" -eq 1 ] && SPLICEAI_PREFILTER_PARTS+=('INFO/SpliceAI!="."')
  [ "$HAS_SPLICEAI_INDEL" -eq 1 ] && SPLICEAI_PREFILTER_PARTS+=('INFO/SpliceAI_indel!="."')
  SPLICEAI_PREFILTER=$(IFS=' || '; echo "${SPLICEAI_PREFILTER_PARTS[*]}")

  echo "[${STEP_NUM}/${TOTAL_STEPS}] Extracting cryptic splice variants (SpliceAI delta >= 0.2)..."
  # SpliceAI INFO field from vcfanno is a pipe-delimited string:
  #   ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
  # bcftools cannot numerically compare sub-fields within a string, so we use
  # awk to parse the SpliceAI value and check if any delta score >= 0.2.
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS -i '${SPLICEAI_PREFILTER}' ${CONTAINER_INPUT} | \
      awk -F'\t' 'BEGIN{OFS=\"\t\"} /^#/{print;next} {
        dominated=0
        n=split(\$8, info_arr, \";\")
        for(i=1;i<=n;i++){
          if(info_arr[i] ~ /^SpliceAI=/){
            sub(/^SpliceAI=/,\"\",info_arr[i])
            split(info_arr[i],sp,\"|\")
            for(j=3;j<=6;j++) if(sp[j]+0>=0.2) dominated=1
          }
          if(info_arr[i] ~ /^SpliceAI_indel=/){
            sub(/^SpliceAI_indel=/,\"\",info_arr[i])
            split(info_arr[i],sp,\"|\")
            for(j=3;j<=6;j++) if(sp[j]+0>=0.2) dominated=1
          }
        }
        if(dominated) print
      }' | \
      bgzip -c > /genome/${SAMPLE}/clinical/${SAMPLE}_spliceai_high.vcf.gz && \
    tabix -p vcf /genome/${SAMPLE}/clinical/${SAMPLE}_spliceai_high.vcf.gz"

  SPLICEAI_COUNT=$(docker run --rm \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_spliceai_high.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "  Found: ${SPLICEAI_COUNT} cryptic splice variants (SpliceAI >= 0.2)"
  SPLICEAI_FILE="/genome/${SAMPLE}/clinical/${SAMPLE}_spliceai_high.vcf.gz"
  STEP_NUM=$((STEP_NUM + 1))
fi

# Step N: High-confidence deleterious missense (REVEL/AlphaMissense)
MISSENSE_COUNT=0
MISSENSE_FILE=""
if [ "$REVEL_OR_AM" -eq 1 ]; then
  echo "[${STEP_NUM}/${TOTAL_STEPS}] Extracting high-confidence deleterious missense variants..."
  # Build filter expression dynamically based on available annotations
  MISSENSE_FILTER=""
  [ "$HAS_REVEL" -eq 1 ] && MISSENSE_FILTER="INFO/REVEL>=0.644"
  if [ "$HAS_ALPHAMISSENSE" -eq 1 ]; then
    if [ -n "$MISSENSE_FILTER" ]; then
      MISSENSE_FILTER="${MISSENSE_FILTER} || INFO/AM_pathogenicity>=0.564"
    else
      MISSENSE_FILTER="INFO/AM_pathogenicity>=0.564"
    fi
  fi

  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view -f PASS -i '${MISSENSE_FILTER}' \
      ${CONTAINER_INPUT} \
      -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_missense_deleterious.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_missense_deleterious.vcf.gz"

  MISSENSE_COUNT=$(docker run --rm \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_missense_deleterious.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "  Found: ${MISSENSE_COUNT} high-confidence deleterious missense variants"
  MISSENSE_FILE="/genome/${SAMPLE}/clinical/${SAMPLE}_missense_deleterious.vcf.gz"
fi

# Final step: Merge into combined clinical VCF
echo "[${TOTAL_STEPS}/${TOTAL_STEPS}] Merging into combined clinical VCF..."
MERGE_FILES="/genome/${SAMPLE}/clinical/${SAMPLE}_high_impact.vcf.gz /genome/${SAMPLE}/clinical/${SAMPLE}_rare_moderate.vcf.gz"
[ -n "$CLINVAR_FILE" ] && MERGE_FILES="${MERGE_FILES} ${CLINVAR_FILE}"
[ -n "$CADD_FILE" ] && MERGE_FILES="${MERGE_FILES} ${CADD_FILE}"
[ -n "$SPLICEAI_FILE" ] && MERGE_FILES="${MERGE_FILES} ${SPLICEAI_FILE}"
[ -n "$MISSENSE_FILE" ] && MERGE_FILES="${MERGE_FILES} ${MISSENSE_FILE}"

docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools concat -a -D \
    ${MERGE_FILES} | \
    bcftools sort -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz"

TOTAL_COUNT=$(docker run --rm \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" 2>/dev/null | wc -l || echo 0)

# Generate summary TSV with annotation scores
echo ""
echo "Generating human-readable summary..."
docker run --rm --user root \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -c "echo -e 'CHROM\tPOS\tREF\tALT\tGT\tIMPACT\tGENE\tCADD_PHRED\tREVEL\tAM_CLASS\tCSQ_EXCERPT' > /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv && \
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
      else if(\$8 ~ /CADD_PHRED/) impact=\"CADD\";
      else if(\$8 ~ /SpliceAI/) impact=\"SPLICEAI\";
      gene=\".\";
      if(match(\$8, /SYMBOL=[^;|]+/)) gene=substr(\$8, RSTART+7, RLENGTH-7);
      cadd=\".\";
      if(match(\$8, /CADD_PHRED=[^;]+/)) cadd=substr(\$8, RSTART+11, RLENGTH-11);
      revel=\".\";
      if(match(\$8, /REVEL=[^;]+/)) revel=substr(\$8, RSTART+6, RLENGTH-6);
      am=\".\";
      if(match(\$8, /AM_class=[^;]+/)) am=substr(\$8, RSTART+9, RLENGTH-9);
      csq=\".\";
      match(\$8, /CSQ=[^;]+/);
      if(RSTART>0) csq=substr(\$8, RSTART, RLENGTH>150?150:RLENGTH);
      print \$1\"\t\"\$2\"\t\"\$4\"\t\"\$5\"\t\"gt\"\t\"impact\"\t\"gene\"\t\"cadd\"\t\"revel\"\t\"am\"\t\"csq;
    }' >> /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv"

# Add gnomAD gene constraint columns if available
if [ "$HAS_CONSTRAINT" -eq 1 ]; then
  echo "Adding gnomAD gene constraint metrics..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${PYTHON_IMAGE}" \
    python3 -c "
import csv, sys

# Load constraint metrics (gene -> {loeuf, pli, mis_z})
constraint = {}
with open('/genome/annotations/gnomad_v4.1_constraint.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row.get('gene', '')
        if not gene or row.get('canonical', '') != 'true':
            continue
        try:
            loeuf = row.get('lof.oe_ci.upper', '.')
            pli = row.get('lof.pLI', '.')
            mis_z = row.get('mis.z_score', '.')
        except (KeyError, ValueError):
            loeuf, pli, mis_z = '.', '.', '.'
        constraint[gene] = (loeuf, pli, mis_z)

# Read TSV, add constraint columns
infile = '/genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv'
outfile = '/genome/${SAMPLE}/clinical/${SAMPLE}_clinical_summary_enriched.tsv'
with open(infile) as fin, open(outfile, 'w') as fout:
    header = fin.readline().rstrip('\n')
    fout.write(header + '\tLOEUF\tpLI\tmis_Z\n')
    for line in fin:
        line = line.rstrip('\n')
        cols = line.split('\t')
        gene = cols[6] if len(cols) > 6 else '.'
        loeuf, pli, mis_z = constraint.get(gene, ('.', '.', '.'))
        fout.write(line + '\t' + str(loeuf) + '\t' + str(pli) + '\t' + str(mis_z) + '\n')

# Replace original with enriched version
import shutil
shutil.move(outfile, infile)
print(f'  Added constraint metrics for genes in summary TSV')
"
fi

echo ""
echo "============================================"
echo "  Clinical filter complete: ${SAMPLE}"
echo "  Total clinically interesting variants: ${TOTAL_COUNT}"
echo "    HIGH impact (LoF):           ${HIGH_COUNT}"
echo "    Rare MODERATE impact:        ${MODERATE_COUNT}"
[ "$HAS_CLINVAR" -eq 1 ] && \
echo "    ClinVar pathogenic/LP:       ${CLINVAR_COUNT}"
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && \
echo "    High CADD non-coding:        ${CADD_COUNT}"
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && \
echo "    Cryptic splice (SpliceAI):   ${SPLICEAI_COUNT}"
[ "$REVEL_OR_AM" -eq 1 ] && \
echo "    Deleterious missense:        ${MISSENSE_COUNT}"
echo ""
echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_clinical.vcf.gz              (combined)"
echo "    ${OUTDIR}/${SAMPLE}_clinical_summary.tsv         (human-readable table)"
echo "    ${OUTDIR}/${SAMPLE}_high_impact.vcf.gz           (HIGH only)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate.vcf.gz         (rare MODERATE only)"
[ "$HAS_CLINVAR" -eq 1 ] && \
echo "    ${OUTDIR}/${SAMPLE}_clinvar_pathogenic.vcf.gz    (ClinVar P/LP only)"
{ [ "$HAS_CADD" -eq 1 ] || [ "$HAS_CADD_INDEL" -eq 1 ]; } && \
echo "    ${OUTDIR}/${SAMPLE}_cadd_high.vcf.gz             (CADD >= 20 non-coding)"
{ [ "$HAS_SPLICEAI" -eq 1 ] || [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; } && \
echo "    ${OUTDIR}/${SAMPLE}_spliceai_high.vcf.gz         (SpliceAI >= 0.2)"
[ "$REVEL_OR_AM" -eq 1 ] && \
echo "    ${OUTDIR}/${SAMPLE}_missense_deleterious.vcf.gz  (REVEL/AlphaMissense)"
echo "============================================"
echo ""
echo "Next: Review ${SAMPLE}_clinical_summary.tsv or load the VCF in IGV/gene.iobio"
