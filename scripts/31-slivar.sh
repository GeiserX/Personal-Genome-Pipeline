#!/usr/bin/env bash
# 31-slivar.sh — Variant prioritization and compound heterozygote detection
# Usage: ./scripts/31-slivar.sh <sample_name>
#
# Prioritizes clinically interesting variants using tiered filters and detects
# compound heterozygote candidates using slivar. Optionally annotates results
# with gnomAD gene constraint metrics (LOEUF, pLI).
#
# Input: vcfanno-enriched VCF (step 30) or VEP-annotated VCF (step 13)
# Output: prioritized VCF, compound het TSV, summary TSV
# Requires: VEP-annotated VCF. Step 30 (vcfanno) recommended for full filtering.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

# Validate sample name to prevent shell injection in bash -c / python3 -c strings
if [[ "$SAMPLE" =~ [^a-zA-Z0-9._-] ]]; then
  echo "ERROR: Sample name contains invalid characters. Use only a-z, A-Z, 0-9, ., _, -" >&2
  exit 1
fi

SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${SAMPLE_DIR}/slivar"

# Input VCF: prefer vcfanno-enriched (step 30), fall back to VEP (step 13)
ANNOTATED_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_annotated.vcf.gz"
VEP_VCF_GZ="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf.gz"
VEP_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf"

# Optional gene constraint data
CONSTRAINT_TSV="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"

INPUT=""
HAS_VCFANNO=0
if [ -f "$ANNOTATED_VCF" ]; then
  INPUT="$ANNOTATED_VCF"
  HAS_VCFANNO=1
elif [ -f "$VEP_VCF_GZ" ]; then
  INPUT="$VEP_VCF_GZ"
elif [ -f "$VEP_VCF" ]; then
  INPUT="$VEP_VCF"
else
  echo "ERROR: No annotated VCF found. Run step 13 (VEP) first."
  echo "  Expected: ${ANNOTATED_VCF} or ${VEP_VCF_GZ} or ${VEP_VCF}"
  exit 1
fi

CONTAINER_INPUT="/genome/${SAMPLE}/vep/$(basename "$INPUT")"

echo "============================================"
echo "  Step 31: Variant Prioritization (slivar)"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${INPUT}"
echo "  vcfanno annotations: $([ "$HAS_VCFANNO" -eq 1 ] && echo 'yes (full filtering)' || echo 'no (VEP-only mode)')"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

mkdir -p "$OUTDIR"

# ── Step 1: Ensure VCF is bgzipped and indexed ────────────────────────
echo "[1/5] Preparing input VCF..."
if [[ "$INPUT" == *.vcf ]] && [ ! -f "$VEP_VCF_GZ" ]; then
  echo "  Compressing VEP VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf -Oz \
      -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  CONTAINER_INPUT="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  echo "  Done."
elif [[ "$INPUT" == *.vcf.gz ]] && [ ! -f "${INPUT}.tbi" ]; then
  echo "  Indexing VCF..."
  docker run --rm --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools index -t "$CONTAINER_INPUT"
  echo "  Done."
else
  echo "  VCF already compressed and indexed."
fi

# ── Step 2: Detect available annotation fields ────────────────────────
echo ""
echo "[2/5] Detecting annotation fields..."

VEP_FIELDS=$(docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools +split-vep -l "$CONTAINER_INPUT" 2>/dev/null || echo "")

if [ -z "$VEP_FIELDS" ]; then
  echo "ERROR: No CSQ annotation found in VCF. Was VEP step 13 run correctly?"
  exit 1
fi

HAS_GNOMAD=0
HAS_CLINVAR=0
echo "$VEP_FIELDS" | grep -q 'gnomADe_AF' && HAS_GNOMAD=1
echo "$VEP_FIELDS" | grep -q 'CLIN_SIG' && HAS_CLINVAR=1

# Check for vcfanno INFO fields
HAS_CADD=0
HAS_CADD_INDEL=0
HAS_REVEL=0
HAS_AM=0
HAS_SPLICEAI=0
HAS_SPLICEAI_INDEL=0
if [ "$HAS_VCFANNO" -eq 1 ]; then
  INFO_HEADER=$(docker run --rm \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -h "$CONTAINER_INPUT" 2>/dev/null | grep '^##INFO' || echo "")
  echo "$INFO_HEADER" | grep -q 'ID=CADD_PHRED,' && HAS_CADD=1
  echo "$INFO_HEADER" | grep -q 'ID=CADD_PHRED_indel,' && HAS_CADD_INDEL=1
  echo "$INFO_HEADER" | grep -q 'ID=REVEL' && HAS_REVEL=1
  echo "$INFO_HEADER" | grep -q 'ID=AM_class' && HAS_AM=1
  echo "$INFO_HEADER" | grep -q 'ID=SpliceAI,' && HAS_SPLICEAI=1
  echo "$INFO_HEADER" | grep -q 'ID=SpliceAI_indel,' && HAS_SPLICEAI_INDEL=1
fi

echo "  VEP CSQ fields: gnomAD=$([ "$HAS_GNOMAD" -eq 1 ] && echo 'yes' || echo 'no'), ClinVar=$([ "$HAS_CLINVAR" -eq 1 ] && echo 'yes' || echo 'no')"
echo "  vcfanno INFO fields: CADD=$([ "$HAS_CADD" -eq 1 ] && echo 'yes' || echo 'no'), REVEL=$([ "$HAS_REVEL" -eq 1 ] && echo 'yes' || echo 'no'), AlphaMissense=$([ "$HAS_AM" -eq 1 ] && echo 'yes' || echo 'no'), SpliceAI=$([ "$HAS_SPLICEAI" -eq 1 ] && echo 'yes' || echo 'no')"
echo ""

# ── Step 3: Variant prioritization with bcftools ──────────────────────
# Uses bcftools +split-vep for CSQ fields and bcftools view -i for INFO fields.
# Three filter tiers: rare_high, rare_moderate_deleterious, clinvar_pathogenic.
echo "[3/5] Prioritizing variants..."

# --- Filter 1: rare_high ---
# PASS + HIGH VEP impact + gnomAD AF < 0.01 (or missing)
echo "  [a] rare_high: PASS + HIGH impact + gnomAD AF < 1%..."
if [ "$HAS_GNOMAD" -eq 1 ]; then
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT,gnomADe_AF -s worst \
        -i 'IMPACT=\"HIGH\" && (gnomADe_AF<0.01 || gnomADe_AF=\".\")' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz"
else
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c IMPACT -s worst \
        -i 'IMPACT=\"HIGH\"' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz"
fi

RARE_HIGH_COUNT=$(docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${RARE_HIGH_COUNT} variants"

# --- Filter 2: rare_moderate_deleterious ---
# PASS + MODERATE impact + gnomAD AF < 0.01 + at least one deleterious predictor
echo "  [b] rare_moderate_deleterious: PASS + MODERATE + rare + deleterious predictors..."

# Build the filter expression depending on available annotations
MODERATE_FILTER='IMPACT="MODERATE" && (gnomADe_AF<0.01 || gnomADe_AF=".")'
if [ "$HAS_GNOMAD" -eq 0 ]; then
  MODERATE_FILTER='IMPACT="MODERATE"'
fi

# If vcfanno annotations are available, add predictor thresholds as a second pass
if [ "$HAS_VCFANNO" -eq 1 ]; then
  # First pass: extract rare MODERATE via split-vep, then second pass: filter on INFO fields
  VEP_COLUMNS="IMPACT"
  [ "$HAS_GNOMAD" -eq 1 ] && VEP_COLUMNS="IMPACT,gnomADe_AF"

  # Build INFO-level predictor filter — only reference tags that exist in the header
  PREDICTOR_PARTS=()
  if [ "$HAS_CADD" -eq 1 ] && [ "$HAS_CADD_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED>=20 || INFO/CADD_PHRED_indel>=20')
  elif [ "$HAS_CADD" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED>=20')
  elif [ "$HAS_CADD_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/CADD_PHRED_indel>=20')
  fi
  [ "$HAS_REVEL" -eq 1 ] && PREDICTOR_PARTS+=('INFO/REVEL>=0.5')
  [ "$HAS_AM" -eq 1 ] && PREDICTOR_PARTS+=('INFO/AM_class="likely_pathogenic"')
  # Note: SpliceAI is a pipe-delimited string, not a numeric field.
  # bcftools cannot numerically compare sub-fields, so we use presence check only.
  # Variants with SpliceAI annotations are included; threshold filtering happens in step 23.
  if [ "$HAS_SPLICEAI" -eq 1 ] && [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI!="." || INFO/SpliceAI_indel!="."')
  elif [ "$HAS_SPLICEAI" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI!="."')
  elif [ "$HAS_SPLICEAI_INDEL" -eq 1 ]; then
    PREDICTOR_PARTS+=('INFO/SpliceAI_indel!="."')
  fi

  if [ ${#PREDICTOR_PARTS[@]} -gt 0 ]; then
    PREDICTOR_EXPR=$(printf ' || %s' "${PREDICTOR_PARTS[@]}")
    PREDICTOR_EXPR="${PREDICTOR_EXPR:4}"  # strip leading ' || '

    docker run --rm --user root \
      --cpus 2 --memory 4g \
      -v "${GENOME_DIR}:/genome" \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
        bcftools +split-vep - -c ${VEP_COLUMNS} -s worst \
          -i '${MODERATE_FILTER}' | \
        bcftools view -i '${PREDICTOR_EXPR}' \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  else
    # No predictors available despite vcfanno — fall back to all rare MODERATE
    docker run --rm --user root \
      --cpus 2 --memory 4g \
      -v "${GENOME_DIR}:/genome" \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
        bcftools +split-vep - -c ${VEP_COLUMNS} -s worst \
          -i '${MODERATE_FILTER}' \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  fi
else
  # No vcfanno — include all rare MODERATE (same as step 23 behavior)
  VEP_COLUMNS="IMPACT"
  [ "$HAS_GNOMAD" -eq 1 ] && VEP_COLUMNS="IMPACT,gnomADe_AF"
  echo "    WARNING: No vcfanno annotations — including all rare MODERATE variants."
  echo "    Run step 30 (vcfanno) for CADD/REVEL/AlphaMissense/SpliceAI filtering."

  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c ${VEP_COLUMNS} -s worst \
        -i '${MODERATE_FILTER}' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
fi

MODERATE_DEL_COUNT=$(docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${MODERATE_DEL_COUNT} variants"

# --- Filter 3: clinvar_pathogenic ---
CLINVAR_COUNT=0
CLINVAR_FILE=""
if [ "$HAS_CLINVAR" -eq 1 ]; then
  echo "  [c] clinvar_pathogenic: ClinVar pathogenic/likely_pathogenic..."
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"

  CLINVAR_COUNT=$(docker run --rm \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz" 2>/dev/null | wc -l || echo 0)
  echo "      Found: ${CLINVAR_COUNT} variants"
  CLINVAR_FILE="/genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"
else
  echo "  [c] clinvar_pathogenic: skipped (CLIN_SIG not in VEP annotations)"
fi

# Merge all tiers into a single prioritized VCF
echo ""
echo "  Merging filter tiers into prioritized VCF..."
MERGE_FILES="/genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
[ -n "$CLINVAR_FILE" ] && MERGE_FILES="${MERGE_FILES} ${CLINVAR_FILE}"

docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools concat -a -D \
    ${MERGE_FILES} | \
    bcftools sort -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz"

PRIORITIZED_COUNT=$(docker run --rm \
  --cpus 2 --memory 2g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Total prioritized variants: ${PRIORITIZED_COUNT}"

# ── Step 4: Compound heterozygote detection with slivar ───────────────
echo ""
echo "[4/5] Detecting compound heterozygote candidates..."

# slivar compound-hets requires:
#   --vcf: input VCF with CSQ annotations
#   --ped: PED file describing sample relationships
#   --allow-non-trios: required for singleton/duo samples (no trio structure)
# It groups heterozygous variants by gene (from CSQ) and outputs VCF to stdout
# with pairs of variants per gene.
# For single-sample unphased data, these are CANDIDATES only.
COMPHET_VCF="${OUTDIR}/${SAMPLE}_compound_hets.vcf.gz"
COMPHET_TSV="${OUTDIR}/${SAMPLE}_compound_hets.tsv"
COMPHET_PED="${OUTDIR}/${SAMPLE}.ped"

# Generate minimal PED file for single sample (no parents, unaffected)
# Format: family_id sample_id father mother sex phenotype
echo -e "${SAMPLE}\t${SAMPLE}\t0\t0\t0\t-9" > "$COMPHET_PED"

if docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "${SLIVAR_IMAGE}" \
  slivar compound-hets \
    --allow-non-trios \
    --vcf "/genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" \
    --ped "/genome/${SAMPLE}/slivar/${SAMPLE}.ped" \
  2>"${OUTDIR}/${SAMPLE}_compound_hets.log" | \
  docker run --rm -i --user root \
    --cpus 2 --memory 2g \
    -v "${GENOME_DIR}:/genome" \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -Oz -o "/genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz"; then

  # Remove any stale TSV from a previous run before processing
  rm -f "$COMPHET_TSV"

  # Count pairs from slivar_comphet INFO field (not VCF record count).
  # slivar compound-hets outputs one record per unique VARIANT, with the
  # slivar_comphet INFO field listing all partner variants. Format:
  #   sample/GENE/PAIR_ID/chrom/pos/ref/alt  (comma-separated when multiple)
  # A gene with N variants produces C(N,2) pairs but only N VCF records.
  COMPHET_PAIRS=0
  COMPHET_GENES=0
  if [ -s "$COMPHET_VCF" ]; then
    COMPHET_STATS=$(docker run --rm \
      --cpus 2 --memory 2g \
      -v "${GENOME_DIR}:/genome" \
      "${BCFTOOLS_IMAGE}" \
      bash -c "bcftools query -f '%INFO/slivar_comphet\n' \
        /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz 2>/dev/null \
      | tr ',' '\n' \
      | awk -F'/' '{pairs[\$3]=1; genes[\$2]=1} END{print length(pairs), length(genes)}'" \
      2>/dev/null || echo "0 0")
    COMPHET_PAIRS=$(echo "$COMPHET_STATS" | awk '{print $1}')
    COMPHET_GENES=$(echo "$COMPHET_STATS" | awk '{print $2}')

    # Export to TSV sorted by gene for human review
    if docker run --rm --user root \
      --cpus 2 --memory 2g \
      -v "${GENOME_DIR}:/genome" \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools +split-vep \
          /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz \
          -f '%SYMBOL\t%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%Consequence[\t%GT]\n' \
          -s worst -d \
        | sort -t\$'\t' -k1,1 -k2,2V -k3,3n \
        > /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets_raw.tsv"; then
      {
        printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\n'
        cat "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
      } > "$COMPHET_TSV"
      rm -f "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
    else
      rm -f "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
      echo "  WARNING: TSV conversion failed; VCF is available for manual inspection."
    fi
  fi
  echo "  Found: ${COMPHET_PAIRS} compound het candidate pairs across ${COMPHET_GENES} genes"
else
  echo "  WARNING: slivar compound-hets failed. See ${OUTDIR}/${SAMPLE}_compound_hets.log"
  echo "  No compound heterozygote candidates found."
  COMPHET_PAIRS=0
fi

if [ ! -s "$COMPHET_TSV" ]; then
  printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\n' > "$COMPHET_TSV"
fi

# ── Step 5: Generate summary TSV with gene constraint ─────────────────
echo ""
echo "[5/5] Generating summary with gene constraint annotations..."

SUMMARY_TSV="${OUTDIR}/${SAMPLE}_slivar_summary.tsv"

# Extract variant info from prioritized VCF into a TSV
if ! docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep \
    /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%Consequence\t%Existing_variation[\t%GT]\n' \
    -s worst -d \
  > /genome/${SAMPLE}/slivar/${SAMPLE}_variants_raw.tsv"; then
  echo "  WARNING: Summary extraction failed. Prioritized VCF is still available."
fi

# Add header and optional gene constraint columns
if [ -f "$CONSTRAINT_TSV" ]; then
  echo "  Joining with gnomAD gene constraint metrics..."
  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${PYTHON_IMAGE}" \
    python3 -c "
import csv, sys

# Load gene constraint data
constraint = {}
try:
    with open('/genome/annotations/gnomad_v4.1_constraint.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('gene', row.get('gene_symbol', row.get('symbol', '')))
            if not gene:
                continue
            loeuf = row.get('oe_lof_upper', row.get('lof.oe_ci.upper', '.'))
            pli = row.get('pLI', row.get('lof.pLI', '.'))
            mis_z = row.get('mis_z', row.get('missense.z_score', '.'))
            constraint[gene] = (loeuf, pli, mis_z)
except Exception as e:
    print(f'WARNING: Could not load constraint file: {e}', file=sys.stderr)

# Process variants
header = 'CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT\tLOEUF\tpLI\tmis_z\tCONSTRAINED'
print(header)

try:
    with open('/genome/${SAMPLE}/slivar/${SAMPLE}_variants_raw.tsv') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            gene = fields[5]
            loeuf, pli, mis_z = constraint.get(gene, ('.', '.', '.'))
            # Mark as constrained if LOEUF < 0.35 or pLI > 0.9
            constrained = 'NO'
            try:
                if loeuf != '.' and float(loeuf) < 0.35:
                    constrained = 'YES'
                elif pli != '.' and float(pli) > 0.9:
                    constrained = 'YES'
            except ValueError:
                pass
            print(f'{line}\t{loeuf}\t{pli}\t{mis_z}\t{constrained}')
except FileNotFoundError:
    pass
" > "$SUMMARY_TSV" 2>/dev/null || true
else
  echo "  Gene constraint file not found (optional): ${CONSTRAINT_TSV}"
  echo "  Generating summary without constraint annotations."
  {
    echo -e "CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT"
    cat "${OUTDIR}/${SAMPLE}_variants_raw.tsv" 2>/dev/null || true
  } > "$SUMMARY_TSV"
fi

# Clean up intermediate file
rm -f "${OUTDIR}/${SAMPLE}_variants_raw.tsv"

# ── Summary ───────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 31 complete: ${SAMPLE}"
echo ""
echo "  Filter results:"
echo "    rare_high (HIGH + AF<1%):          ${RARE_HIGH_COUNT}"
echo "    rare_moderate_deleterious:          ${MODERATE_DEL_COUNT}"
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    clinvar_pathogenic:                 ${CLINVAR_COUNT}"
fi
echo "    ────────────────────────────────────"
echo "    Total prioritized (deduplicated):   ${PRIORITIZED_COUNT}"
echo ""
echo "  Compound heterozygote candidates:     ${COMPHET_PAIRS} pairs (${COMPHET_GENES} genes)"
echo ""

# Highlight constrained genes if constraint data was used
if [ -f "$CONSTRAINT_TSV" ] && [ -f "$SUMMARY_TSV" ]; then
  CONSTRAINED_COUNT=$(tail -n +2 "$SUMMARY_TSV" 2>/dev/null | awk -F'\t' '$NF=="YES"' | wc -l | tr -d ' ')
  if [ "$CONSTRAINED_COUNT" -gt 0 ]; then
    echo "  Variants in constrained genes:        ${CONSTRAINED_COUNT}"
    echo "  (LOEUF < 0.35 or pLI > 0.9 — loss-of-function intolerant)"
    echo ""
    echo "  Top constrained gene hits:"
    tail -n +2 "$SUMMARY_TSV" | awk -F'\t' '$NF=="YES" {print "    "$6" ("$5") "$7}' | sort -u | head -10
    echo ""
  fi
fi

echo "  Output files:"
echo "    ${OUTDIR}/${SAMPLE}_prioritized.vcf.gz         (all prioritized variants)"
echo "    ${OUTDIR}/${SAMPLE}_compound_hets.vcf.gz       (compound het VCF)"
echo "    ${OUTDIR}/${SAMPLE}_compound_hets.tsv           (compound het candidates)"
echo "    ${OUTDIR}/${SAMPLE}_slivar_summary.tsv          (summary + gene constraint)"
echo "    ${OUTDIR}/${SAMPLE}_rare_high.vcf.gz            (HIGH impact tier)"
echo "    ${OUTDIR}/${SAMPLE}_rare_moderate_del.vcf.gz    (MODERATE + deleterious)"
if [ "$HAS_CLINVAR" -eq 1 ]; then
echo "    ${OUTDIR}/${SAMPLE}_clinvar_path.vcf.gz         (ClinVar P/LP)"
fi
echo "============================================"
echo ""
echo "NOTE: Compound het candidates from single-sample unphased data are"
echo "  not confirmed. Phased data (trio or read-backed) is needed to"
echo "  distinguish true compound hets from variants on the same haplotype."
echo ""
echo "Next: Review ${SAMPLE}_slivar_summary.tsv or load prioritized VCF in IGV/gene.iobio"
