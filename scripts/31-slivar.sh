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

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"


SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
OUTDIR="${SAMPLE_DIR}/slivar"

# Input VCF: prefer vcfanno-enriched (step 30), fall back to VEP (step 13)
ANNOTATED_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_annotated.vcf.gz"
VEP_VCF_GZ="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf.gz"
VEP_VCF="${SAMPLE_DIR}/vep/${SAMPLE}_vep.vcf"

# Optional gene constraint data
CONSTRAINT_TSV="${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv"

# A derived file is used only when it is newer than what it was built from,
# so a re-run of step 13 is never hidden behind an older _vep.vcf.gz or an
# older step 30 output.
VEP_SRC=""
if [ -f "$VEP_VCF" ] && { [ ! -f "$VEP_VCF_GZ" ] || [ "$VEP_VCF" -nt "$VEP_VCF_GZ" ]; }; then
  VEP_SRC="$VEP_VCF"
elif [ -f "$VEP_VCF_GZ" ]; then
  VEP_SRC="$VEP_VCF_GZ"
fi
INPUT=""
HAS_VCFANNO=0
if [ -f "$ANNOTATED_VCF" ] && { [ -z "$VEP_SRC" ] || [ "$ANNOTATED_VCF" -nt "$VEP_SRC" ]; }; then
  INPUT="$ANNOTATED_VCF"
  HAS_VCFANNO=1
elif [ -n "$VEP_SRC" ]; then
  INPUT="$VEP_SRC"
  if [ -f "$ANNOTATED_VCF" ]; then
    echo "NOTICE: ${ANNOTATED_VCF} is older than ${VEP_SRC}; ignoring it. Run step 30 again for the score filters."
  fi
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
# The report reads the summary; a failed run must not leave an older one behind.
rm -f "${OUTDIR}/${SAMPLE}_slivar_summary.tsv"

# ── Step 1: Ensure VCF is bgzipped and indexed ────────────────────────
echo "[1/5] Preparing input VCF..."
if [ "$INPUT" = "$VEP_VCF" ]; then
  echo "  Compressing VEP VCF..."
  rm -f "${VEP_VCF_GZ:?}" "${VEP_VCF_GZ:?}.tbi"
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -c "bcftools view /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf -Oz \
      -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp && \
      mv /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz.tmp /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz && \
      bcftools index -f -t /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  INPUT="$VEP_VCF_GZ"
  CONTAINER_INPUT="/genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz"
  echo "  Done."
elif [ ! -f "${INPUT}.tbi" ] || [ "$INPUT" -nt "${INPUT}.tbi" ]; then
  echo "  Indexing VCF..."
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools index -f -t "$CONTAINER_INPUT"
  echo "  Done."
else
  echo "  VCF already compressed and indexed."
fi

# ── Step 2: Detect available annotation fields ────────────────────────
echo ""
echo "[2/5] Detecting annotation fields..."

VEP_FIELDS=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools +split-vep -l "$CONTAINER_INPUT" 2>/dev/null || echo "")

if [ -z "$VEP_FIELDS" ]; then
  echo "ERROR: No CSQ annotation found in VCF. Was VEP step 13 run correctly?"
  exit 1
fi

HAS_CLINVAR=0
printf '%s\n' "$VEP_FIELDS" | awk -F'\t' '$2 == "CLIN_SIG" {f = 1} END {exit !f}' && HAS_CLINVAR=1

# Rarity: VEP's MAX_AF (highest frequency in any 1000 Genomes or gnomAD
# exome/genome population), or gnomADe_AF and gnomADg_AF when MAX_AF is
# absent. Exome frequency alone would call a variant common in genomes but
# absent from exomes rare. The same rule as step 23.
FREQ_COLS=""
FREQ_EXPR=""
FREQ_NAME=""
if printf '%s\n' "$VEP_FIELDS" | awk -F'\t' '$2 == "MAX_AF" {f = 1} END {exit !f}'; then
  FREQ_COLS="MAX_AF:Float"
  FREQ_EXPR='(MAX_AF<0.01 || MAX_AF=".")'
  FREQ_NAME="MAX_AF"
else
  for f in gnomADe_AF gnomADg_AF; do
    printf '%s\n' "$VEP_FIELDS" | awk -F'\t' -v f="$f" '$2 == f {x = 1} END {exit !x}' || continue
    FREQ_COLS="${FREQ_COLS:+${FREQ_COLS},}${f}:Float"
    FREQ_EXPR="${FREQ_EXPR:+${FREQ_EXPR} && }(${f}<0.01 || ${f}=\".\")"
    FREQ_NAME="${FREQ_NAME:+${FREQ_NAME} and }${f}"
  done
fi

# Check for vcfanno INFO fields
HAS_CADD=0
HAS_CADD_INDEL=0
HAS_REVEL=0
HAS_AM=0
HAS_SPLICEAI=0
HAS_SPLICEAI_INDEL=0
if [ "$HAS_VCFANNO" -eq 1 ]; then
  INFO_HEADER=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -h "$CONTAINER_INPUT" 2>/dev/null | grep '^##INFO' || echo "")
  echo "$INFO_HEADER" | grep -q 'ID=CADD_PHRED,' && HAS_CADD=1
  echo "$INFO_HEADER" | grep -q 'ID=CADD_PHRED_indel,' && HAS_CADD_INDEL=1
  echo "$INFO_HEADER" | grep -q 'ID=REVEL' && HAS_REVEL=1
  echo "$INFO_HEADER" | grep -q 'ID=AM_class' && HAS_AM=1
  echo "$INFO_HEADER" | grep -q 'ID=SpliceAI,' && HAS_SPLICEAI=1
  echo "$INFO_HEADER" | grep -q 'ID=SpliceAI_indel,' && HAS_SPLICEAI_INDEL=1
fi

echo "  VEP CSQ fields: frequency=${FREQ_NAME:-none}, ClinVar=$([ "$HAS_CLINVAR" -eq 1 ] && echo 'yes' || echo 'no')"
if [ -z "$FREQ_EXPR" ]; then
  echo "  NOTICE: no MAX_AF, gnomADe_AF or gnomADg_AF in the VEP output: the rare tiers are not filtered by frequency."
fi
echo "  vcfanno INFO fields: CADD=$([ "$HAS_CADD" -eq 1 ] && echo 'yes' || echo 'no'), REVEL=$([ "$HAS_REVEL" -eq 1 ] && echo 'yes' || echo 'no'), AlphaMissense=$([ "$HAS_AM" -eq 1 ] && echo 'yes' || echo 'no'), SpliceAI=$([ "$HAS_SPLICEAI" -eq 1 ] && echo 'yes' || echo 'no')"
echo ""

# ── Step 3: Variant prioritization with bcftools ──────────────────────
# Uses bcftools +split-vep for CSQ fields and bcftools view -i for INFO fields.
# Three filter tiers: rare_high, rare_moderate_deleterious, clinvar_pathogenic.
echo "[3/5] Prioritizing variants..."

# --- Filter 1: rare_high ---
# PASS + HIGH VEP impact + rare (see FREQ_EXPR above)
RARE_COLS="IMPACT${FREQ_COLS:+,${FREQ_COLS}}"
RARE_AND="${FREQ_EXPR:+ && ${FREQ_EXPR}}"
echo "  [a] rare_high: PASS + HIGH impact + ${FREQ_NAME:-no frequency filter}..."
run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
    bcftools +split-vep - -c '${RARE_COLS}' -s worst \
      -i 'IMPACT=\"HIGH\"${RARE_AND}' \
      -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz && \
    bcftools index -f -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz"

RARE_HIGH_COUNT=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_high.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${RARE_HIGH_COUNT} variants"

# --- Filter 2: rare_moderate_deleterious ---
# PASS + MODERATE impact + gnomAD AF < 0.01 + at least one deleterious predictor
echo "  [b] rare_moderate_deleterious: PASS + MODERATE + rare + deleterious predictors..."

# Build the filter expression depending on available annotations
MODERATE_FILTER="IMPACT=\"MODERATE\"${RARE_AND}"
VEP_COLUMNS="$RARE_COLS"

# If vcfanno annotations are available, add predictor thresholds as a second pass
if [ "$HAS_VCFANNO" -eq 1 ]; then
  # First pass: extract rare MODERATE via split-vep, then second pass: filter on INFO fields

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

    run_in --cpus 2 --memory 4g \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
        bcftools +split-vep - -c '${VEP_COLUMNS}' -s worst \
          -i '${MODERATE_FILTER}' | \
        bcftools view -i '${PREDICTOR_EXPR}' \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  else
    # No predictors available despite vcfanno — fall back to all rare MODERATE
    run_in --cpus 2 --memory 4g \
      "${BCFTOOLS_IMAGE}" \
      bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
        bcftools +split-vep - -c '${VEP_COLUMNS}' -s worst \
          -i '${MODERATE_FILTER}' \
          -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
        bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
  fi
else
  # No vcfanno — include all rare MODERATE (same as step 23 behavior)
  echo "    WARNING: No vcfanno annotations — including all rare MODERATE variants."
  echo "    Run step 30 (vcfanno) for CADD/REVEL/AlphaMissense/SpliceAI filtering."

  run_in --cpus 2 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c '${VEP_COLUMNS}' -s worst \
        -i '${MODERATE_FILTER}' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz && \
      bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz"
fi

MODERATE_DEL_COUNT=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_rare_moderate_del.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "      Found: ${MODERATE_DEL_COUNT} variants"

# --- Filter 3: clinvar_pathogenic ---
CLINVAR_COUNT=0
CLINVAR_FILE=""
if [ "$HAS_CLINVAR" -eq 1 ]; then
  echo "  [c] clinvar_pathogenic: ClinVar pathogenic/likely_pathogenic..."
  run_in --cpus 2 --memory 4g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools view -f PASS ${CONTAINER_INPUT} | \
      bcftools +split-vep - -c CLIN_SIG \
        -i 'CLIN_SIG~\"pathogenic\" && CLIN_SIG!~\"conflicting\"' \
        -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_clinvar_path.vcf.gz"

  CLINVAR_COUNT=$(run_in \
    --cpus 2 --memory 2g \
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

run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools concat -a -D \
    ${MERGE_FILES} | \
    bcftools sort -Oz -o /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz && \
    bcftools index -t /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz"

PRIORITIZED_COUNT=$(run_in \
  --cpus 2 --memory 2g \
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

COMPHET_LOG="${OUTDIR}/${SAMPLE}_compound_hets.log"
# Remove stale outputs from a previous run so a failure cannot leave old results behind
rm -f "$COMPHET_VCF" "$COMPHET_TSV"

# A failure here (wrong image, slivar crash, broken input) stops the step with
# slivar's own error. It must never be reported as "no candidates found".
if ! run_in --cpus 2 --memory 4g \
  "${SLIVAR_IMAGE}" \
  slivar compound-hets \
    --allow-non-trios \
    --vcf "/genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz" \
    --ped "/genome/${SAMPLE}/slivar/${SAMPLE}.ped" \
  2>"$COMPHET_LOG" | \
  run_in -i    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -Oz -o "/genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz"; then
  echo "ERROR: slivar compound-hets failed (image: ${SLIVAR_IMAGE})." >&2
  echo "  Output of the failed command (also in ${COMPHET_LOG}):" >&2
  sed 's/^/    /' "$COMPHET_LOG" >&2 || true
  rm -f "$COMPHET_VCF"
  exit 1
fi

# Count pairs from slivar_comphet INFO field (not VCF record count).
# slivar compound-hets outputs one record per unique VARIANT, with the
# slivar_comphet INFO field listing all partner variants. Format:
#   sample/GENE/PAIR_ID/chrom/pos/ref/alt  (comma-separated when multiple)
# A gene with N variants produces C(N,2) pairs but only N VCF records.
COMPHET_PAIRS=0
COMPHET_GENES=0
COMPHET_RECORDS=$(run_in \
  --cpus 2 --memory 2g \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -H "/genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz" | wc -l | tr -d ' ')
if [ "$COMPHET_RECORDS" -gt 0 ]; then
  COMPHET_STATS=$(run_in \
    --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools query -f '%INFO/slivar_comphet\n' \
      /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz \
    | tr ',' '\n' \
    | awk -F'/' '{pairs[\$3]=1; genes[\$2]=1} END{print length(pairs), length(genes)}'")
  COMPHET_PAIRS=$(echo "$COMPHET_STATS" | awk '{print $1}')
  COMPHET_GENES=$(echo "$COMPHET_STATS" | awk '{print $2}')

  # Export to TSV sorted by gene for human review
  run_in --cpus 2 --memory 2g \
    "${BCFTOOLS_IMAGE}" \
    bash -o pipefail -c "bcftools +split-vep \
        /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets.vcf.gz \
        -f '%SYMBOL\t%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%Consequence[\t%GT]\n' \
        -s worst -d \
      | sort -t\$'\t' -k1,1 -k2,2V -k3,3n \
      > /genome/${SAMPLE}/slivar/${SAMPLE}_compound_hets_raw.tsv"
  {
    printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\n'
    cat "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
  } > "${COMPHET_TSV}.tmp"
  mv "${COMPHET_TSV}.tmp" "$COMPHET_TSV"
  rm -f "${OUTDIR}/${SAMPLE}_compound_hets_raw.tsv"
  echo "  Found: ${COMPHET_PAIRS} compound het candidate pairs across ${COMPHET_GENES} genes"
else
  echo "  slivar ran and returned no records: no compound heterozygote candidates found."
fi

if [ ! -s "$COMPHET_TSV" ]; then
  printf 'GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\n' > "$COMPHET_TSV"
fi

# ── Step 5: Generate summary TSV with gene constraint ─────────────────
echo ""
echo "[5/5] Generating summary with gene constraint annotations..."

SUMMARY_TSV="${OUTDIR}/${SAMPLE}_slivar_summary.tsv"

# Extract variant info from prioritized VCF into a TSV
if ! run_in --cpus 2 --memory 4g \
  "${BCFTOOLS_IMAGE}" \
  bash -o pipefail -c "bcftools +split-vep \
    /genome/${SAMPLE}/slivar/${SAMPLE}_prioritized.vcf.gz \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%IMPACT\t%SYMBOL\t%Consequence\t%Existing_variation[\t%GT]\n' \
    -s worst -d \
  > /genome/${SAMPLE}/slivar/${SAMPLE}_variants_raw.tsv"; then
  echo "ERROR: Summary extraction from ${OUTDIR}/${SAMPLE}_prioritized.vcf.gz failed." >&2
  rm -f "${OUTDIR}/${SAMPLE}_variants_raw.tsv"
  exit 1
fi

# Add header and optional gene constraint columns (bin/constraint_join.awk,
# the loader step 23 and the Nextflow slivar module use: canonical rows,
# mis.z_score, the Ensembl row over the RefSeq one; it exits non-zero when rows
# carry genes and not one matches the file).
{
  printf 'CHROM\tPOS\tREF\tALT\tIMPACT\tSYMBOL\tConsequence\tExisting_variation\tGT\n'
  cat "${OUTDIR}/${SAMPLE}_variants_raw.tsv"
} > "${SUMMARY_TSV}.raw"
if [ -f "$CONSTRAINT_TSV" ]; then
  echo "  Joining with gnomAD gene constraint metrics..."
  awk -f "${PGP_ROOT}/bin/constraint_join.awk" gene_col=SYMBOL constrained=1 \
    "$CONSTRAINT_TSV" "${SUMMARY_TSV}.raw" > "${SUMMARY_TSV}.tmp"
  rm -f "${SUMMARY_TSV:?}.raw"
else
  echo "  Gene constraint file not found (optional): ${CONSTRAINT_TSV}"
  echo "  Generating summary without constraint annotations."
  mv "${SUMMARY_TSV}.raw" "${SUMMARY_TSV}.tmp"
fi
mv "${SUMMARY_TSV}.tmp" "$SUMMARY_TSV"

# Clean up intermediate file
rm -f "${OUTDIR}/${SAMPLE}_variants_raw.tsv"

# ── Summary ───────────────────────────────────────────────────────────
echo ""
echo "============================================"
echo "  Step 31 complete: ${SAMPLE}"
echo ""
echo "  Filter results:"
echo "    rare_high (HIGH, rare):            ${RARE_HIGH_COUNT}"
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
