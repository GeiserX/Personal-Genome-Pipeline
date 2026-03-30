#!/usr/bin/env bash
# 25-prs.sh — Calculate Polygenic Risk Scores using plink2
# Usage: ./scripts/25-prs.sh <sample_name>
#
# Downloads PGS Catalog scoring files for common conditions and calculates
# polygenic risk scores from your VCF. PRS are NOT diagnostic — they estimate
# relative genetic predisposition compared to population averages.
#
# IMPORTANT: Raw PRS scores from a single sample are NOT directly interpretable.
# They only become meaningful when compared against a population distribution.
# Most GWAS-derived scores also have ancestry bias (European-centric).
# Treat these as exploratory, not clinical.
#
# Requires: VCF from step 3
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/prs"
SCORING_DIR="${GENOME_DIR}/prs_scores"
mkdir -p "$OUTDIR" "$SCORING_DIR"

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}"
  echo "  Run step 3 (DeepVariant) first."
  exit 1
fi

echo "============================================"
echo "  Step 25: Polygenic Risk Scores"
echo "  Tool: plink2"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${VCF}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Key PGS Catalog scoring files for common conditions
# These are well-validated, large-scale GWAS-derived scores
declare -A PGS_IDS
PGS_IDS=(
  ["coronary_artery_disease"]="PGS000018"
  ["type_2_diabetes"]="PGS000014"
  ["breast_cancer"]="PGS000004"
  ["prostate_cancer"]="PGS000662"
  ["atrial_fibrillation"]="PGS000016"
  ["alzheimers_disease"]="PGS000334"
  ["body_mass_index"]="PGS000027"
  ["schizophrenia"]="PGS000738"
  ["inflammatory_bowel_disease"]="PGS000020"
  ["colorectal_cancer"]="PGS000055"
)

# Download scoring files from PGS Catalog
echo "[1/3] Downloading PGS Catalog scoring files..."
for CONDITION in "${!PGS_IDS[@]}"; do
  PGS_ID="${PGS_IDS[$CONDITION]}"
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"

  if [ -f "$SCORE_FILE" ]; then
    echo "  [OK] ${CONDITION} (${PGS_ID}) — already downloaded"
  else
    echo "  Downloading ${CONDITION} (${PGS_ID})..."
    wget -q -O "$SCORE_FILE" \
      "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${PGS_ID}/ScoringFiles/Harmonized/${PGS_ID}_hmPOS_GRCh38.txt.gz" 2>/dev/null || \
    wget -q -O "$SCORE_FILE" \
      "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/${PGS_ID}/ScoringFiles/${PGS_ID}.txt.gz" 2>/dev/null || {
      echo "    WARNING: Could not download ${PGS_ID}. Skipping."
      rm -f "$SCORE_FILE"
    }
  fi
done

echo ""
echo "[2/3] Converting VCF to plink2 format..."

# Convert VCF to plink2 binary format for scoring
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  pgscatalog/plink2:2.00a5.10 \
  plink2 \
    --vcf "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --make-pgen \
    --out "/genome/${SAMPLE}/prs/${SAMPLE}" \
    --threads 4 \
    --memory 6000 \
    --set-all-var-ids '@:#' \
    --new-id-max-allele-len 100

echo ""
echo "[3/3] Calculating polygenic risk scores..."

# Process each scoring file
RESULTS_FILE="${OUTDIR}/${SAMPLE}_prs_summary.tsv"
echo -e "Condition\tPGS_ID\tScore\tVariants_Used\tVariants_Total" > "$RESULTS_FILE"

for CONDITION in "${!PGS_IDS[@]}"; do
  PGS_ID="${PGS_IDS[$CONDITION]}"
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"

  if [ ! -f "$SCORE_FILE" ]; then
    continue
  fi

  echo "  Scoring: ${CONDITION} (${PGS_ID})..."

  # Extract scoring columns from PGS Catalog format
  # PGS files have: rsID/chr_name, chr_position, effect_allele, effect_weight
  # Convert to plink2 --score format: variant_id, allele, weight
  FORMATTED="${OUTDIR}/${PGS_ID}_formatted.tsv"
  zcat "$SCORE_FILE" 2>/dev/null | grep -v "^#" | \
    awk -F'\t' 'NR==1 {
      for(i=1;i<=NF;i++) {
        if($i=="chr_name") chr_col=i;
        if($i=="chr_position") pos_col=i;
        if($i=="effect_allele") ea_col=i;
        if($i=="effect_weight") ew_col=i;
        if($i=="hm_chr") chr_col=i;
        if($i=="hm_pos") pos_col=i;
      }
      next
    }
    chr_col && pos_col && ea_col && ew_col {
      chr=$chr_col; pos=$pos_col; ea=$ea_col; ew=$ew_col;
      if(chr!="" && pos!="" && ea!="" && ew!="") {
        # Add chr prefix if missing to match GRCh38 VCF contig names
        if(chr !~ /^chr/) chr="chr"chr;
        printf "%s:%s\t%s\t%s\n", chr, pos, ea, ew;
      }
    }' > "$FORMATTED" 2>/dev/null || true

  TOTAL_VARS=$(wc -l < "$FORMATTED" 2>/dev/null || echo 0)

  if [ "$TOTAL_VARS" -eq 0 ]; then
    echo "    WARNING: Could not parse scoring file for ${PGS_ID}. Skipping."
    continue
  fi

  # Run plink2 --score
  docker run --rm --user root \
    --cpus 4 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    pgscatalog/plink2:2.00a5.10 \
    plink2 \
      --pfile "/genome/${SAMPLE}/prs/${SAMPLE}" \
      --score "/genome/${SAMPLE}/prs/${PGS_ID}_formatted.tsv" 1 2 3 \
        header-read \
        ignore-dup-ids \
        no-mean-imputation \
      --out "/genome/${SAMPLE}/prs/${PGS_ID}" \
      --threads 4 \
      --memory 3000 2>/dev/null || true

  # Extract score from .sscore file
  SSCORE="${OUTDIR}/${PGS_ID}.sscore"
  if [ -f "$SSCORE" ]; then
    SCORE=$(awk 'NR==2 {print $NF}' "$SSCORE" 2>/dev/null || echo "N/A")
    USED_VARS=$(awk 'NR==2 {print $(NF-1)}' "$SSCORE" 2>/dev/null || echo "N/A")
    echo -e "${CONDITION}\t${PGS_ID}\t${SCORE}\t${USED_VARS}\t${TOTAL_VARS}" >> "$RESULTS_FILE"
    echo "    Score: ${SCORE} (${USED_VARS}/${TOTAL_VARS} variants matched)"
  else
    echo "    WARNING: No score produced for ${PGS_ID}"
  fi
done

echo ""
echo "============================================"
echo "  Polygenic Risk Scores complete: ${SAMPLE}"
echo ""
echo "  Summary: ${RESULTS_FILE}"
cat "$RESULTS_FILE" | column -t 2>/dev/null || cat "$RESULTS_FILE"
echo ""
echo "============================================"
echo ""
echo "IMPORTANT: PRS are NOT diagnostic. They estimate relative genetic"
echo "predisposition. A high PRS does NOT mean you will develop the condition."
echo "Many factors (lifestyle, environment, other genes) are not captured."
echo ""
echo "These scores are most meaningful when compared against population"
echo "distributions, which requires a reference panel (not included)."
echo "See docs/25-prs.md for interpretation guidance."
