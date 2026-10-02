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
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

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

# PGS Catalog scores, one per line: "<PGS ID>|<trait_reported>".
# The label is the trait exactly as https://www.pgscatalog.org/rest/score/<PGS ID>
# reports it, so a wrong ID cannot be printed under the wrong disease.
PGS_SCORES=(
  "PGS000018|Coronary artery disease"
  "PGS000014|Type 2 diabetes (T2D)"
  "PGS000004|Breast cancer"
  "PGS000662|Prostate cancer"
  "PGS000016|Atrial fibrillation"
  "PGS000334|Late-onset Alzheimer’s disease"
  "PGS000027|Body mass index (BMI)"
  "PGS000017|Inflammatory bowel disease"
  "PGS000055|Colorectal cancer"
)

# Only GRCh38-harmonised scoring files are used. The author-reported files are
# often GRCh37 or rsID-only, and scoring them against a GRCh38 VCF gives a
# number that looks valid but is not, so there is no fallback to them.
PGS_BASE_URL=${PGS_BASE_URL:-https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores}

# Prints the genome build recorded in a scoring file's #HmPOS_build header (empty if none).
hm_build() {
  { gzip -cd "$1" 2>/dev/null || true; } | awk -F= '/^#HmPOS_build=/ {b=$2} !/^#/ {exit} END {print b}'
}

# Download scoring files from PGS Catalog
echo "[1/3] Downloading PGS Catalog scoring files..."
for ENTRY in "${PGS_SCORES[@]}"; do
  PGS_ID="${ENTRY%%|*}"
  CONDITION="${ENTRY#*|}"
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"

  # A cached file from an older version of this script may be the author-reported
  # build; drop it and fetch the harmonised one.
  if [ -f "$SCORE_FILE" ] && [ "$(hm_build "$SCORE_FILE")" != "GRCh38" ]; then
    echo "  Cached ${PGS_ID} is not a GRCh38-harmonised file; downloading it again."
    rm -f "$SCORE_FILE"
  fi

  if [ -f "$SCORE_FILE" ]; then
    echo "  [OK] ${CONDITION} (${PGS_ID}) — already downloaded"
    continue
  fi

  echo "  Downloading ${CONDITION} (${PGS_ID})..."
  URL="${PGS_BASE_URL}/${PGS_ID}/ScoringFiles/Harmonized/${PGS_ID}_hmPOS_GRCh38.txt.gz"
  if ! wget -q -O "${SCORE_FILE}.part" "$URL"; then
    rm -f "${SCORE_FILE}.part"
    echo "ERROR: Could not download the GRCh38-harmonised scoring file for ${PGS_ID}:" >&2
    echo "  ${URL}" >&2
    exit 1
  fi
  BUILD=$(hm_build "${SCORE_FILE}.part")
  if [ "$BUILD" != "GRCh38" ]; then
    rm -f "${SCORE_FILE}.part"
    echo "ERROR: ${PGS_ID} scoring file has #HmPOS_build='${BUILD}', expected GRCh38. Refusing to score it." >&2
    exit 1
  fi
  mv "${SCORE_FILE}.part" "$SCORE_FILE"
done

echo ""
echo "[2/3] Converting VCF to plink2 format..."

# Convert VCF to plink2 binary format for scoring
run_in  --cpus 4 --memory 8g \
  "${PLINK2_IMAGE}" \
  plink2 \
    --vcf "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --make-pgen \
    --out "/genome/${SAMPLE}/prs/${SAMPLE}" \
    --threads 4 \
    --memory 6000 \
    --set-all-var-ids '@:#' \
    --new-id-max-allele-len 100 \
    --chr 1-22 \
    --allow-extra-chr \
    --output-chr chrM

echo ""
echo "[3/3] Calculating polygenic risk scores..."

HOMREF_NOTE="hom-ref sites are absent from this VCF, so the score is biased; not comparable to published distributions"

# Process each scoring file
RESULTS_FILE="${OUTDIR}/${SAMPLE}_prs_summary.tsv"
echo -e "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total" > "$RESULTS_FILE"

for ENTRY in "${PGS_SCORES[@]}"; do
  PGS_ID="${ENTRY%%|*}"
  CONDITION="${ENTRY#*|}"
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"

  echo "  Scoring: ${CONDITION} (${PGS_ID})..."

  # Convert the harmonised PGS Catalog file to plink2 --score input:
  # chr:pos variant ID (GRCh38 hm_chr/hm_pos), effect allele, weight.
  # Rows the catalog could not map to GRCh38 have an empty hm_pos and are skipped.
  FORMATTED="${OUTDIR}/${PGS_ID}_formatted.tsv"
  gzip -cd "$SCORE_FILE" | \
    awk -F'\t' '/^#/ {next}
    !hdr {
      for(i=1;i<=NF;i++) {
        if($i=="hm_chr") chr_col=i;
        if($i=="hm_pos") pos_col=i;
        if($i=="effect_allele") ea_col=i;
        if($i=="effect_weight") ew_col=i;
      }
      hdr=1
      next
    }
    chr_col && pos_col && ea_col && ew_col {
      chr=$chr_col; pos=$pos_col; ea=$ea_col; ew=$ew_col;
      if(chr!="" && pos!="" && ea!="" && ew!="") {
        # Add chr prefix if missing to match GRCh38 VCF contig names
        if(chr !~ /^chr/) chr="chr"chr;
        key=chr":"pos"\t"ea;
        if(!(key in seen)) {
          seen[key]=1;
          printf "%s:%s\t%s\t%s\n", chr, pos, ea, ew;
        }
      }
    }' > "$FORMATTED"

  TOTAL_VARS=$(wc -l < "$FORMATTED" | tr -d ' ')

  if [ "$TOTAL_VARS" -eq 0 ]; then
    echo "ERROR: No GRCh38 hm_chr/hm_pos/effect_allele/effect_weight rows in ${SCORE_FILE}" >&2
    exit 1
  fi

  # Remove any score left by an earlier run, so a failed plink2 run cannot be
  # reported with an old number.
  SSCORE="${OUTDIR}/${PGS_ID}.sscore"
  rm -f "$SSCORE" "${OUTDIR}/${PGS_ID}.log"

  # cols=+scoresums adds SCORE1_SUM: the plain weighted sum. The default
  # SCORE1_AVG divides by the alleles present in this VCF, which differs per sample.
  if ! run_in    --cpus 4 --memory 4g \
    "${PLINK2_IMAGE}" \
    plink2 \
      --pfile "/genome/${SAMPLE}/prs/${SAMPLE}" \
      --score "/genome/${SAMPLE}/prs/${PGS_ID}_formatted.tsv" 1 2 3 \
        ignore-dup-ids \
        no-mean-imputation \
        cols=+scoresums \
      --out "/genome/${SAMPLE}/prs/${PGS_ID}" \
      --threads 4 \
      --memory 3000 \
      --allow-extra-chr; then
    # None of the score's variants is in this VCF: report that, not a failure.
    if grep -q 'No valid variants' "${OUTDIR}/${PGS_ID}.log" 2>/dev/null; then
      echo -e "${CONDITION}\t${PGS_ID}\tNA\t0\t${TOTAL_VARS}" >> "$RESULTS_FILE"
      echo "    No variant of ${PGS_ID} is present in this VCF; no score."
      continue
    fi
    echo "ERROR: plink2 --score failed for ${PGS_ID}; see ${OUTDIR}/${PGS_ID}.log" >&2
    exit 1
  fi

  if [ ! -s "$SSCORE" ]; then
    echo "ERROR: plink2 exited 0 but wrote no ${SSCORE}" >&2
    exit 1
  fi

  # Read columns by name. ALLELE_CT counts the alleles scored (two per matched
  # variant on the autosomes), so ALLELE_CT/2 is the number of matched variants.
  read -r SCORE USED_VARS < <(awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) col[$i]=i; next }
    NR==2 {
      if(!("SCORE1_SUM" in col) || !("ALLELE_CT" in col)) { print "MISSING MISSING"; exit }
      print $col["SCORE1_SUM"], $col["ALLELE_CT"]/2
    }' "$SSCORE")
  if [ "$SCORE" = "MISSING" ] || [ -z "$SCORE" ]; then
    echo "ERROR: ${SSCORE} has no SCORE1_SUM or ALLELE_CT column" >&2
    exit 1
  fi
  echo -e "${CONDITION}\t${PGS_ID}\t${SCORE}\t${USED_VARS}\t${TOTAL_VARS}" >> "$RESULTS_FILE"
  echo "    Score (sum): ${SCORE} (${USED_VARS}/${TOTAL_VARS} variants matched)"
  echo "    ${HOMREF_NOTE}"
done

echo ""
echo "============================================"
echo "  Polygenic Risk Scores complete: ${SAMPLE}"
echo ""
echo "  Summary: ${RESULTS_FILE}"
column -t -s $'\t' "$RESULTS_FILE" 2>/dev/null || cat "$RESULTS_FILE"
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
