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
# Requires: VCF from step 3. With step 3's gVCF next to it, the score
# positions are genotyped from the gVCF, so a site where you match the
# reference is a real 0/0 (a dosage of 2 for a reference effect allele)
# instead of a missing site. Without it the sum leaves those sites out.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
GVCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz"
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
  # The PGS Catalog publishes an md5 next to every scoring file.
  if ! fetch "$URL" "$SCORE_FILE" md5 "${URL}.md5"; then
    echo "ERROR: Could not download the GRCh38-harmonised scoring file for ${PGS_ID}:" >&2
    echo "  ${URL}" >&2
    exit 1
  fi
  BUILD=$(hm_build "$SCORE_FILE")
  if [ "$BUILD" != "GRCh38" ]; then
    rm -f "$SCORE_FILE"
    echo "ERROR: ${PGS_ID} scoring file has #HmPOS_build='${BUILD}', expected GRCh38. Refusing to score it." >&2
    exit 1
  fi
done

# Each score as plink2 --score input: chr:pos variant ID (GRCh38 hm_chr/hm_pos),
# effect allele, weight. Rows the catalog could not map to GRCh38 have an
# empty hm_pos and are skipped.
for ENTRY in "${PGS_SCORES[@]}"; do
  PGS_ID="${ENTRY%%|*}"
  gzip -cd "${SCORING_DIR}/${PGS_ID}.txt.gz" | \
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
    }' > "${OUTDIR}/${PGS_ID}_formatted.tsv"
  if [ ! -s "${OUTDIR}/${PGS_ID}_formatted.tsv" ]; then
    echo "ERROR: No GRCh38 hm_chr/hm_pos/effect_allele/effect_weight rows in ${SCORING_DIR}/${PGS_ID}.txt.gz" >&2
    exit 1
  fi
done

echo ""
if [ -f "$GVCF" ] && [ -f "${GVCF}.tbi" ]; then
  INPUT_KIND=gvcf
  echo "[2/3] Genotyping the score positions from the gVCF (${GVCF})..."
  # Every score position, once, as CHROM<TAB>POS for bcftools -R/-T, and the
  # effect alleles of each position for the step below.
  SITES="${OUTDIR}/score_sites.tsv"
  ALLELES="${OUTDIR}/score_alleles.tsv"
  FORMATTED_FILES=()
  for ENTRY in "${PGS_SCORES[@]}"; do FORMATTED_FILES+=("${OUTDIR}/${ENTRY%%|*}_formatted.tsv"); done
  cat "${FORMATTED_FILES[@]}" | awk -F'\t' '{split($1, a, ":"); print a[1] "\t" a[2] "\t" $2}' \
    | sort -u -k1,1 -k2,2n -k3,3 > "$ALLELES"
  cut -f1,2 "$ALLELES" | uniq > "$SITES"
  # gvcf2vcf expands each reference block overlapping a score position into
  # one 0/0 record per base, with the base from the reference; -T keeps the
  # score positions, --trim-alt-alleles drops the <*> allele, and a no-call
  # (./.) is dropped, so a position without coverage stays missing.
  # A 0/0 record then has ALT '.', which no effect allele can match: the awk
  # sets ALT to the position's first effect allele that is not the reference,
  # so a score whose effect allele is not the reference matches it with a
  # dosage of 0 and counts as matched.
  SITES_VCF="${OUTDIR}/${SAMPLE}_score_sites.vcf"
  rm -f "${SITES_VCF}.tmp"
  # shellcheck disable=SC2016  # $1 to $3 belong to the inner bash
  run_in --cpus 2 --memory 8g \
    "${BCFTOOLS_IMAGE}" \
    bash -euo pipefail -c '
      bcftools convert --gvcf2vcf -f "$1" -R "$2" -Ou "$3" \
        | bcftools view -T "$2" --trim-alt-alleles -i "GT!=\"mis\"" -Ov' \
    _ "${REF_FASTA_C}" "$(cpath "$SITES")" "/genome/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" \
  | awk -F'\t' -v OFS='\t' -v alleles="$ALLELES" '
      BEGIN { while ((getline l < alleles) > 0) { split(l, f, "\t"); k = f[1] ":" f[2]; ea[k] = ea[k] " " f[3] } }
      /^#/ { print; next }
      $5 == "." {
        n = split(ea[$1 ":" $2], c, " ")
        for (i = 1; i <= n; i++) if (c[i] != $4) { $5 = c[i]; break }
      }
      { print }' > "${SITES_VCF}.tmp"
  mv -f "${SITES_VCF}.tmp" "$SITES_VCF"
  PLINK_VCF="/genome/${SAMPLE}/prs/${SAMPLE}_score_sites.vcf"
else
  INPUT_KIND=vcf
  echo "[2/3] Converting VCF to plink2 format (no gVCF: sites where you match the reference are missing)..."
  PLINK_VCF="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
fi

# Convert to plink2 binary format for scoring
run_in --cpus 4 --memory 8g \
  "${PLINK2_IMAGE}" \
  plink2 \
    --vcf "$PLINK_VCF" \
    --make-pgen \
    --out "/genome/${SAMPLE}/prs/${SAMPLE}" \
    --threads 4 \
    --memory 6000 \
    --set-all-var-ids '@:#' \
    --new-id-max-allele-len 100 \
    --chr 1-22 \
    --allow-extra-chr \
    --output-chr chrM
if [ "$INPUT_KIND" = gvcf ]; then
  rm -f "${OUTDIR}/${SAMPLE}_score_sites.vcf" "${OUTDIR}/score_sites.tsv" "${OUTDIR}/score_alleles.tsv"
fi

echo ""
echo "[3/3] Calculating polygenic risk scores..."

HOMREF_NOTE="hom-ref sites are absent from this VCF, so the score is biased; not comparable to published distributions"

# Process each scoring file. Matched_Pct is Variants_Matched / Variants_Total;
# Input says what was scored: gvcf (hom-ref sites included) or vcf (variant
# sites only, the biased sum).
RESULTS_FILE="${OUTDIR}/${SAMPLE}_prs_summary.tsv"
echo -e "Condition\tPGS_ID\tScore_SUM\tVariants_Matched\tVariants_Total\tMatched_Pct\tInput" > "$RESULTS_FILE"

for ENTRY in "${PGS_SCORES[@]}"; do
  PGS_ID="${ENTRY%%|*}"
  CONDITION="${ENTRY#*|}"
  SCORE_FILE="${SCORING_DIR}/${PGS_ID}.txt.gz"

  echo "  Scoring: ${CONDITION} (${PGS_ID})..."

  FORMATTED="${OUTDIR}/${PGS_ID}_formatted.tsv"
  TOTAL_VARS=$(wc -l < "$FORMATTED" | tr -d ' ')

  # Remove any score left by an earlier run, so a failed plink2 run cannot be
  # reported with an old number.
  SSCORE="${OUTDIR}/${PGS_ID}.sscore"
  rm -f "$SSCORE" "${OUTDIR}/${PGS_ID}.log"

  # cols=+scoresums adds SCORE1_SUM: the plain weighted sum. The default
  # SCORE1_AVG divides by the alleles present in this VCF, which differs per sample.
  if ! run_in --cpus 4 --memory 4g \
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
      echo -e "${CONDITION}\t${PGS_ID}\tNA\t0\t${TOTAL_VARS}\t0.0\t${INPUT_KIND}" >> "$RESULTS_FILE"
      echo "    No variant of ${PGS_ID} is present in this input; no score."
      if [ "$INPUT_KIND" = vcf ]; then
        echo "    ${HOMREF_NOTE}"
      fi
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
  PCT=$(awk -v u="$USED_VARS" -v t="$TOTAL_VARS" 'BEGIN { printf "%.1f", 100 * u / t }')
  echo -e "${CONDITION}\t${PGS_ID}\t${SCORE}\t${USED_VARS}\t${TOTAL_VARS}\t${PCT}\t${INPUT_KIND}" >> "$RESULTS_FILE"
  echo "    Score (sum): ${SCORE}"
  echo "    Matched: ${USED_VARS} of ${TOTAL_VARS} score variants (${PCT}%)"
  if awk -v p="$PCT" 'BEGIN { exit !(p < 50) }'; then
    echo "    WARNING: under half of the score's variants were genotyped; the sum is not comparable to published distributions."
  fi
  if [ "$INPUT_KIND" = vcf ]; then
    echo "    ${HOMREF_NOTE}"
  fi
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
