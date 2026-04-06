#!/usr/bin/env bash
# generate-report.sh — Create a plain-text summary report of all pipeline results
# Usage: ./scripts/generate-report.sh <sample_name>
#
# Scans all output directories and produces a readable summary.
# Run after completing all (or some) pipeline steps.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
SAMPLE_DIR="${GENOME_DIR}/${SAMPLE}"
REPORT="${SAMPLE_DIR}/${SAMPLE}_report.txt"

# Check sample directory exists
if [ ! -d "$SAMPLE_DIR" ]; then
  echo "ERROR: Sample directory not found: ${SAMPLE_DIR}" >&2
  exit 1
fi

{
echo "================================================================================"
echo "  GENOMICS ANALYSIS REPORT"
echo "  Sample: ${SAMPLE}"
echo "  Generated: $(date -u '+%Y-%m-%d %H:%M UTC')"
echo "================================================================================"
echo ""

# ---------- Variant Calling (Step 3) ----------
VCF="${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz"
if [ -f "$VCF" ]; then
  echo "## Variant Calling (DeepVariant)"
  echo "---"
  TOTAL=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools stats "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
    grep '^SN.*number of records' | awk '{print $NF}' || echo "N/A")
  PASS=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" 2>/dev/null | \
    grep -cv '^#' || echo "N/A")
  echo "  Total variants: ${TOTAL}"
  echo "  PASS variants:  ${PASS}"
  echo ""
fi

# ---------- ClinVar Screen (Step 6) ----------
CLINVAR_DIR="${SAMPLE_DIR}/clinvar"
if [ -d "$CLINVAR_DIR" ] && [ -f "${CLINVAR_DIR}/isec/0002.vcf" ]; then
  echo "## ClinVar Pathogenic Screen"
  echo "---"
  HITS=$(grep -cv '^#' "${CLINVAR_DIR}/isec/0002.vcf" 2>/dev/null || echo "0")
  echo "  Pathogenic/Likely Pathogenic hits: ${HITS}"
  if [ "$HITS" -gt 0 ]; then
    echo ""
    echo "  Genes affected:"
    grep -v '^#' "${CLINVAR_DIR}/isec/0002.vcf" 2>/dev/null | head -20 | while IFS=$'\t' read -r chr pos id ref alt rest; do
      echo "    ${chr}:${pos} ${ref}>${alt} (${id})"
    done
  fi
  echo ""
fi

# ---------- PharmCAT (Step 7) ----------
# PharmCAT writes reports alongside the VCF in vcf/
PHARMCAT_REPORT=""
for DIR in "${SAMPLE_DIR}/vcf" "${SAMPLE_DIR}/pharmcat"; do
  [ -d "$DIR" ] || continue
  PHARMCAT_REPORT=$(find "$DIR" -maxdepth 1 -name "*.report.json" 2>/dev/null | head -1)
  [ -n "$PHARMCAT_REPORT" ] && break
done
if [ -n "$PHARMCAT_REPORT" ]; then
  echo "## Pharmacogenomics (PharmCAT)"
  echo "---"
  echo "  Report: $(basename "$PHARMCAT_REPORT")"
  PHARMCAT_HTML=$(find "$(dirname "$PHARMCAT_REPORT")" -maxdepth 1 -name "*.report.html" 2>/dev/null | head -1)
  if [ -n "$PHARMCAT_HTML" ]; then
    echo "  HTML report: $(basename "$PHARMCAT_HTML") (open in browser for full details)"
  fi
  echo ""
fi

# ---------- Manta SVs (Step 4) ----------
MANTA_DIR="${SAMPLE_DIR}/manta"
if [ -d "$MANTA_DIR" ]; then
  SV_VCF=$(find "$MANTA_DIR" -name "diploidSV.vcf.gz" 2>/dev/null | head -1)
  if [ -n "$SV_VCF" ]; then
    echo "## Structural Variants (Manta)"
    echo "---"
    SV_COUNT=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
      bcftools view "${SV_VCF/#$GENOME_DIR//genome}" 2>/dev/null | grep -cv '^#' || echo "N/A")
    echo "  Total SVs: ${SV_COUNT}"
    echo ""
  fi
fi

# ---------- ExpansionHunter (Step 9) ----------
EH_DIR="${SAMPLE_DIR}/expansion_hunter"
if [ -d "$EH_DIR" ]; then
  EH_JSON=$(find "$EH_DIR" -name "*_eh.json" 2>/dev/null | head -1)
  if [ -n "$EH_JSON" ]; then
    echo "## Repeat Expansions (ExpansionHunter)"
    echo "---"
    # EH JSON uses locus names as top-level keys; count objects with a Genotype field
    LOCI=$(jq '[to_entries[] | select(.value | type=="object" and has("Genotype"))] | length' "$EH_JSON" 2>/dev/null || grep -c '"Genotype"' "$EH_JSON" 2>/dev/null || echo "0")
    echo "  Loci tested: ${LOCI}"
    echo "  (See interpreting-results.md for disease thresholds)"
    echo ""
  fi
fi

# ---------- TelomereHunter (Step 10) ----------
TEL_DIR="${SAMPLE_DIR}/telomere"
if [ -d "$TEL_DIR" ]; then
  SUMMARY=$(find "$TEL_DIR" -name "*_summary.tsv" 2>/dev/null | head -1)
  if [ -n "$SUMMARY" ]; then
    echo "## Telomere Length (TelomereHunter)"
    echo "---"
    TEL_CONTENT=$(awk -F'\t' 'NR==2 {print $11}' "$SUMMARY" 2>/dev/null || echo "N/A")
    echo "  Telomere content: ${TEL_CONTENT}"
    echo ""
  fi
fi

# ---------- ROH (Step 11) ----------
ROH_FILE="${SAMPLE_DIR}/vcf/${SAMPLE}_roh.txt"
if [ -f "$ROH_FILE" ]; then
  echo "## Runs of Homozygosity"
  echo "---"
  LARGE_ROH=$(grep '^RG' "$ROH_FILE" 2>/dev/null | awk '$3 !~ /chrX|chrY/ && $6 > 5000000' | wc -l || echo "0")
  echo "  Autosomal ROH > 5MB: ${LARGE_ROH}"
  if [ "$LARGE_ROH" -gt 0 ]; then
    grep '^RG' "$ROH_FILE" | awk '$3 !~ /chrX|chrY/ && $6 > 5000000 {printf "    %s:%s-%s  %.1fMB\n", $3,$4,$5,$6/1e6}'
  fi
  echo ""
fi

# ---------- Mito Haplogroup (Step 12) ----------
HAPLO_FILE="${SAMPLE_DIR}/mito/${SAMPLE}_haplogroup.txt"
if [ -f "$HAPLO_FILE" ]; then
  echo "## Mitochondrial Haplogroup"
  echo "---"
  HAPLO=$(awk -F'\t' 'NR==2 {gsub(/"/, "", $2); print $2}' "$HAPLO_FILE" 2>/dev/null || echo "N/A")
  echo "  Haplogroup: ${HAPLO}"
  echo ""
fi

# ---------- CPSR (Step 17) ----------
CPSR_DIR="${SAMPLE_DIR}/cpsr"
if [ -d "$CPSR_DIR" ]; then
  CPSR_HTML=$(find "$CPSR_DIR" -name "*.cpsr.*.html" 2>/dev/null | head -1)
  if [ -n "$CPSR_HTML" ]; then
    echo "## Cancer Predisposition Screening (CPSR)"
    echo "---"
    echo "  HTML report: $(basename "$CPSR_HTML") (open in browser)"
    CPSR_TSV=$(find "$CPSR_DIR" -name "*.snvs_indels.tiers.tsv" 2>/dev/null | head -1)
    if [ -n "$CPSR_TSV" ]; then
      echo "  Classification breakdown:"
      tail -n +2 "$CPSR_TSV" | cut -f88 | sort | uniq -c | sort -rn | while read -r count class; do
        echo "    ${class}: ${count}"
      done
    fi
    echo ""
  fi
fi

# ---------- CNVnator (Step 18) ----------
CNV_FILE="${SAMPLE_DIR}/cnvnator/${SAMPLE}_cnvs.txt"
if [ -f "$CNV_FILE" ]; then
  echo "## Copy Number Variants (CNVnator)"
  echo "---"
  TOTAL_CNV=$(wc -l < "$CNV_FILE")
  SIG_CNV=$(awk '$5 < 0.01' "$CNV_FILE" | wc -l)
  DEL_CNV=$(grep -c '^deletion' "$CNV_FILE" || echo "0")
  DUP_CNV=$(grep -c '^duplication' "$CNV_FILE" || echo "0")
  echo "  Total CNVs: ${TOTAL_CNV} (${DEL_CNV} deletions, ${DUP_CNV} duplications)"
  echo "  Significant (e-val < 0.01): ${SIG_CNV}"
  echo ""
fi

# ---------- Delly (Step 19) ----------
DELLY_VCF="${SAMPLE_DIR}/delly/${SAMPLE}_sv.vcf.gz"
if [ -f "$DELLY_VCF" ]; then
  echo "## Structural Variants (Delly)"
  echo "---"
  TOTAL_DELLY=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view "/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz" 2>/dev/null | grep -cv '^#' || echo "N/A")
  PASS_DELLY=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS "/genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz" 2>/dev/null | grep -cv '^#' || echo "N/A")
  echo "  Total SVs: ${TOTAL_DELLY}"
  echo "  PASS SVs: ${PASS_DELLY}"
  echo ""
fi

# ---------- Mitochondrial (Step 20) ----------
MITO_VCF="${SAMPLE_DIR}/mito/${SAMPLE}_chrM_filtered.vcf.gz"
if [ -f "$MITO_VCF" ]; then
  echo "## Mitochondrial Variants (Mutect2)"
  echo "---"
  TOTAL_MITO=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" 2>/dev/null | grep -cv '^#' || echo "N/A")
  HETERO=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -f PASS "/genome/${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz" 2>/dev/null | \
    docker run --rm -i staphb/bcftools:1.21 bcftools query -f '[%AF]\n' 2>/dev/null | \
    awk '$1 < 0.95' | wc -l || echo "N/A")
  echo "  PASS variants: ${TOTAL_MITO}"
  echo "  Heteroplasmic (AF < 0.95): ${HETERO}"
  echo ""
fi

# ---------- SV Consensus Merge (Step 22) ----------
SV_MERGE_VCF="${SAMPLE_DIR}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz"
if [ -f "$SV_MERGE_VCF" ]; then
  echo "## SV Consensus Merge"
  echo "---"
  CONSENSUS_COUNT=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
  echo "  Consensus SVs (2+ callers): ${CONSENSUS_COUNT}"
  echo ""
fi

# ---------- Clinical Filter (Step 23) ----------
CLINICAL_VCF="${SAMPLE_DIR}/clinical/${SAMPLE}_clinical.vcf.gz"
if [ -f "$CLINICAL_VCF" ]; then
  echo "## Clinical Variant Filter"
  echo "---"
  CLINICAL_COUNT=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
    bcftools view -H "/genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz" 2>/dev/null | wc -l || echo "N/A")
  echo "  Clinical variants: ${CLINICAL_COUNT}"
  echo ""
fi

# ---------- PRS (Step 25) ----------
PRS_SUMMARY="${SAMPLE_DIR}/prs/${SAMPLE}_prs_summary.tsv"
if [ -f "$PRS_SUMMARY" ] && [ "$(wc -l < "$PRS_SUMMARY")" -gt 1 ]; then
  echo "## Polygenic Risk Scores"
  echo "---"
  tail -n +2 "$PRS_SUMMARY" | while IFS=$'\t' read -r CONDITION _PGS_ID SCORE USED TOTAL; do
    printf "  %-35s %s (%s/%s variants)\n" "$CONDITION" "$SCORE" "$USED" "$TOTAL"
  done
  echo ""
  echo "  NOTE: Raw PRS scores are NOT directly interpretable without a"
  echo "  population reference panel. See docs/25-prs.md."
  echo ""
fi

# ---------- CPIC (Step 27) ----------
CPIC_REPORT="${SAMPLE_DIR}/cpic/${SAMPLE}_cpic_recommendations.txt"
if [ -f "$CPIC_REPORT" ]; then
  echo "## CPIC Drug-Gene Recommendations"
  echo "---"
  echo "  Report: ${CPIC_REPORT}"
  AFFECTED=$(grep -c 'Affected drugs:' "$CPIC_REPORT" 2>/dev/null || echo "0")
  echo "  Genes with non-standard phenotypes requiring drug adjustments: ${AFFECTED}"
  echo ""
fi

# ---------- Steps Not Run ----------
echo "## Steps Not Run"
echo "---"
NOT_RUN=""
[ ! -f "${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz" ] && NOT_RUN="${NOT_RUN}  - Variant Calling (step 3)\n"
[ ! -d "${SAMPLE_DIR}/clinvar" ] && NOT_RUN="${NOT_RUN}  - ClinVar Screen (step 6)\n"
PHARMCAT_FOUND=0; find "${SAMPLE_DIR}/vcf" -maxdepth 1 -name "*.report.html" 2>/dev/null | grep -q . && PHARMCAT_FOUND=1; find "${SAMPLE_DIR}/pharmcat" -maxdepth 1 -name "*.report.html" 2>/dev/null | grep -q . && PHARMCAT_FOUND=1; [ "$PHARMCAT_FOUND" -eq 0 ] && NOT_RUN="${NOT_RUN}  - PharmCAT (step 7)\n"
[ ! -d "${SAMPLE_DIR}/manta" ] && NOT_RUN="${NOT_RUN}  - Manta SVs (step 4)\n"
[ ! -d "${SAMPLE_DIR}/expansion_hunter" ] && NOT_RUN="${NOT_RUN}  - ExpansionHunter (step 9)\n"
[ ! -d "${SAMPLE_DIR}/cpsr" ] && NOT_RUN="${NOT_RUN}  - CPSR (step 17)\n"
[ ! -d "${SAMPLE_DIR}/vep" ] && NOT_RUN="${NOT_RUN}  - VEP Annotation (step 13)\n"
[ ! -f "${SAMPLE_DIR}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz" ] && NOT_RUN="${NOT_RUN}  - SV Consensus Merge (step 22)\n"
[ ! -f "${SAMPLE_DIR}/clinical/${SAMPLE}_clinical.vcf.gz" ] && NOT_RUN="${NOT_RUN}  - Clinical Filter (step 23)\n"
[ ! -f "${SAMPLE_DIR}/prs/${SAMPLE}_prs_summary.tsv" ] && NOT_RUN="${NOT_RUN}  - PRS (step 25)\n"
[ ! -d "${SAMPLE_DIR}/ancestry" ] && NOT_RUN="${NOT_RUN}  - Ancestry PCA (step 26)\n"
[ ! -f "${SAMPLE_DIR}/cpic/${SAMPLE}_cpic_recommendations.txt" ] && NOT_RUN="${NOT_RUN}  - CPIC Recommendations (step 27)\n"
if [ -z "$NOT_RUN" ]; then
  echo "  All major steps completed."
else
  printf "%b" "$NOT_RUN"
fi
echo ""

echo "================================================================================"
echo "  DISCLAIMER: This is NOT a clinical report. Discuss findings with a"
echo "  qualified healthcare professional before making medical decisions."
echo "================================================================================"

} | tee "$REPORT"

echo ""
echo "Report saved to: ${REPORT}"
