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
PHARMCAT_DIR="${SAMPLE_DIR}/pharmcat"
if [ -d "$PHARMCAT_DIR" ]; then
  REPORT_JSON=$(find "$PHARMCAT_DIR" -name "*.report.json" 2>/dev/null | head -1)
  if [ -n "$REPORT_JSON" ]; then
    echo "## Pharmacogenomics (PharmCAT)"
    echo "---"
    echo "  Report: $(basename "$REPORT_JSON")"
    HTML=$(find "$PHARMCAT_DIR" -name "*.report.html" 2>/dev/null | head -1)
    if [ -n "$HTML" ]; then
      echo "  HTML report: $(basename "$HTML") (open in browser for full details)"
    fi
    echo ""
  fi
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
    LOCI=$(grep -c '"LocusId"' "$EH_JSON" 2>/dev/null || echo "0")
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
    TEL_CONTENT=$(awk -F'\t' 'NR==2 {print $2}' "$SUMMARY" 2>/dev/null || echo "N/A")
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
  HAPLO=$(cat "$HAPLO_FILE" 2>/dev/null || echo "N/A")
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

# ---------- Steps Not Run ----------
echo "## Steps Not Run"
echo "---"
NOT_RUN=""
[ ! -f "${SAMPLE_DIR}/vcf/${SAMPLE}.vcf.gz" ] && NOT_RUN="${NOT_RUN}  - Variant Calling (step 3)\n"
[ ! -d "${SAMPLE_DIR}/clinvar" ] && NOT_RUN="${NOT_RUN}  - ClinVar Screen (step 6)\n"
[ ! -d "${SAMPLE_DIR}/pharmcat" ] && NOT_RUN="${NOT_RUN}  - PharmCAT (step 7)\n"
[ ! -d "${SAMPLE_DIR}/manta" ] && NOT_RUN="${NOT_RUN}  - Manta SVs (step 4)\n"
[ ! -d "${SAMPLE_DIR}/expansion_hunter" ] && NOT_RUN="${NOT_RUN}  - ExpansionHunter (step 9)\n"
[ ! -d "${SAMPLE_DIR}/cpsr" ] && NOT_RUN="${NOT_RUN}  - CPSR (step 17)\n"
[ ! -d "${SAMPLE_DIR}/vep" ] && NOT_RUN="${NOT_RUN}  - VEP Annotation (step 13)\n"
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
