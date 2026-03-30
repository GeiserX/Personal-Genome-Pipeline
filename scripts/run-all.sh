#!/usr/bin/env bash
# run-all.sh — Run the complete genomics analysis pipeline for one sample
# Usage: ./run-all.sh <sample_name> <sex: male|female>
#
# Assumes:
# - FASTQ files at $GENOME_DIR/<sample>/fastq/ OR
# - BAM already exists at $GENOME_DIR/<sample>/aligned/<sample>_sorted.bam
# - Reference genome at $GENOME_DIR/reference/
# - ClinVar database at $GENOME_DIR/clinvar/
#
# Steps run in parallel where possible. Total time: ~6-12 hours on 16 cores.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <sex: male|female>}
SEX=${2:?Usage: $0 <sample_name> <sex: male|female>}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

PIPELINE_START=$(date +%s)

echo "============================================"
echo "  Genomics Pipeline — Full Analysis"
echo "  Sample: ${SAMPLE}, Sex: ${SEX}"
echo "  Data: ${GENOME_DIR}/${SAMPLE}/"
echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""

# Pre-flight check
echo "[Pre-flight] Validating setup..."
if ! "${SCRIPT_DIR}/validate-setup.sh" "${SAMPLE}" 2>/dev/null; then
  echo "WARNING: Setup validation reported issues. Proceeding anyway..."
  echo "         Run ./scripts/validate-setup.sh ${SAMPLE} for details."
  echo ""
fi

# Phase 1: Alignment (if FASTQ exists but BAM doesn't)
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
if [ ! -f "$BAM" ]; then
  echo "[Phase 1] Alignment — FASTQ to sorted BAM..."
  bash "${SCRIPT_DIR}/02-alignment.sh" "$SAMPLE"
else
  echo "[Phase 1] BAM already exists, skipping alignment."
fi
echo ""

# Phase 2: Variant Calling (if VCF doesn't exist)
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
if [ ! -f "$VCF" ]; then
  echo "[Phase 2] Variant calling — DeepVariant..."
  bash "${SCRIPT_DIR}/03-deepvariant.sh" "$SAMPLE"
else
  echo "[Phase 2] VCF already exists, skipping variant calling."
fi
echo ""

# Phase 3: Parallel analyses (all independent after BAM + VCF exist)
echo "[Phase 3] Running parallel analyses..."
echo ""

# --- Group A: Quick jobs (minutes each) ---
echo "  Starting quick analyses..."

echo "  [A1] ClinVar screen..."
bash "${SCRIPT_DIR}/06-clinvar-screen.sh" "$SAMPLE" &
PID_CLINVAR=$!

echo "  [A2] PharmCAT pharmacogenomics..."
bash "${SCRIPT_DIR}/07-pharmacogenomics.sh" "$SAMPLE" &
PID_PHARMCAT=$!

echo "  [A3] ROH analysis..."
bash "${SCRIPT_DIR}/11-roh-analysis.sh" "$SAMPLE" &
PID_ROH=$!

echo "  [A4] Mito haplogroup..."
bash "${SCRIPT_DIR}/12-mito-haplogroup.sh" "$SAMPLE" &
PID_HAPLO=$!

echo "  [A5] indexcov coverage QC..."
bash "${SCRIPT_DIR}/16-indexcov.sh" "$SAMPLE" "$SEX" &
PID_INDEXCOV=$!

echo "  [A6] Imputation prep..."
bash "${SCRIPT_DIR}/14-imputation-prep.sh" "$SAMPLE" &
PID_IMPUTATION=$!

echo "  [A7] HLA typing (T1K)..."
bash "${SCRIPT_DIR}/08-hla-typing.sh" "$SAMPLE" &
PID_HLA=$!

# --- Group B: Medium jobs (10-60 minutes each) ---
echo "  Starting medium analyses..."

echo "  [B1] Manta structural variants..."
bash "${SCRIPT_DIR}/04-manta.sh" "$SAMPLE" &
PID_MANTA=$!

echo "  [B2] ExpansionHunter STR screening..."
bash "${SCRIPT_DIR}/09-expansion-hunter.sh" "$SAMPLE" "$SEX" &
PID_EH=$!

echo "  [B3] TelomereHunter telomere length..."
bash "${SCRIPT_DIR}/10-telomere-hunter.sh" "$SAMPLE" &
PID_TH=$!

echo "  [B4] MToolBox mitochondrial analysis..."
bash "${SCRIPT_DIR}/20-mtoolbox.sh" "$SAMPLE" &
PID_MTOOLBOX=$!

echo "  [B5] CPSR cancer predisposition..."
bash "${SCRIPT_DIR}/17-cpsr.sh" "$SAMPLE" &
PID_CPSR=$!

# Wait for quick jobs
wait $PID_CLINVAR $PID_PHARMCAT $PID_ROH $PID_HAPLO $PID_INDEXCOV $PID_IMPUTATION $PID_HLA 2>/dev/null || true
echo ""
echo "  Quick analyses complete."

# --- Group C: Heavy jobs (2-4 hours each) ---
# These are CPU+RAM intensive — run sequentially or limit parallelism
echo "  Starting heavy analyses..."

echo "  [C1] VEP functional annotation..."
bash "${SCRIPT_DIR}/13-vep-annotation.sh" "$SAMPLE" &
PID_VEP=$!

echo "  [C2] CNVnator depth-based CNV calling..."
bash "${SCRIPT_DIR}/18-cnvnator.sh" "$SAMPLE" &
PID_CNVNATOR=$!

echo "  [C3] Delly structural variant calling..."
bash "${SCRIPT_DIR}/19-delly.sh" "$SAMPLE" &
PID_DELLY=$!

# Wait for Manta before running duphold and AnnotSV
wait $PID_MANTA 2>/dev/null || true
echo "  Manta complete. Running SV post-processing..."

echo "  [B6] duphold SV quality annotation..."
bash "${SCRIPT_DIR}/15-duphold.sh" "$SAMPLE" &
PID_DUPHOLD=$!

echo "  [B7] AnnotSV structural variant annotation..."
bash "${SCRIPT_DIR}/05-annotsv.sh" "$SAMPLE" &
PID_ANNOTSV=$!

# Wait for everything
wait $PID_EH $PID_TH $PID_MTOOLBOX $PID_CPSR 2>/dev/null || true
wait $PID_VEP $PID_CNVNATOR $PID_DELLY 2>/dev/null || true
wait $PID_DUPHOLD $PID_ANNOTSV 2>/dev/null || true

# Phase 4: Post-processing (uses outputs from Phase 3)
echo ""
echo "[Phase 4] Running post-processing steps..."

echo "  [D1] CYP2D6 star alleles (Cyrius)..."
bash "${SCRIPT_DIR}/21-cyrius.sh" "$SAMPLE" 2>/dev/null &
PID_CYRIUS=$!

echo "  [D2] SV consensus merge..."
bash "${SCRIPT_DIR}/22-survivor-merge.sh" "$SAMPLE" 2>/dev/null &
PID_SURVIVOR=$!

echo "  [D3] Clinical variant filter..."
bash "${SCRIPT_DIR}/23-clinical-filter.sh" "$SAMPLE" 2>/dev/null &
PID_CLINICAL=$!

echo "  [D4] Polygenic Risk Scores..."
bash "${SCRIPT_DIR}/25-prs.sh" "$SAMPLE" 2>/dev/null &
PID_PRS=$!

echo "  [D5] Ancestry PCA..."
bash "${SCRIPT_DIR}/26-ancestry.sh" "$SAMPLE" 2>/dev/null &
PID_ANCESTRY=$!

echo "  [D6] CPIC drug-gene recommendations..."
bash "${SCRIPT_DIR}/27-cpic-lookup.sh" "$SAMPLE" 2>/dev/null &
PID_CPIC=$!

wait $PID_CYRIUS $PID_SURVIVOR $PID_CLINICAL $PID_PRS $PID_ANCESTRY $PID_CPIC 2>/dev/null || true
echo "  Post-processing complete."

echo "  [D7] HTML summary report..."
bash "${SCRIPT_DIR}/24-html-report.sh" "$SAMPLE" 2>/dev/null || true

# Generate summary report
echo ""
echo "[Report] Generating summary report..."
bash "${SCRIPT_DIR}/generate-report.sh" "$SAMPLE" 2>/dev/null || echo "  (report generation had warnings — check output manually)"

PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

echo ""
echo "============================================"
echo "  Pipeline complete for: ${SAMPLE}"
echo "  All results in: ${GENOME_DIR}/${SAMPLE}/"
echo "  Total runtime: ${HOURS}h ${MINUTES}m"
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""
echo "Key outputs:"
echo "  HTML Report:    ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"
echo "  Text Report:    ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt"
echo "  VCF:            ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
echo "  ClinVar hits:   ${GENOME_DIR}/${SAMPLE}/clinvar/"
echo "  PharmCAT:       ${GENOME_DIR}/${SAMPLE}/pharmcat/"
echo "  CYP2D6:         ${GENOME_DIR}/${SAMPLE}/cyrius/"
echo "  CPIC drugs:     ${GENOME_DIR}/${SAMPLE}/cpic/"
echo "  Clinical VCF:   ${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_clinical.vcf.gz"
echo "  SV consensus:   ${GENOME_DIR}/${SAMPLE}/sv_merged/"
echo "  PRS scores:     ${GENOME_DIR}/${SAMPLE}/prs/"
echo "  Ancestry PCA:   ${GENOME_DIR}/${SAMPLE}/ancestry/"
echo "  CPSR report:    ${GENOME_DIR}/${SAMPLE}/cpsr/"
echo ""
echo "Next steps:"
echo "  1. Open the PharmCAT HTML report in a browser — it's the most actionable output"
echo "  2. Review ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt for a quick summary"
echo "  3. See docs/interpreting-results.md for help understanding your results"
