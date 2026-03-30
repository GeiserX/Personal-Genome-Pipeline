#!/usr/bin/env bash
# run-all.sh — Run the complete genomics analysis pipeline for one sample
# Usage: ./run-all.sh <sample_name> <sex: male|female>
#
# Assumes:
# - BAM already exists at $GENOMA_DIR/<sample>/aligned/<sample>_sorted.bam
# - VCF already exists at $GENOMA_DIR/<sample>/vcf/<sample>.vcf.gz
# - Manta SV VCF exists at $GENOMA_DIR/<sample>/manta/results/variants/diploidSV.vcf.gz
# - GRCh38 reference + ClinVar at $GENOMA_DIR/reference/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> <sex: male|female>}
SEX=${2:?Usage: $0 <sample_name> <sex: male|female>}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export GENOMA_DIR=${GENOMA_DIR:?Set GENOMA_DIR to your genomics data root}

echo "============================================"
echo "  Medical Genomics Pipeline"
echo "  Sample: ${SAMPLE}, Sex: ${SEX}"
echo "  Data: ${GENOMA_DIR}/${SAMPLE}/"
echo "============================================"
echo ""

# These can run in parallel (independent of each other)
echo "[1/9] AnnotSV — Structural variant annotation..."
bash "${SCRIPT_DIR}/05-annotsv.sh" "$SAMPLE" &
PID_ANNOTSV=$!

echo "[2/9] ClinVar screen — Known pathogenic variants..."
bash "${SCRIPT_DIR}/06-clinvar-screen.sh" "$SAMPLE" &
PID_CLINVAR=$!

echo "[3/9] ExpansionHunter — STR expansion screening..."
bash "${SCRIPT_DIR}/09-expansion-hunter.sh" "$SAMPLE" "$SEX" &
PID_EH=$!

echo "[4/9] TelomereHunter — Telomere length estimation..."
bash "${SCRIPT_DIR}/10-telomere-hunter.sh" "$SAMPLE" &
PID_TH=$!

echo "[5/9] ROH analysis — Consanguinity screening..."
bash "${SCRIPT_DIR}/11-roh-analysis.sh" "$SAMPLE" &
PID_ROH=$!

# Wait for quick jobs
wait $PID_CLINVAR $PID_ROH $PID_EH
echo ""
echo "[6/9] Imputation prep — MIS-ready VCFs..."
bash "${SCRIPT_DIR}/14-imputation-prep.sh" "$SAMPLE"

echo "[7/9] HLA typing — T1K Class I+II..."
bash "${SCRIPT_DIR}/08-hla-typing.sh" "$SAMPLE"

echo "[8/9] VEP annotation — Functional impact prediction..."
bash "${SCRIPT_DIR}/13-vep-annotation.sh" "$SAMPLE"

# Wait for remaining parallel jobs
wait $PID_ANNOTSV $PID_TH

echo ""
echo "============================================"
echo "  Pipeline complete for: ${SAMPLE}"
echo "  All results in: ${GENOMA_DIR}/${SAMPLE}/"
echo "============================================"
