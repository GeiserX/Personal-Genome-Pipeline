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

# Phase 0.5: fastp QC + trimming (if FASTQ exists but BAM doesn't)
BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
R1="${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz"
if [ ! -f "$BAM" ] && [ -f "$R1" ] && [ "${SKIP_TRIM:-false}" != "true" ]; then
  echo "[Phase 0.5] fastp QC + adapter trimming..."
  bash "${SCRIPT_DIR}/01b-fastp-qc.sh" "$SAMPLE"
elif [ "${SKIP_TRIM:-false}" = "true" ]; then
  echo "[Phase 0.5] fastp skipped (SKIP_TRIM=true)."
  export FASTQ_SUBDIR=fastq
fi
echo ""

# Phase 1: Alignment (if FASTQ exists but BAM doesn't)
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

# Phase 2b: Extra callers (optional, for benchmarking)
EXTRA_CALLERS=${EXTRA_CALLERS:-""}
if [ -n "$EXTRA_CALLERS" ]; then
  echo "[Phase 2b] Running extra variant callers: ${EXTRA_CALLERS}"
  PHASE2B_PIDS=()
  IFS=',' read -ra CALLERS <<< "$EXTRA_CALLERS"
  for CALLER in "${CALLERS[@]}"; do
    CALLER=$(echo "$CALLER" | tr -d ' ')
    case "$CALLER" in
      gatk)
        echo "  Starting GATK HaplotypeCaller..."
        bash "${SCRIPT_DIR}/03a-gatk-haplotypecaller.sh" "$SAMPLE" &
        PHASE2B_PIDS+=($!)
        ;;
      freebayes)
        echo "  Starting FreeBayes..."
        bash "${SCRIPT_DIR}/03b-freebayes.sh" "$SAMPLE" &
        PHASE2B_PIDS+=($!)
        ;;
      strelka2)
        echo "  Starting Strelka2..."
        echo "  NOTE: Strelka2 is using the default minimap2 BAM. For best SNP precision,"
        echo "        align with BWA-MEM2 first, then run: ALIGN_DIR=aligned_bwamem2 ./scripts/03c-strelka2-germline.sh $SAMPLE"
        bash "${SCRIPT_DIR}/03c-strelka2-germline.sh" "$SAMPLE" &
        PHASE2B_PIDS+=($!)
        ;;
      octopus)
        echo "  Starting Octopus..."
        bash "${SCRIPT_DIR}/03d-octopus.sh" "$SAMPLE" &
        PHASE2B_PIDS+=($!)
        ;;
      *)
        echo "  WARNING: Unknown caller '${CALLER}'. Skipping."
        ;;
    esac
  done
  PHASE2B_FAIL=0
  for PID in "${PHASE2B_PIDS[@]}"; do
    wait "$PID" 2>/dev/null || PHASE2B_FAIL=$((PHASE2B_FAIL + 1))
  done
  if [ "$PHASE2B_FAIL" -gt 0 ]; then
    echo "  WARNING: ${PHASE2B_FAIL} extra caller(s) failed."
  else
    echo "  Extra callers complete."
  fi
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

echo "  [A8] mosdepth coverage stats..."
bash "${SCRIPT_DIR}/16b-mosdepth.sh" "$SAMPLE" &
PID_MOSDEPTH=$!

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

# Wait for quick jobs, counting failures
PHASE3_FAIL=0
for PID in $PID_CLINVAR $PID_PHARMCAT $PID_ROH $PID_HAPLO $PID_INDEXCOV $PID_MOSDEPTH $PID_IMPUTATION $PID_HLA; do
  wait "$PID" 2>/dev/null || PHASE3_FAIL=$((PHASE3_FAIL + 1))
done
echo ""
if [ "$PHASE3_FAIL" -gt 0 ]; then
  echo "  WARNING: ${PHASE3_FAIL} quick analysis step(s) failed."
else
  echo "  Quick analyses complete."
fi

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

echo "  [C4] GRIDSS assembly-based SV calling..."
bash "${SCRIPT_DIR}/04b-gridss.sh" "$SAMPLE" &
PID_GRIDSS=$!

# Wait for Manta before running duphold and AnnotSV
PID_DUPHOLD=""
PID_ANNOTSV=""
if wait "$PID_MANTA" 2>/dev/null; then
  echo "  Manta complete. Running SV post-processing..."

  echo "  [B6] duphold SV quality annotation..."
  bash "${SCRIPT_DIR}/15-duphold.sh" "$SAMPLE" &
  PID_DUPHOLD=$!

  echo "  [B7] AnnotSV structural variant annotation..."
  bash "${SCRIPT_DIR}/05-annotsv.sh" "$SAMPLE" &
  PID_ANNOTSV=$!
else
  PHASE3_FAIL=$((PHASE3_FAIL + 1))
  echo "  WARNING: Manta failed — skipping duphold and AnnotSV."
fi

# Wait for remaining Phase 3 jobs
for PID in $PID_EH $PID_TH $PID_MTOOLBOX $PID_CPSR $PID_VEP $PID_CNVNATOR $PID_DELLY $PID_GRIDSS $PID_DUPHOLD $PID_ANNOTSV; do
  [ -z "$PID" ] && continue
  wait "$PID" 2>/dev/null || PHASE3_FAIL=$((PHASE3_FAIL + 1))
done

if [ "$PHASE3_FAIL" -gt 0 ]; then
  echo ""
  echo "  WARNING: ${PHASE3_FAIL} Phase 3 step(s) had errors. Check individual step output above."
fi

# Phase 4: Post-processing (uses outputs from Phase 3)
echo ""
echo "[Phase 4] Running post-processing steps..."

# Post-processing steps: log errors instead of swallowing them
POST_LOG="${GENOME_DIR}/${SAMPLE}/post_processing.log"
: > "$POST_LOG"

echo "  [D1] CYP2D6 star alleles (Cyrius) [experimental]..."
bash "${SCRIPT_DIR}/21-cyrius.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_CYRIUS=$!

echo "  [D2] SV consensus merge [experimental]..."
bash "${SCRIPT_DIR}/22-survivor-merge.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_SURVIVOR=$!

echo "  [D3] Clinical variant filter..."
bash "${SCRIPT_DIR}/23-clinical-filter.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_CLINICAL=$!

echo "  [D4] Polygenic Risk Scores [exploratory]..."
bash "${SCRIPT_DIR}/25-prs.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_PRS=$!

echo "  [D5] Ancestry PCA [experimental]..."
bash "${SCRIPT_DIR}/26-ancestry.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_ANCESTRY=$!

echo "  [D6] CPIC drug-gene recommendations..."
bash "${SCRIPT_DIR}/27-cpic-lookup.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_CPIC=$!

echo "  [D7] Somatic variant calling (Mutect2 tumor-only) [experimental]..."
bash "${SCRIPT_DIR}/29-mutect2-somatic.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 &
PID_SOMATIC=$!

PHASE4_FAIL=0
for PID in $PID_CYRIUS $PID_SURVIVOR $PID_CLINICAL $PID_PRS $PID_ANCESTRY $PID_CPIC $PID_SOMATIC; do
  wait "$PID" 2>/dev/null || PHASE4_FAIL=$((PHASE4_FAIL + 1))
done
if [ "$PHASE4_FAIL" -gt 0 ]; then
  echo "  WARNING: ${PHASE4_FAIL} post-processing step(s) had errors."
  echo "  See: ${POST_LOG}"
else
  echo "  Post-processing complete."
fi

# Phase 4b: Benchmarking (optional)
BENCHMARK=${BENCHMARK:-false}
if [ "$BENCHMARK" = "true" ] || [ "$BENCHMARK" = "1" ]; then
  # Count available caller VCFs (need at least 2 for pairwise comparison)
  CALLER_COUNT=0
  for d in vcf vcf_gatk vcf_freebayes vcf_octopus; do
    [ -f "${GENOME_DIR}/${SAMPLE}/${d}/${SAMPLE}.vcf.gz" ] && CALLER_COUNT=$((CALLER_COUNT + 1))
  done
  # Strelka2 writes to a different path
  [ -f "${GENOME_DIR}/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz" ] && CALLER_COUNT=$((CALLER_COUNT + 1))
  if [ "$CALLER_COUNT" -ge 2 ]; then
    echo ""
    echo "  [D8] Variant caller benchmarking (${CALLER_COUNT} caller VCFs found)..."
    bash "${SCRIPT_DIR}/benchmark-variants.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 || { echo "  WARNING: Benchmarking failed. See ${POST_LOG}"; PHASE4_FAIL=$((PHASE4_FAIL + 1)); }
  else
    echo ""
    echo "  [D8] Skipping benchmarking: only ${CALLER_COUNT} caller VCF(s) found (need 2+)."
    echo "  Set EXTRA_CALLERS=gatk,freebayes,strelka2 to run alternative callers in Phase 2b."
  fi
fi

REPORT_FAIL=0

echo "  [D9] HTML summary report..."
bash "${SCRIPT_DIR}/24-html-report.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 || { echo "  WARNING: HTML report generation failed. See ${POST_LOG}"; REPORT_FAIL=$((REPORT_FAIL + 1)); }

echo "  [D10] MultiQC aggregated QC report..."
bash "${SCRIPT_DIR}/28-multiqc.sh" "$SAMPLE" >> "$POST_LOG" 2>&1 || { echo "  WARNING: MultiQC report generation failed. See ${POST_LOG}"; REPORT_FAIL=$((REPORT_FAIL + 1)); }

# Generate summary report
echo ""
echo "[Report] Generating summary report..."
bash "${SCRIPT_DIR}/generate-report.sh" "$SAMPLE" 2>/dev/null || { echo "  (report generation had warnings — check output manually)"; REPORT_FAIL=$((REPORT_FAIL + 1)); }

PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

TOTAL_FAIL=$((${PHASE2B_FAIL:-0} + PHASE3_FAIL + PHASE4_FAIL + REPORT_FAIL))

echo ""
echo "============================================"
if [ "$TOTAL_FAIL" -gt 0 ]; then
  echo "  Pipeline finished with errors for: ${SAMPLE}"
  echo "  ${TOTAL_FAIL} step(s) failed (Phase 2b: ${PHASE2B_FAIL:-0}, Phase 3: ${PHASE3_FAIL}, Phase 4: ${PHASE4_FAIL}, Reports: ${REPORT_FAIL})"
else
  echo "  Pipeline complete for: ${SAMPLE}"
fi
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
echo "  PharmCAT:       ${GENOME_DIR}/${SAMPLE}/vcf/ (PharmCAT reports alongside VCF)"
echo "  CYP2D6:         ${GENOME_DIR}/${SAMPLE}/cyrius/"
echo "  CPIC drugs:     ${GENOME_DIR}/${SAMPLE}/cpic/"
echo "  Clinical VCF:   ${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz"
echo "  SV consensus:   ${GENOME_DIR}/${SAMPLE}/sv_merged/"
echo "  PRS scores:     ${GENOME_DIR}/${SAMPLE}/prs/"
echo "  Ancestry PCA:   ${GENOME_DIR}/${SAMPLE}/ancestry/"
echo "  Somatic calls:  ${GENOME_DIR}/${SAMPLE}/somatic/"
echo "  CPSR report:    ${GENOME_DIR}/${SAMPLE}/cpsr/"
echo ""
echo "Next steps:"
echo "  1. Open the PharmCAT HTML report in a browser — it's the most actionable output"
echo "  2. Review ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt for a quick summary"
echo "  3. See docs/interpreting-results.md for help understanding your results"

exit "$TOTAL_FAIL"
