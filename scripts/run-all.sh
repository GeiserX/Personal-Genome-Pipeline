#!/usr/bin/env bash
# run-all.sh — Run the complete genomics analysis pipeline for one sample
# Usage: ./run-all.sh <sample_name> <sex: male|female>
#
# Opt-in steps (off by default):
#   GRIDSS=true      step 4b, needs a classic BWA index and ~32 GB RAM
#   ANCESTRY=true    step 26, downloads ~1 GB and only counts shared SNPs on one sample
#   IMPUTATION=true  step 14, per-chromosome VCFs for an imputation server upload
#   SOMATIC=true     step 29, tumor-only Mutect2 (high false-positive rate)
#   EXTRA_CALLERS=gatk,freebayes,strelka2,octopus   alternative callers
#   BENCHMARK=true   compare caller VCFs
#
# Every step writes its log to $GENOME_DIR/<sample>/logs/. The run ends with a
# table of each step as ok, skipped or failed, and exits 1 only when a step
# failed. Steps whose data is not installed are reported as skipped: VEP
# (step 13, and with it steps 30, 23 and 31), CNVpytor, CPSR, pypgx, AnnotSV.
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
case "$SEX" in
  male|female) ;;
  *)
    echo "Usage: $0 <sample_name> <sex: male|female>" >&2
    echo "ERROR: sex must be 'male' or 'female', got '${SEX}'." >&2
    exit 2
    ;;
esac
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

export GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

# Concurrency control — limit parallel Docker containers to prevent host oversubscription
# Default: half the CPU count, clamped to [4, 12]
_detect_max_jobs() {
  local cpus
  cpus=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)
  local max=$(( cpus / 2 ))
  [ "$max" -lt 4 ] && max=4
  [ "$max" -gt 12 ] && max=12
  echo "$max"
}
MAX_JOBS=${MAX_JOBS:-$(_detect_max_jobs)}

_throttle() {
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
    wait -n 2>/dev/null || sleep 2
  done
}

# Step bookkeeping: one entry per step, reported in the final table.
LOG_DIR="${GENOME_DIR}/${SAMPLE}/logs"
mkdir -p "$LOG_DIR"
STEP_NAMES=()
STEP_LOGS=()
STEP_PIDS=()
STEP_RESULTS=()
LAST_STEP=0

_record() {  # <name> <log> <pid> <result>
  STEP_NAMES+=("$1")
  STEP_LOGS+=("$2")
  STEP_PIDS+=("$3")
  STEP_RESULTS+=("$4")
  LAST_STEP=$(( ${#STEP_NAMES[@]} - 1 ))
}

# Start a step in the background, its output in logs/<log>.log
_launch() {  # <name> <log> <script> [args...]
  local name="$1" log="${LOG_DIR}/$2.log"
  shift 2
  echo "  [${name}] started (log: ${log})"
  _throttle
  bash "$@" > "$log" 2>&1 &
  _record "$name" "$log" "$!" "running"
}

# Run a step in the foreground, its output in logs/<log>.log
_run() {  # <name> <log> <script> [args...]
  local name="$1" log="${LOG_DIR}/$2.log"
  shift 2
  echo "  [${name}] running (log: ${log})"
  if bash "$@" > "$log" 2>&1; then
    _record "$name" "$log" "" "ok"
  else
    _record "$name" "$log" "" "failed"
    echo "  [${name}] FAILED. See ${log}"
  fi
}

_skip() {  # <name> <reason>
  echo "  [${1}] skipped (${2})"
  _record "$1" "-" "" "skipped (${2})"
}

# Wait for one launched step; returns 0 if it succeeded
_wait_step() {  # <index>
  local i="$1"
  if [ "${STEP_RESULTS[$i]}" = "running" ]; then
    if wait "${STEP_PIDS[$i]}"; then
      STEP_RESULTS[i]="ok"
    else
      STEP_RESULTS[i]="failed"
      echo "  [${STEP_NAMES[$i]}] FAILED. See ${STEP_LOGS[$i]}"
    fi
  fi
  [ "${STEP_RESULTS[$i]}" = "ok" ]
}

_wait_all() {
  local i
  [ "${#STEP_NAMES[@]}" -gt 0 ] || return 0
  for i in "${!STEP_NAMES[@]}"; do
    _wait_step "$i" || true
  done
}

_enabled() {  # <VAR value>: true for "true" or "1"
  [ "$1" = "true" ] || [ "$1" = "1" ]
}

PIPELINE_START=$(date +%s)

echo "============================================"
echo "  Personal Genome Pipeline — Full Analysis"
echo "  Sample: ${SAMPLE}, Sex: ${SEX}"
echo "  Data: ${GENOME_DIR}/${SAMPLE}/"
echo "  Max parallel jobs: ${MAX_JOBS} (override: MAX_JOBS=N)"
echo "  Started: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""

# Pre-flight check — abort on failure unless explicitly skipped
echo "[Pre-flight] Validating setup..."
if [ "${SKIP_VALIDATION:-false}" = "true" ]; then
  echo "  Skipping validation (SKIP_VALIDATION=true)."
elif ! "${SCRIPT_DIR}/validate-setup.sh" "${SAMPLE}"; then
  echo ""
  echo "ERROR: Setup validation failed. Fix the issues above before running the pipeline."
  echo "  To bypass: SKIP_VALIDATION=true ./scripts/run-all.sh $SAMPLE $SEX"
  exit 1
fi

# What produces this run's outputs (images and their digests, database
# releases, the pipeline commit), and when the run started. The reports read
# logs/run_status.tsv to tell this run's results from an earlier run's: a
# step that is skipped or fails here leaves its older output on disk.
GENOME_DIR="$GENOME_DIR" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" \
  || echo "WARNING: could not write ${GENOME_DIR}/${SAMPLE}/run_manifest.tsv; the reports write one later."
RUN_STATUS="${LOG_DIR}/run_status.tsv"
{
  printf '# run-all.sh: when this run started and how each step ended (bin/collect_summary.py reads it)\n'
  printf 'meta\tstarted_epoch\t%s\n' "$PIPELINE_START"
  printf 'meta\tstarted_utc\t%s\n' "$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  printf 'meta\tdeclared_sex\t%s\n' "$SEX"
} > "$RUN_STATUS"
echo ""

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

# Phase 1: Alignment, unless a BAM with its index that passes samtools
# quickcheck is there (step 02 renames its BAM into place only after both).
if [ -f "$BAM" ] && [ -f "${BAM}.bai" ] \
   && run_in "$SAMTOOLS_IMAGE" samtools quickcheck "$(cpath "$BAM")"; then
  echo "[Phase 1] BAM and its index already exist, skipping alignment."
else
  echo "[Phase 1] Alignment — FASTQ to sorted, duplicate-marked BAM..."
  bash "${SCRIPT_DIR}/02-alignment.sh" "$SAMPLE"
fi
echo ""

# Phase 2: Variant calling, unless the VCF and its index are there. The sex
# makes DeepVariant call chrX and chrY haploid for a male sample.
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
if [ -f "$VCF" ] && [ -f "${VCF}.tbi" ]; then
  echo "[Phase 2] VCF and its index already exist, skipping variant calling."
  if [ ! -f "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz" ]; then
    echo "  NOTE: no gVCF next to it (an older run). PharmCAT and PRS read only the variant sites;"
    echo "  remove the VCF to call again and get ${SAMPLE}.g.vcf.gz as well."
  fi
else
  echo "[Phase 2] Variant calling — DeepVariant..."
  bash "${SCRIPT_DIR}/03-deepvariant.sh" "$SAMPLE" "$SEX"
fi

# Phase 2b: Extra callers (optional, for benchmarking)
EXTRA_CALLERS=${EXTRA_CALLERS:-""}
if [ -n "$EXTRA_CALLERS" ]; then
  echo "[Phase 2b] Running extra variant callers: ${EXTRA_CALLERS}"
  IFS=',' read -ra CALLERS <<< "$EXTRA_CALLERS"
  for CALLER in "${CALLERS[@]}"; do
    CALLER=$(echo "$CALLER" | tr -d ' ')
    case "$CALLER" in
      gatk)
        _launch "03a GATK HaplotypeCaller" 03a_gatk "${SCRIPT_DIR}/03a-gatk-haplotypecaller.sh" "$SAMPLE"
        ;;
      freebayes)
        _launch "03b FreeBayes" 03b_freebayes "${SCRIPT_DIR}/03b-freebayes.sh" "$SAMPLE"
        ;;
      strelka2)
        echo "  NOTE: Strelka2 is using the default minimap2 BAM. For best SNP precision,"
        echo "        align with BWA-MEM2 first, then run: ALIGN_DIR=aligned_bwamem2 ./scripts/03c-strelka2-germline.sh $SAMPLE"
        _launch "03c Strelka2" 03c_strelka2 "${SCRIPT_DIR}/03c-strelka2-germline.sh" "$SAMPLE"
        ;;
      octopus)
        _launch "03d Octopus" 03d_octopus "${SCRIPT_DIR}/03d-octopus.sh" "$SAMPLE"
        ;;
      *)
        echo "  WARNING: Unknown caller '${CALLER}'. Skipping."
        ;;
    esac
  done
  _wait_all
  echo "  Extra callers finished."
fi
echo ""

# Phase 3: Parallel analyses (all independent after BAM + VCF exist)
echo "[Phase 3] Running parallel analyses..."
echo ""

# --- Group A: Quick jobs (minutes each) ---
echo "  Starting quick analyses..."
QUICK_FIRST=${#STEP_NAMES[@]}
_launch "06 ClinVar screen" 06_clinvar "${SCRIPT_DIR}/06-clinvar-screen.sh" "$SAMPLE"
_launch "07 PharmCAT" 07_pharmcat "${SCRIPT_DIR}/07-pharmacogenomics.sh" "$SAMPLE"
_launch "11 ROH" 11_roh "${SCRIPT_DIR}/11-roh-analysis.sh" "$SAMPLE"
_launch "12 Mito haplogroup" 12_haplogroup "${SCRIPT_DIR}/12-mito-haplogroup.sh" "$SAMPLE"
_launch "16 indexcov" 16_indexcov "${SCRIPT_DIR}/16-indexcov.sh" "$SAMPLE" "$SEX"
_launch "16b mosdepth" 16b_mosdepth "${SCRIPT_DIR}/16b-mosdepth.sh" "$SAMPLE"
if _enabled "${IMPUTATION:-false}"; then
  _launch "14 Imputation prep" 14_imputation "${SCRIPT_DIR}/14-imputation-prep.sh" "$SAMPLE"
else
  echo "  [14 Imputation prep] off (opt-in: IMPUTATION=true)"
fi
_launch "08 HLA typing (T1K)" 08_hla "${SCRIPT_DIR}/08-hla-typing.sh" "$SAMPLE"
QUICK_LAST=$LAST_STEP

# --- Group B: Medium jobs (10-60 minutes each) ---
echo "  Starting medium analyses..."
_launch "04 Manta" 04_manta "${SCRIPT_DIR}/04-manta.sh" "$SAMPLE"
IDX_MANTA=$LAST_STEP
_launch "09 ExpansionHunter" 09_expansionhunter "${SCRIPT_DIR}/09-expansion-hunter.sh" "$SAMPLE" "$SEX"
IDX_EH=$LAST_STEP
_launch "10 TelomereHunter" 10_telomerehunter "${SCRIPT_DIR}/10-telomere-hunter.sh" "$SAMPLE"
_launch "20 Mito variants (Mutect2)" 20_mito "${SCRIPT_DIR}/20-mtoolbox.sh" "$SAMPLE"

# CPSR and pypgx need optional data; without it they are skipped, not failed
if [ -f "${GENOME_DIR}/vep_cache/homo_sapiens/${PCGR_VEP_CACHE_RELEASE}_GRCh38/info.txt" ] && [ -d "${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}/data" ]; then
  _launch "17 CPSR" 17_cpsr "${SCRIPT_DIR}/17-cpsr.sh" "$SAMPLE"
else
  _skip "17 CPSR" "data not installed: VEP ${PCGR_VEP_CACHE_RELEASE} cache and PCGR bundle, see docs/17-cpsr.md"
fi
if [ -d "${GENOME_DIR}/reference/pypgx-bundle" ]; then
  _launch "32 pypgx" 32_pypgx "${SCRIPT_DIR}/32-pypgx.sh" "$SAMPLE"
else
  _skip "32 pypgx" "data not installed: reference/pypgx-bundle, see docs/32-pypgx.md"
fi

# Wait for the quick jobs before starting the heavy ones
for i in $(seq "$QUICK_FIRST" "$QUICK_LAST"); do
  _wait_step "$i" || true
done
echo "  Quick analyses finished."

# --- Group C: Heavy jobs (2-4 hours each) ---
echo "  Starting heavy analyses..."
# VEP and CNVpytor need data setup.sh does not download; without it they are
# skipped, not failed (and step 13 does not start a 26 GB download mid-run)
if [ -f "${GENOME_DIR}/vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38/info.txt" ]; then
  _launch "13 VEP" 13_vep "${SCRIPT_DIR}/13-vep-annotation.sh" "$SAMPLE"
else
  _skip "13 VEP" "data not installed: vep_cache/homo_sapiens/${VEP_CACHE_RELEASE}_GRCh38, see docs/13-vep-annotation.md"
fi
IDX_VEP=$LAST_STEP
if [ -s "${GENOME_DIR}/reference/cnvpytor/gc_hg38.pytor" ]; then
  _launch "18 CNVpytor" 18_cnvpytor "${SCRIPT_DIR}/18-cnvpytor.sh" "$SAMPLE"
else
  _skip "18 CNVpytor" "data not installed: reference/cnvpytor, see docs/18-cnvpytor.md"
fi
_launch "19 Delly" 19_delly "${SCRIPT_DIR}/19-delly.sh" "$SAMPLE"
if _enabled "${GRIDSS:-false}"; then
  _launch "04b GRIDSS" 04b_gridss "${SCRIPT_DIR}/04b-gridss.sh" "$SAMPLE"
else
  echo "  [04b GRIDSS] off (opt-in: GRIDSS=true; needs a classic BWA index, see docs/04b-gridss.md)"
fi

# duphold and AnnotSV need Manta's VCF
if _wait_step "$IDX_MANTA"; then
  _launch "15 duphold" 15_duphold "${SCRIPT_DIR}/15-duphold.sh" "$SAMPLE"
  if [ -d "${GENOME_DIR}/annotsv_annotations/Annotations_Human/Genes/GRCh38" ]; then
    _launch "05 AnnotSV" 05_annotsv "${SCRIPT_DIR}/05-annotsv.sh" "$SAMPLE"
  else
    _skip "05 AnnotSV" "data not installed: annotsv_annotations, run setup.sh"
  fi
else
  _skip "15 duphold" "Manta failed"
  _skip "05 AnnotSV" "Manta failed"
fi

# Stranger needs ExpansionHunter's output
if _wait_step "$IDX_EH"; then
  _launch "09b Stranger" 09b_stranger "${SCRIPT_DIR}/09b-stranger.sh" "$SAMPLE"
else
  _skip "09b Stranger" "ExpansionHunter failed"
fi

_wait_all
echo "  Phase 3 finished."

# Phase 4: Post-processing (uses outputs from Phase 3)
echo ""
echo "[Phase 4] Running post-processing steps..."

# Steps 30, 23 and 31 read step 13's VEP output: run them only when step 13
# succeeded in this run, so they never work on an old or missing file
VEP_RESULT="${STEP_RESULTS[$IDX_VEP]}"
VEP_RESULT="${VEP_RESULT%% *}"

# vcfanno annotation enrichment (must complete before clinical filter)
# Adds CADD, SpliceAI, REVEL, AlphaMissense scores to VEP VCF
if [ "$VEP_RESULT" = "ok" ]; then
  _run "30 vcfanno" 30_vcfanno "${SCRIPT_DIR}/30-vcfanno.sh" "$SAMPLE"
else
  _skip "30 vcfanno" "needs VEP, step 13 ${VEP_RESULT}"
fi

_launch "21 Cyrius CYP2D6 [experimental]" 21_cyrius "${SCRIPT_DIR}/21-cyrius.sh" "$SAMPLE"
_launch "22 SV consensus merge [experimental]" 22_survivor "${SCRIPT_DIR}/22-survivor-merge.sh" "$SAMPLE"
if [ "$VEP_RESULT" = "ok" ]; then
  _launch "23 Clinical filter" 23_clinical "${SCRIPT_DIR}/23-clinical-filter.sh" "$SAMPLE"
else
  _skip "23 Clinical filter" "needs VEP, step 13 ${VEP_RESULT}"
fi
_launch "25 PRS [exploratory]" 25_prs "${SCRIPT_DIR}/25-prs.sh" "$SAMPLE"
if _enabled "${ANCESTRY:-false}"; then
  _launch "26 Ancestry SNP intersection [experimental]" 26_ancestry "${SCRIPT_DIR}/26-ancestry.sh" "$SAMPLE"
else
  echo "  [26 Ancestry] off (opt-in: ANCESTRY=true; on one sample it only counts shared SNPs)"
fi
_launch "27 CPIC lookup" 27_cpic "${SCRIPT_DIR}/27-cpic-lookup.sh" "$SAMPLE"
if [ "$VEP_RESULT" = "ok" ]; then
  _launch "31 slivar" 31_slivar "${SCRIPT_DIR}/31-slivar.sh" "$SAMPLE"
else
  _skip "31 slivar" "needs VEP, step 13 ${VEP_RESULT}"
fi

# Mutect2 somatic (tumor-only) is opt-in due to high false-positive rate.
if _enabled "${SOMATIC:-false}"; then
  _launch "29 Somatic (Mutect2 tumor-only) [experimental]" 29_somatic "${SCRIPT_DIR}/29-mutect2-somatic.sh" "$SAMPLE"
else
  echo "  [29 Somatic] off (opt-in: SOMATIC=true; high false-positive rate)"
fi

_wait_all
echo "  Post-processing finished."

# Phase 4b: Benchmarking (optional)
if _enabled "${BENCHMARK:-false}"; then
  # Count available caller VCFs (need at least 2 for pairwise comparison)
  CALLER_COUNT=0
  for d in vcf vcf_gatk vcf_freebayes vcf_octopus; do
    [ -f "${GENOME_DIR}/${SAMPLE}/${d}/${SAMPLE}.vcf.gz" ] && CALLER_COUNT=$((CALLER_COUNT + 1))
  done
  # Strelka2 writes to a different path
  [ -f "${GENOME_DIR}/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz" ] && CALLER_COUNT=$((CALLER_COUNT + 1))
  if [ "$CALLER_COUNT" -ge 2 ]; then
    echo ""
    echo "  Variant caller benchmarking (${CALLER_COUNT} caller VCFs found)..."
    _run "Benchmark callers" benchmark "${SCRIPT_DIR}/benchmark-variants.sh" "$SAMPLE"
  else
    _skip "Benchmark callers" "only ${CALLER_COUNT} caller VCF found, need 2+; set EXTRA_CALLERS=gatk,freebayes,strelka2"
  fi
fi

# Every step's result in this run, for the reports (written before they run).
for i in "${!STEP_NAMES[@]}"; do
  printf 'step\t%s\t%s\n' "${STEP_NAMES[$i]%% *}" "${STEP_RESULTS[$i]}" >> "$RUN_STATUS"
done

# Rewrite the manifest now that every step has run: an image a step pulled
# after the start (possible with SKIP_VALIDATION=true) gets its digest.
GENOME_DIR="$GENOME_DIR" bash "${PGP_ROOT}/bin/write_manifest.sh" "$SAMPLE" run-all.sh "$SEX" \
  || echo "WARNING: could not refresh ${GENOME_DIR}/${SAMPLE}/run_manifest.tsv; it keeps the digests from the start of the run."

_run "24 HTML report" 24_html_report "${SCRIPT_DIR}/24-html-report.sh" "$SAMPLE"
_run "28 MultiQC" 28_multiqc "${SCRIPT_DIR}/28-multiqc.sh" "$SAMPLE"
_run "Summary report" generate_report "${SCRIPT_DIR}/generate-report.sh" "$SAMPLE"

# Aggregate the post-processing logs into one file for easy review
POST_LOG="${GENOME_DIR}/${SAMPLE}/post_processing.log"
: > "$POST_LOG"
# Only logs of steps this run started: a log left by an earlier run with
# other options must not end up in it.
for logf in "${STEP_LOGS[@]}"; do
  case "${logf##*/}" in
    2[0-9]_*.log|3[0-9]_*.log|benchmark.log|generate_report.log) ;;
    *) continue ;;
  esac
  [ -f "$logf" ] || continue
  echo "=== $(basename "$logf") ===" >> "$POST_LOG"
  cat "$logf" >> "$POST_LOG"
  echo "" >> "$POST_LOG"
done

PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
HOURS=$(( ELAPSED / 3600 ))
MINUTES=$(( (ELAPSED % 3600) / 60 ))

# Final table: every step, its result and its log
OK_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0
echo ""
echo "============================================"
printf '  %-45s %-10s %s\n' "Step" "Result" "Log"
for i in "${!STEP_NAMES[@]}"; do
  RESULT="${STEP_RESULTS[$i]}"
  case "$RESULT" in
    ok) OK_COUNT=$((OK_COUNT + 1)) ;;
    failed) FAIL_COUNT=$((FAIL_COUNT + 1)) ;;
    skipped*) SKIP_COUNT=$((SKIP_COUNT + 1)) ;;
  esac
  if [ "${RESULT%% *}" = "skipped" ]; then
    printf '  %-45s %-10s %s\n' "${STEP_NAMES[$i]}" "skipped" "${RESULT#skipped }"
  else
    printf '  %-45s %-10s %s\n' "${STEP_NAMES[$i]}" "$RESULT" "${STEP_LOGS[$i]}"
  fi
done
echo "============================================"
if [ "$FAIL_COUNT" -gt 0 ]; then
  echo "  Pipeline finished with errors for: ${SAMPLE}"
  echo "  ${OK_COUNT} ok, ${SKIP_COUNT} skipped, ${FAIL_COUNT} failed"
  echo "  Failed:"
  for i in "${!STEP_NAMES[@]}"; do
    if [ "${STEP_RESULTS[$i]}" = "failed" ]; then
      echo "    ${STEP_NAMES[$i]}: ${STEP_LOGS[$i]}"
    fi
  done
else
  echo "  Pipeline complete for: ${SAMPLE}"
  echo "  ${OK_COUNT} ok, ${SKIP_COUNT} skipped, 0 failed"
fi
echo "  All results in: ${GENOME_DIR}/${SAMPLE}/"
echo "  Step logs in:   ${LOG_DIR}/"
echo "  Total runtime: ${HOURS}h ${MINUTES}m"
echo "  Finished: $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================"
echo ""
echo "Key outputs:"
echo "  HTML Report:    ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html"
echo "  Text Report:    ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt"
echo "  Summary JSON:   ${GENOME_DIR}/${SAMPLE}/summary.json (both reports are rendered from it)"
echo "  Run manifest:   ${GENOME_DIR}/${SAMPLE}/run_manifest.tsv"
echo "  VCF:            ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
echo "  ClinVar hits:   ${GENOME_DIR}/${SAMPLE}/clinvar/"
echo "  PharmCAT:       ${GENOME_DIR}/${SAMPLE}/vcf/ (PharmCAT reports alongside VCF)"
echo "  CYP2D6:         ${GENOME_DIR}/${SAMPLE}/cyrius/"
echo "  CPIC drugs:     ${GENOME_DIR}/${SAMPLE}/cpic/"
echo "  Clinical VCF:   ${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz"
echo "  SV consensus:   ${GENOME_DIR}/${SAMPLE}/sv_merged/"
echo "  PRS scores:     ${GENOME_DIR}/${SAMPLE}/prs/"
echo "  pypgx PGx:      ${GENOME_DIR}/${SAMPLE}/pypgx/"
echo "  Slivar:         ${GENOME_DIR}/${SAMPLE}/slivar/"
echo "  CPSR report:    ${GENOME_DIR}/${SAMPLE}/cpsr/"
echo ""
echo "Next steps:"
echo "  1. Open the PharmCAT HTML report in a browser — it's the most actionable output"
echo "  2. Review ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.txt for a quick summary"
echo "  3. See docs/interpreting-results.md for help understanding your results"

exit "$(( FAIL_COUNT > 0 ? 1 : 0 ))"
