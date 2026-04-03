#!/usr/bin/env bash
# benchmark-variants.sh — Compare variant calls across multiple callers
# Supports pairwise concordance (bcftools isec) and truth set benchmarking (hap.py)
# Input: VCFs from vcf/, vcf_gatk/, vcf_freebayes/ directories
# Output: comparison tables in $GENOME_DIR/<sample>/benchmark/
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [--truth <vcf> --regions <bed>]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
shift

# --- Parse optional flags ---
TRUTH_VCF=""
REGIONS_BED=""
while [ $# -gt 0 ]; do
  case "$1" in
    --truth)
      TRUTH_VCF="${2:?--truth requires a VCF path}"
      shift 2
      ;;
    --regions)
      REGIONS_BED="${2:?--regions requires a BED path}"
      shift 2
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      echo "Usage: $0 <sample_name> [--truth <vcf> --regions <bed>]" >&2
      exit 1
      ;;
  esac
done

BENCHMARK_DIR="${GENOME_DIR}/${SAMPLE}/benchmark"
SUMMARY="${BENCHMARK_DIR}/summary.txt"
TSV="${BENCHMARK_DIR}/comparison.tsv"
REF="${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta"

# --- Discover available caller VCFs ---
declare -a CALLER_NAMES=()
declare -a CALLER_VCFS=()

DV_VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
GATK_VCF="${GENOME_DIR}/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
FB_VCF="${GENOME_DIR}/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"
SK_VCF="${GENOME_DIR}/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz"

if [ -f "$DV_VCF" ]; then
  CALLER_NAMES+=("DeepVariant")
  CALLER_VCFS+=("$DV_VCF")
fi
if [ -f "$GATK_VCF" ]; then
  CALLER_NAMES+=("GATK")
  CALLER_VCFS+=("$GATK_VCF")
fi
if [ -f "$FB_VCF" ]; then
  CALLER_NAMES+=("FreeBayes")
  CALLER_VCFS+=("$FB_VCF")
fi
if [ -f "$SK_VCF" ]; then
  CALLER_NAMES+=("Strelka2")
  CALLER_VCFS+=("$SK_VCF")
fi

NUM_CALLERS=${#CALLER_VCFS[@]}

echo "=== Variant Caller Benchmarking: ${SAMPLE} ==="
echo "Found ${NUM_CALLERS} caller VCF(s):"
for i in $(seq 0 $((NUM_CALLERS - 1))); do
  echo "  - ${CALLER_NAMES[$i]}: ${CALLER_VCFS[$i]}"
done

if [ "$NUM_CALLERS" -lt 2 ] && [ -z "$TRUTH_VCF" ]; then
  echo "" >&2
  echo "ERROR: Need at least 2 caller VCFs for pairwise comparison." >&2
  echo "Found VCFs in these locations:" >&2
  echo "  vcf/           (DeepVariant) — run scripts/03-deepvariant.sh" >&2
  echo "  vcf_gatk/      (GATK)       — run GATK HaplotypeCaller" >&2
  echo "  vcf_freebayes/ (FreeBayes)  — run FreeBayes" >&2
  echo "" >&2
  echo "Or use --truth <vcf> --regions <bed> for single-caller truth set benchmarking." >&2
  exit 1
fi

if [ -n "$TRUTH_VCF" ] && [ "$NUM_CALLERS" -lt 1 ]; then
  echo "ERROR: Need at least 1 caller VCF for truth set benchmarking." >&2
  exit 1
fi

# Validate index files exist for all caller VCFs
for i in $(seq 0 $((NUM_CALLERS - 1))); do
  vcf="${CALLER_VCFS[$i]}"
  if [ ! -f "${vcf}.tbi" ]; then
    echo "ERROR: Index not found: ${vcf}.tbi" >&2
    echo "Run: bcftools index -t ${vcf}" >&2
    exit 1
  fi
done

# Validate truth mode inputs
if [ -n "$TRUTH_VCF" ]; then
  for f in "$TRUTH_VCF" "$REF"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: File not found: ${f}" >&2
      exit 1
    fi
  done
  if [ -n "$REGIONS_BED" ] && [ ! -f "$REGIONS_BED" ]; then
    echo "ERROR: Regions BED not found: ${REGIONS_BED}" >&2
    exit 1
  fi
fi

mkdir -p "$BENCHMARK_DIR"

# =============================================================================
# TRUTH SET MODE (hap.py)
# =============================================================================
if [ -n "$TRUTH_VCF" ]; then
  echo ""
  echo "=== Truth Set Benchmarking (hap.py) ==="
  echo "Truth VCF: ${TRUTH_VCF}"
  [ -n "$REGIONS_BED" ] && echo "Confident regions: ${REGIONS_BED}"
  echo ""
  echo "WARNING: Truth set benchmarking is only meaningful when the query VCF was"
  echo "generated from the SAME biological sample as the truth set. If the truth set"
  echo "is HG002, the query VCF must also come from HG002 sequencing data. Running"
  echo "hap.py against a non-matching sample produces meaningless precision/recall."
  echo ""

  # Build the TSV header
  printf "Caller\tSNP_TP\tSNP_FP\tSNP_FN\tSNP_Precision\tSNP_Recall\tSNP_F1\tINDEL_TP\tINDEL_FP\tINDEL_FN\tINDEL_Precision\tINDEL_Recall\tINDEL_F1\n" > "$TSV"

  {
    echo "================================================================================"
    echo "  VARIANT CALLER TRUTH SET BENCHMARK"
    echo "  Sample: ${SAMPLE}"
    echo "  Truth:  $(basename "$TRUTH_VCF")"
    echo "  Generated: $(date -u '+%Y-%m-%d %H:%M UTC')"
    echo "================================================================================"
    echo ""
  } > "$SUMMARY"

  TRUTH_CONTAINER_PATH="${TRUTH_VCF/#$GENOME_DIR//genome}"
  REGIONS_FLAG=""
  if [ -n "$REGIONS_BED" ]; then
    REGIONS_CONTAINER_PATH="${REGIONS_BED/#$GENOME_DIR//genome}"
    REGIONS_FLAG="-f ${REGIONS_CONTAINER_PATH}"
  fi

  for i in $(seq 0 $((NUM_CALLERS - 1))); do
    CALLER="${CALLER_NAMES[$i]}"
    VCF="${CALLER_VCFS[$i]}"
    VCF_CONTAINER="/genome/${SAMPLE}/${VCF//${GENOME_DIR}\/${SAMPLE}\//}"
    VCF_CONTAINER="${VCF/#$GENOME_DIR//genome}"
    PREFIX="/genome/${SAMPLE}/benchmark/happy_${CALLER}"

    echo ""
    echo "=== Running hap.py: ${CALLER} vs truth ==="

    # shellcheck disable=SC2086
    docker run --rm \
      --cpus 4 --memory 8g \
      --user root \
      -v "${GENOME_DIR}:/genome" \
      jmcdani20/hap.py:v0.3.12 \
      /opt/hap.py/bin/hap.py \
        "${TRUTH_CONTAINER_PATH}" \
        "${VCF_CONTAINER}" \
        -r /genome/reference/Homo_sapiens_assembly38.fasta \
        ${REGIONS_FLAG} \
        -o "${PREFIX}" \
        --engine=vcfeval

    # Parse hap.py summary.csv
    HAPPY_CSV="${BENCHMARK_DIR}/happy_${CALLER}.summary.csv"
    if [ ! -f "$HAPPY_CSV" ]; then
      echo "WARNING: hap.py summary not found for ${CALLER}, skipping" >&2
      continue
    fi

    # Extract SNP and INDEL metrics from summary.csv
    # Columns: Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score
    SNP_LINE=$(awk -F',' '$1=="SNP" && $2=="PASS"' "$HAPPY_CSV" || true)
    INDEL_LINE=$(awk -F',' '$1=="Indel" && $2=="PASS"' "$HAPPY_CSV" || true)

    SNP_TP=$(echo "$SNP_LINE" | awk -F',' '{print $4}')
    SNP_FP=$(echo "$SNP_LINE" | awk -F',' '{print $7}')
    SNP_FN=$(echo "$SNP_LINE" | awk -F',' '{print $5}')
    SNP_PREC=$(echo "$SNP_LINE" | awk -F',' '{print $11}')
    SNP_RECALL=$(echo "$SNP_LINE" | awk -F',' '{print $10}')
    SNP_F1=$(echo "$SNP_LINE" | awk -F',' '{print $13}')

    INDEL_TP=$(echo "$INDEL_LINE" | awk -F',' '{print $4}')
    INDEL_FP=$(echo "$INDEL_LINE" | awk -F',' '{print $7}')
    INDEL_FN=$(echo "$INDEL_LINE" | awk -F',' '{print $5}')
    INDEL_PREC=$(echo "$INDEL_LINE" | awk -F',' '{print $11}')
    INDEL_RECALL=$(echo "$INDEL_LINE" | awk -F',' '{print $10}')
    INDEL_F1=$(echo "$INDEL_LINE" | awk -F',' '{print $13}')

    # Default empty fields to N/A
    for var in SNP_TP SNP_FP SNP_FN SNP_PREC SNP_RECALL SNP_F1 \
               INDEL_TP INDEL_FP INDEL_FN INDEL_PREC INDEL_RECALL INDEL_F1; do
      [ -z "${!var}" ] && printf -v "$var" '%s' 'N/A'
    done

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$CALLER" "$SNP_TP" "$SNP_FP" "$SNP_FN" "$SNP_PREC" "$SNP_RECALL" "$SNP_F1" \
      "$INDEL_TP" "$INDEL_FP" "$INDEL_FN" "$INDEL_PREC" "$INDEL_RECALL" "$INDEL_F1" >> "$TSV"
  done

  # Print results table
  {
    echo "--- SNP Accuracy ---"
    printf "%-14s %10s %10s %10s %10s %10s %10s\n" "Caller" "TP" "FP" "FN" "Precision" "Recall" "F1"
    printf "%-14s %10s %10s %10s %10s %10s %10s\n" "--------------" "----------" "----------" "----------" "----------" "----------" "----------"
    while IFS=$'\t' read -r caller stp sfp sfn sprec srec sf1 _itp _ifp _ifn _iprec _irec _if1; do
      printf "%-14s %10s %10s %10s %10s %10s %10s\n" "$caller" "$stp" "$sfp" "$sfn" "$sprec" "$srec" "$sf1"
    done < <(tail -n +2 "$TSV")
    echo ""
    echo "--- INDEL Accuracy ---"
    printf "%-14s %10s %10s %10s %10s %10s %10s\n" "Caller" "TP" "FP" "FN" "Precision" "Recall" "F1"
    printf "%-14s %10s %10s %10s %10s %10s %10s\n" "--------------" "----------" "----------" "----------" "----------" "----------" "----------"
    while IFS=$'\t' read -r caller _stp _sfp _sfn _sprec _srec _sf1 itp ifp ifn iprec irec if1; do
      printf "%-14s %10s %10s %10s %10s %10s %10s\n" "$caller" "$itp" "$ifp" "$ifn" "$iprec" "$irec" "$if1"
    done < <(tail -n +2 "$TSV")
  } | tee -a "$SUMMARY"

# =============================================================================
# PAIRWISE CONCORDANCE MODE (bcftools isec)
# =============================================================================
else
  echo ""
  echo "=== Pairwise Concordance (bcftools isec) ==="

  # Build the TSV header
  printf "CallerA\tCallerB\tShared\tA_Unique\tB_Unique\tJaccard\n" > "$TSV"

  {
    echo "================================================================================"
    echo "  VARIANT CALLER PAIRWISE CONCORDANCE"
    echo "  Sample: ${SAMPLE}"
    echo "  Generated: $(date -u '+%Y-%m-%d %H:%M UTC')"
    echo "================================================================================"
    echo ""
  } > "$SUMMARY"

  [ -n "${INTERVALS:-}" ] && echo "Restricted to regions: ${INTERVALS}"

  ISEC_REGIONS_FLAG=""
  if [ -n "${INTERVALS:-}" ]; then
    ISEC_REGIONS_FLAG="-r ${INTERVALS}"
  fi

  for i in $(seq 0 $((NUM_CALLERS - 1))); do
    for j in $(seq $((i + 1)) $((NUM_CALLERS - 1))); do
      CALLER_A="${CALLER_NAMES[$i]}"
      CALLER_B="${CALLER_NAMES[$j]}"
      VCF_A="${CALLER_VCFS[$i]/#$GENOME_DIR//genome}"
      VCF_B="${CALLER_VCFS[$j]/#$GENOME_DIR//genome}"
      ISEC_DIR="/genome/${SAMPLE}/benchmark/isec_${CALLER_A}_vs_${CALLER_B}"
      ISEC_HOST="${BENCHMARK_DIR}/isec_${CALLER_A}_vs_${CALLER_B}"

      echo ""
      echo "=== Comparing ${CALLER_A} vs ${CALLER_B} ==="

      # Clean previous run if present
      rm -rf "$ISEC_HOST"

      # Normalize both VCFs (decompose MNPs, left-align indels) for fair comparison
      NORM_A="/genome/${SAMPLE}/benchmark/.norm_a.vcf.gz"
      NORM_B="/genome/${SAMPLE}/benchmark/.norm_b.vcf.gz"
      docker run --rm \
        --cpus 2 --memory 4g \
        --user root \
        -v "${GENOME_DIR}:/genome" \
        staphb/bcftools:1.21 \
        bash -c "
          bcftools norm -m-both -f /genome/reference/Homo_sapiens_assembly38.fasta ${VCF_A} -Oz -o ${NORM_A} && bcftools index -t ${NORM_A} &&
          bcftools norm -m-both -f /genome/reference/Homo_sapiens_assembly38.fasta ${VCF_B} -Oz -o ${NORM_B} && bcftools index -t ${NORM_B}
        "

      # shellcheck disable=SC2086
      docker run --rm \
        --cpus 2 --memory 4g \
        --user root \
        -v "${GENOME_DIR}:/genome" \
        staphb/bcftools:1.21 \
        bcftools isec -p "${ISEC_DIR}" \
          -f PASS \
          ${ISEC_REGIONS_FLAG} \
          "${NORM_A}" "${NORM_B}"

      # Clean up normalized temp files
      rm -f "${BENCHMARK_DIR}/.norm_a.vcf.gz"* "${BENCHMARK_DIR}/.norm_b.vcf.gz"*

      # Count variants in each output file
      # 0000.vcf = unique to A
      # 0001.vcf = unique to B
      # 0002.vcf = shared (from A's perspective)
      # 0003.vcf = shared (from B's perspective)
      A_UNIQUE=$(grep -cv '^#' "${ISEC_HOST}/0000.vcf" 2>/dev/null || echo "0")
      B_UNIQUE=$(grep -cv '^#' "${ISEC_HOST}/0001.vcf" 2>/dev/null || echo "0")
      SHARED=$(grep -cv '^#' "${ISEC_HOST}/0002.vcf" 2>/dev/null || echo "0")

      # Jaccard = shared / (A_unique + B_unique + shared)
      DENOM=$((A_UNIQUE + B_UNIQUE + SHARED))
      if [ "$DENOM" -gt 0 ]; then
        JACCARD=$(awk "BEGIN {printf \"%.4f\", ${SHARED} / ${DENOM}}")
      else
        JACCARD="0.0000"
      fi

      echo "  ${CALLER_A} unique: ${A_UNIQUE}"
      echo "  ${CALLER_B} unique: ${B_UNIQUE}"
      echo "  Shared:             ${SHARED}"
      echo "  Jaccard index:      ${JACCARD}"

      printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$CALLER_A" "$CALLER_B" "$SHARED" "$A_UNIQUE" "$B_UNIQUE" "$JACCARD" >> "$TSV"
    done
  done

  # Print results table
  {
    echo ""
    printf "%-14s %-14s %12s %12s %12s %10s\n" "Caller A" "Caller B" "Shared" "A Unique" "B Unique" "Jaccard"
    printf "%-14s %-14s %12s %12s %12s %10s\n" "--------------" "--------------" "------------" "------------" "------------" "----------"
    while IFS=$'\t' read -r ca cb shared au bu jaccard; do
      printf "%-14s %-14s %12s %12s %12s %10s\n" "$ca" "$cb" "$shared" "$au" "$bu" "$jaccard"
    done < <(tail -n +2 "$TSV")
    echo ""
    echo "NOTE: Jaccard index = shared / (A_unique + B_unique + shared)."
    echo "Computed on PASS variants after normalization (bcftools norm -m-both)."
    echo "A value of 1.0 means perfect agreement; typical WGS callers agree on ~95%+ of PASS SNPs."
  } | tee -a "$SUMMARY"
fi

echo ""
echo "=== Benchmarking complete ==="
echo "Summary: ${SUMMARY}"
echo "TSV:     ${TSV}"
