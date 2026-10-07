#!/usr/bin/env bash
# benchmark-variants.sh — Compare variant calls across multiple callers
# Supports pairwise concordance (bcftools isec) and truth set benchmarking (hap.py)
# Input: VCFs from vcf/, vcf_gatk/, vcf_freebayes/ directories
# Output: comparison tables in $GENOME_DIR/<sample>/benchmark/
#
# Truth sets: --truth VCF (with --regions BED) for any truth set, or
# --giab v4.2.1|v5.0q for one of GIAB's two HG002 GRCh38 small-variant
# benchmarks, downloaded into GENOME_DIR/giab/ on first use and checked by md5:
#   v4.2.1  NISTv4.2.1, chr1-22, mapping-based (the long-standing benchmark)
#   v5.0q   the draft assembly-based set from the T2T HG002 Q100 assembly,
#           with chrX, chrY and harder regions; hap.py gets --gender male
# With --giab, --regions keeps only the set's benchmark regions inside that BED
# (one chromosome, or a slice). Only an HG002 sequencing run gives meaningful
# numbers against either set.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name> [--truth <vcf> --regions <bed> | --giab v4.2.1|v5.0q [--regions <bed>]]}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
shift

# --- Parse optional flags ---
TRUTH_VCF=""
REGIONS_BED=""
GIAB=""
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
    --giab)
      GIAB="${2:?--giab requires v4.2.1 or v5.0q}"
      shift 2
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      echo "Usage: $0 <sample_name> [--truth <vcf> --regions <bed> | --giab v4.2.1|v5.0q [--regions <bed>]]" >&2
      exit 1
      ;;
  esac
done

# --- GIAB HG002 truth sets (--giab) ---
# Each set's VCF, its index and its benchmark-regions BED, with their md5:
# v5.0q's from GIAB's checksum.md5, v4.2.1's (no checksum file) as downloaded
# on 2026-10-07. An explicit --truth or --regions wins over the set's file.
HAPPY_EXTRA=()
if [ -n "$GIAB" ]; then
  [ -z "$TRUTH_VCF" ] || { echo "ERROR: use --giab or --truth, not both." >&2; exit 1; }
  GIAB_BASE=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son
  case "$GIAB" in
    v4.2.1)
      GIAB_URL="${GIAB_BASE}/NISTv4.2.1/GRCh38"
      GIAB_FILES=("HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz dc750b3807d4af1f7ffec852e9c2f771"
                  "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi 121e2975fb3ff0317ae6a684d0ce6f2f"
                  "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed 97265e922a97c69a0391cf3f92a89b8b") ;;
    v5.0q)
      GIAB_URL="${GIAB_BASE}/v5.0q"
      GIAB_FILES=("HG002_GRCh38_v5.0q_smvar.vcf.gz c71acc71069bf7cd7f51eb8fb0c1a1ab"
                  "HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi 58cab2a06e29b74a5bdc190ead079083"
                  "HG002_GRCh38_v5.0q_smvar.benchmark.bed 3366858af85875cbb3c579412ffa4c56") ;;
    *) echo "ERROR: --giab must be v4.2.1 or v5.0q, got '${GIAB}'" >&2; exit 1 ;;
  esac
  for entry in "${GIAB_FILES[@]}"; do
    f=${entry% *}
    if [ ! -s "${GENOME_DIR}/giab/${f}" ]; then
      echo "Downloading GIAB HG002 ${GIAB}: ${f}"
      fetch "${GIAB_URL}/${f}" "${GENOME_DIR}/giab/${f}" md5 "${entry#* }"
    fi
  done
  TRUTH_VCF="${GENOME_DIR}/giab/${GIAB_FILES[0]% *}"
  SET_BED="${GENOME_DIR}/giab/${GIAB_FILES[2]% *}"
  if [ -n "$REGIONS_BED" ]; then
    # --regions with --giab: the set's benchmark regions inside that BED only.
    [ -f "$REGIONS_BED" ] || { echo "ERROR: Regions BED not found: ${REGIONS_BED}" >&2; exit 1; }
    CUT="${GENOME_DIR}/${SAMPLE}/benchmark/giab_${GIAB}_regions.bed"
    mkdir -p "$(dirname "$CUT")"
    awk 'BEGIN {OFS = "\t"} FNR == NR { if ($0 !~ /^(#|track|browser)/ && NF >= 3) { n[$1]++; s[$1, n[$1]] = $2 + 0; e[$1, n[$1]] = $3 + 0 }; next }
      { for (i = 1; i <= n[$1]; i++) { a = ($2 + 0 > s[$1, i] ? $2 + 0 : s[$1, i]); b = ($3 + 0 < e[$1, i] ? $3 + 0 : e[$1, i]); if (a < b) print $1, a, b } }' \
      "$REGIONS_BED" "$SET_BED" | sort -k1,1 -k2,2n > "$CUT"
    if [ ! -s "$CUT" ]; then
      echo "ERROR: no ${GIAB} benchmark region lies inside ${REGIONS_BED}." >&2
      exit 1
    fi
    echo "Benchmark regions of ${GIAB} inside ${REGIONS_BED}: $(wc -l < "$CUT" | tr -d ' ') intervals"
    REGIONS_BED=$CUT
  else
    REGIONS_BED=$SET_BED
  fi
  # HG002 is male; hap.py does not infer the sex, and v5.0q has chrX and chrY.
  HAPPY_EXTRA=(--gender male)
  echo "Truth set: GIAB HG002 ${GIAB} (${TRUTH_VCF}), regions ${REGIONS_BED}"
fi

BENCHMARK_DIR="${GENOME_DIR}/${SAMPLE}/benchmark"
SUMMARY="${BENCHMARK_DIR}/summary.txt"
TSV="${BENCHMARK_DIR}/comparison.tsv"
REF="$REF_FASTA"

# --- Path validation helper ---
# Canonicalize GENOME_DIR and verify paths are within it (prevents sibling prefix matches)
# Uses python3 (macOS) with readlink -f fallback (GNU/Linux) — realpath -e is not portable.
_resolve() {
  python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$1" 2>/dev/null \
    || readlink -f "$1" 2>/dev/null \
    || echo "$1"
}
GENOME_DIR_REAL=$(_resolve "$GENOME_DIR")
[ -d "$GENOME_DIR_REAL" ] || { echo "ERROR: GENOME_DIR does not exist: ${GENOME_DIR}" >&2; exit 1; }

within_genome_dir() {
  local p
  p=$(_resolve "$1")
  [ -e "$p" ] || return 1
  [[ "$p" == "$GENOME_DIR_REAL" || "$p" == "$GENOME_DIR_REAL"/* ]]
}

# --- Discover available caller VCFs ---
declare -a CALLER_NAMES=()
declare -a CALLER_VCFS=()

DV_VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
GATK_VCF="${GENOME_DIR}/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz"
FB_VCF="${GENOME_DIR}/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz"
SK_VCF="${GENOME_DIR}/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz"
OCT_VCF="${GENOME_DIR}/${SAMPLE}/vcf_octopus/${SAMPLE}.vcf.gz"

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
if [ -f "$OCT_VCF" ]; then
  CALLER_NAMES+=("Octopus")
  CALLER_VCFS+=("$OCT_VCF")
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
  echo "  vcf_gatk/      (GATK)       — run scripts/03a-gatk-haplotypecaller.sh" >&2
  echo "  vcf_freebayes/ (FreeBayes)  — run scripts/03b-freebayes.sh" >&2
  echo "  vcf_strelka2/  (Strelka2)   — run scripts/03c-strelka2-germline.sh" >&2
  echo "  vcf_octopus/   (Octopus)    — run scripts/03d-octopus.sh" >&2
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

# Validate reference FASTA (needed by both pairwise normalization and truth mode)
for f in "$REF" "${REF}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    exit 1
  fi
done

# Validate regions BED (used by both modes)
if [ -n "$REGIONS_BED" ]; then
  if [ ! -f "$REGIONS_BED" ]; then
    echo "ERROR: Regions BED not found: ${REGIONS_BED}" >&2
    exit 1
  fi
  if ! within_genome_dir "$REGIONS_BED"; then
    echo "ERROR: --regions path must be under GENOME_DIR (${GENOME_DIR})." >&2
    echo "The container only mounts GENOME_DIR as /genome." >&2
    exit 1
  fi
fi

# Validate truth mode inputs
if [ -n "$TRUTH_VCF" ]; then
  for f in "$TRUTH_VCF" "${TRUTH_VCF}.tbi"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: File not found: ${f}" >&2
      exit 1
    fi
  done
  if ! within_genome_dir "$TRUTH_VCF"; then
    echo "ERROR: --truth path must be under GENOME_DIR (${GENOME_DIR})." >&2
    echo "The container only mounts GENOME_DIR as /genome." >&2
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
  if [ -z "$REGIONS_BED" ]; then
    echo "WARNING: No --regions BED provided. For GIAB truth sets, you should provide"
    echo "the high-confidence regions BED (e.g., HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed)."
    echo "Without it, variants outside confident regions inflate FP/FN counts and make"
    echo "precision/recall harder to interpret."
    echo ""
  fi

  # Build the TSV header
  printf "Caller\tSNP_TP\tSNP_FP\tSNP_FN\tSNP_Precision\tSNP_Recall\tSNP_F1\tINDEL_TP\tINDEL_FP\tINDEL_FN\tINDEL_Precision\tINDEL_Recall\tINDEL_F1\n" > "$TSV"

  {
    echo "================================================================================"
    echo "  VARIANT CALLER TRUTH SET BENCHMARK"
    echo "  Sample: ${SAMPLE}"
    echo "  Truth:  $(basename "$TRUTH_VCF")${GIAB:+ (GIAB HG002 ${GIAB})}"
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
    VCF_CONTAINER="${VCF/#$GENOME_DIR//genome}"
    PREFIX="/genome/${SAMPLE}/benchmark/happy_${CALLER}"

    echo ""
    echo "=== Running hap.py: ${CALLER} vs truth ==="

    # --root: the hap.py image has not been shown to run as an unprivileged user.
    # shellcheck disable=SC2086
    run_in --root \
      --cpus 4 --memory 8g \
      "${HAPPY_IMAGE}" \
      /opt/hap.py/bin/hap.py \
        "${TRUTH_CONTAINER_PATH}" \
        "${VCF_CONTAINER}" \
        -r "${REF_FASTA_C}" \
        ${REGIONS_FLAG} \
        ${HAPPY_EXTRA[@]+"${HAPPY_EXTRA[@]}"} \
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
    INDEL_LINE=$(awk -F',' '$1=="INDEL" && $2=="PASS"' "$HAPPY_CSV" || true)

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

  ISEC_REGIONS_FLAG=""
  if [ -n "$REGIONS_BED" ]; then
    REGIONS_BED_CONTAINER="${REGIONS_BED/#$GENOME_DIR//genome}"
    ISEC_REGIONS_FLAG="-R ${REGIONS_BED_CONTAINER}"
    echo "Restricted to regions BED: ${REGIONS_BED}"
  elif [ -n "${INTERVALS:-}" ]; then
    ISEC_REGIONS_FLAG="-r ${INTERVALS}"
    echo "Restricted to region: ${INTERVALS}"
  fi

  # Normalise each caller's VCF once (split multiallelics, left-align indels)
  # for a fair comparison; every pair below reads the normalised copies.
  declare -a NORM_VCFS=()
  for i in $(seq 0 $((NUM_CALLERS - 1))); do
    NORM="/genome/${SAMPLE}/benchmark/.norm_${CALLER_NAMES[$i]}.vcf.gz"
    echo "Normalising ${CALLER_NAMES[$i]}..."
    run_in \
      --cpus 2 --memory 4g \
      "${BCFTOOLS_IMAGE}" \
      bash -euo pipefail -c \
        'bcftools norm -m-both -f "$3" "$1" -Oz -o "$2" && bcftools index -f -t "$2"' \
        _ "${CALLER_VCFS[$i]/#$GENOME_DIR//genome}" "$NORM" "${REF_FASTA_C}"
    NORM_VCFS+=("$NORM")
  done

  for i in $(seq 0 $((NUM_CALLERS - 1))); do
    # Not `seq i+1 N-1`: BSD seq (macOS) counts down when i+1 > N-1.
    for ((j = i + 1; j < NUM_CALLERS; j++)); do
      CALLER_A="${CALLER_NAMES[$i]}"
      CALLER_B="${CALLER_NAMES[$j]}"
      ISEC_DIR="/genome/${SAMPLE}/benchmark/isec_${CALLER_A}_vs_${CALLER_B}"
      ISEC_HOST="${BENCHMARK_DIR}/isec_${CALLER_A}_vs_${CALLER_B}"

      echo ""
      echo "=== Comparing ${CALLER_A} vs ${CALLER_B} ==="

      # Clean previous run if present
      rm -rf "$ISEC_HOST"

      NORM_A="${NORM_VCFS[$i]}"
      NORM_B="${NORM_VCFS[$j]}"

      # shellcheck disable=SC2086
      run_in \
        --cpus 2 --memory 4g \
        "${BCFTOOLS_IMAGE}" \
        bcftools isec -p "${ISEC_DIR}" \
          -f .,PASS \
          ${ISEC_REGIONS_FLAG} \
          "${NORM_A}" "${NORM_B}"

      # Count variants in each output file
      # 0000.vcf = unique to A
      # 0001.vcf = unique to B
      # 0002.vcf = shared (from A's perspective)
      # 0003.vcf = shared (from B's perspective)
      A_UNIQUE=$(grep -c '^[^#]' "${ISEC_HOST}/0000.vcf" 2>/dev/null) || A_UNIQUE=0
      B_UNIQUE=$(grep -c '^[^#]' "${ISEC_HOST}/0001.vcf" 2>/dev/null) || B_UNIQUE=0
      SHARED=$(grep -c '^[^#]' "${ISEC_HOST}/0002.vcf" 2>/dev/null) || SHARED=0

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

  # Clean up the normalised copies
  rm -f "${BENCHMARK_DIR}"/.norm_*.vcf.gz "${BENCHMARK_DIR}"/.norm_*.vcf.gz.tbi

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
    echo "Computed on PASS and unfiltered (.) variants after normalization (bcftools norm -m-both)."
    echo "A value of 1.0 means perfect agreement; typical WGS callers agree on ~95%+ of SNPs."
  } | tee -a "$SUMMARY"
fi

echo ""
echo "=== Benchmarking complete ==="
echo "Summary: ${SUMMARY}"
echo "TSV:     ${TSV}"
