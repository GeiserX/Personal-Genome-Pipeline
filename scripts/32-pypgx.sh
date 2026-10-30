#!/usr/bin/env bash
# pypgx — Comprehensive pharmacogenomic star allele calling with SV detection
# Input: BAM + VCF from alignment/variant calling steps
# Output: Per-gene star allele calls and a consolidated summary TSV (step 27
#         compares it with PharmCAT; step 36 compares its CYP2D6 with Cyrius)
#
# CYP2D6 copy number comes from read depth. Before pypgx runs, mosdepth
# measures the depth over CYP2D6 and its flanks (bin/cyp2d6_depth_check.py).
# When the reads there are multi-mapped (a BAM aligned to a reference with
# ALT contigs), the summary's CYP2D6 row says Indeterminate instead of the
# call, and <sample>_cyp2d6_depth_check.tsv says why.
#
# The summary (bin/pgx_outside_calls.py pypgx-summary) also looks up, for each
# star allele called, the variants that define it in the VCF, and marks a call
# 'low-confidence: check the reads' when one of them has a GQ below
# PYPGX_GQ_FLOOR (default 20). GSTT1, which GRCh38 has only on an ALT contig,
# is 'not callable on this reference' on a BAM without that contig. The gene
# list and the mode of each gene are assets/pgx_genes.tsv.
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"

BAM="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTPUT_DIR="${GENOME_DIR}/${SAMPLE}/pypgx"

echo "=== pypgx Pharmacogenomics: ${SAMPLE} ==="
echo "Input BAM: ${BAM}"
echo "Input VCF: ${VCF}"
echo "Output:    ${OUTPUT_DIR}/"

# Validate inputs
for f in "$BAM" "${BAM}.bai" "$VCF" "${VCF}.tbi"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: File not found: ${f}" >&2
    if [ "$f" = "${BAM}.bai" ]; then
      echo "  Generate BAM index with: samtools index ${BAM}" >&2
    elif [ "$f" = "${VCF}.tbi" ]; then
      echo "  Generate tabix index with: bcftools index -t ${VCF}" >&2
    fi
    exit 1
  fi
done

mkdir -p "$OUTPUT_DIR"
# A run that stops anywhere below must not leave the last run's summary and
# depth check behind for step 36 to read.
rm -f "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv" "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv.partial" \
  "${OUTPUT_DIR}/${SAMPLE}_cyp2d6_depth_check.tsv"

# This step feeds step 36's outside calls for PharmCAT. Remove the ones made
# from an earlier result, so step 07 never reads a call this run has not
# confirmed; step 36 writes them again.
rm -f "${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_outside_calls.tsv" \
  "${GENOME_DIR}/${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv"

# Validate pypgx-bundle (required for Beagle phasing panels and CNV models)
PYPGX_BUNDLE="${GENOME_DIR}/reference/pypgx-bundle"
if [ ! -d "$PYPGX_BUNDLE" ]; then
  echo "ERROR: pypgx-bundle not found at ${PYPGX_BUNDLE}" >&2
  echo "  Download it (370 MB, one-time) with:" >&2
  echo "  cd ${GENOME_DIR}/reference && git clone --branch ${PYPGX_BUNDLE_VERSION} --depth 1 https://github.com/sbslee/pypgx-bundle.git" >&2
  exit 1
fi
# The bundle must be the tag that matches the pypgx image: with another tag
# every gene fails. PYPGX_BUNDLE_VERSION in versions.env names it.
if ! command -v git >/dev/null 2>&1; then
  echo "ERROR: git is needed to check the pypgx-bundle tag at ${PYPGX_BUNDLE}; install git and run again." >&2
  exit 1
fi
BUNDLE_TAG=$(git -c safe.directory="$PYPGX_BUNDLE" -C "$PYPGX_BUNDLE" describe --tags 2>/dev/null || true)
if [ "$BUNDLE_TAG" != "$PYPGX_BUNDLE_VERSION" ]; then
  echo "ERROR: pypgx-bundle at ${PYPGX_BUNDLE} is '${BUNDLE_TAG:-not a git checkout of a tag}', but ${PYPGX_IMAGE} needs ${PYPGX_BUNDLE_VERSION}." >&2
  echo "  Replace it with:" >&2
  echo "  git clone --branch ${PYPGX_BUNDLE_VERSION} --depth 1 https://github.com/sbslee/pypgx-bundle.git ${PYPGX_BUNDLE}" >&2
  exit 1
fi

# The genes and how each is called come from assets/pgx_genes.tsv (one row per
# gene; PGX_GENES names another table). bam rows get copy number from read
# depth plus the VCF, vcf rows the VCF only. tests/test_pypgx_genes.sh checks
# the bam rows against pypgx's own gene table (SV=TRUE) at the pinned tag.
PGX_GENES=${PGX_GENES:-${PGP_ROOT}/assets/pgx_genes.tsv}
if [ ! -s "$PGX_GENES" ]; then
  echo "ERROR: gene table not found or empty: ${PGX_GENES}" >&2
  exit 1
fi
# A malformed row stops the step: a gene name goes into a command line, and a
# typo in a mode would drop a gene without a word.
BAD=$(awk -F'\t' '/^#/ || $1 == "gene" || NF == 0 {next}
  NF != 6 || $1 !~ /^[A-Za-z0-9-]+$/ || ($2 != "bam" && $2 != "vcf") || ($3 != "yes" && $3 != "no") \
     || ($3 == "no" && $4 !~ /^[A-Za-z0-9_.-]+$/) || ($3 == "no" && $4 == ".") || seen[$1]++ {print NR": "$0}' "$PGX_GENES")
if [ -n "$BAD" ]; then
  echo "ERROR: ${PGX_GENES} has rows that are not gene<TAB>bam|vcf<TAB>yes|no<TAB>alt_contig<TAB>second_caller<TAB>pharmcat_outside_call, or a gene twice:" >&2
  echo "$BAD" >&2
  exit 1
fi
genes_of() { awk -F'\t' -v m="$1" '!/^#/ && $1 != "gene" && $2 == m {printf "%s%s", s, $1; s = " "}' "$PGX_GENES"; }
BAM_GENES=$(genes_of bam)
VCF_GENES=$(genes_of vcf)
if [ -z "${BAM_GENES}${VCF_GENES}" ]; then
  echo "ERROR: ${PGX_GENES} lists no gene" >&2
  exit 1
fi

GQ_FLOOR=${PYPGX_GQ_FLOOR:-20}
[[ "$GQ_FLOOR" =~ ^[0-9]+$ ]] || { echo "ERROR: PYPGX_GQ_FLOOR must be a whole number, got '${GQ_FLOOR}'" >&2; exit 1; }

# CYP2D6 depth check: mosdepth over CYP2D6 and its flanks, all reads and
# MAPQ >= 1, judged by bin/cyp2d6_depth_check.py (the regions come from it too).
# Beside the verdict, `measure` counts the MAPQ >= 20 depth and the gene body
# without REP6 from the BAM (pysam, in the pypgx image, as the PYPGX module does).
DEPTH_DIR="${OUTPUT_DIR}/cyp2d6_depth"
CHECK="${OUTPUT_DIR}/${SAMPLE}_cyp2d6_depth_check.tsv"
SHARES="${DEPTH_DIR}/mapq_shares.tsv"
mkdir -p "$DEPTH_DIR"
rm -f "$CHECK" "$SHARES"
echo "CYP2D6 depth check..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py bed > "${DEPTH_DIR}/regions.bed"
for Q in 0 1; do
  run_in --cpus 2 --memory 2g "${MOSDEPTH_IMAGE}" \
    mosdepth -n -c chr22 -t 2 -Q "$Q" -b "$(cpath "${DEPTH_DIR}/regions.bed")" \
      "$(cpath "${DEPTH_DIR}/q${Q}")" "$(cpath "$BAM")"
done
MEASURE=()
if run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYPGX_IMAGE}" \
     python3 /pgp-bin/cyp2d6_depth_check.py measure --bam "$(cpath "$BAM")" --out "$(cpath "$SHARES")" \
   && [ -s "$SHARES" ]; then
  MEASURE=(--measure "$(cpath "$SHARES")")
else
  echo "WARNING: the MAPQ >= 20 and gene-body depth at CYP2D6 was not measured; the check says 'not measured'"
  rm -f "$SHARES"
fi
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py check \
    --all "$(cpath "${DEPTH_DIR}/q0.regions.bed.gz")" \
    --mapq1 "$(cpath "${DEPTH_DIR}/q1.regions.bed.gz")" \
    "${MEASURE[@]}" --out "$(cpath "$CHECK")"
# A check that wrote nothing counts as failed: CYP2D6 is then Indeterminate.
DEPTH_STATUS=$(awk -F'\t' '$1 == "status" {print $2}' "$CHECK" 2>/dev/null || true)

# A gene whose callable_without_alt is no lies on an ALT contig of GRCh38
# (GSTT1 on chr22_KI270879v1_alt). A BAM aligned to a reference without ALT
# contigs has no such contig, and depth preparation then fails for every bam
# gene. Leave such a gene out in that case: the summary says it is not
# callable on this reference (a property of the reference, not a tool error).
NOT_CALLABLE=() CALL_BAM_GENES=$BAM_GENES
mapfile -t NEED < <(awk -F'\t' '!/^#/ && $1 != "gene" && $2 == "bam" && $3 == "no" {print $1 "\t" $4}' "$PGX_GENES")
if [ "${#NEED[@]}" -gt 0 ]; then
  CONTIGS=()
  for n in "${NEED[@]}"; do CONTIGS+=("${n#*$'\t'}"); done
  rc=0
  MISSING=$(run_in "${PYPGX_IMAGE}" python3 -c \
    "import pysam, sys; refs = set(pysam.AlignmentFile(sys.argv[1]).references); print(' '.join(c for c in sys.argv[2:] if c not in refs))" \
    "$(cpath "$BAM")" "${CONTIGS[@]}") || rc=$?
  if [ "$rc" -ne 0 ]; then
    echo "ERROR: could not read the contigs of ${BAM} (exit ${rc})" >&2
    exit 1
  fi
  for n in "${NEED[@]}"; do
    g=${n%%$'\t'*} c=${n#*$'\t'}
    case " ${MISSING} " in
      *" ${c} "*)
        echo "NOTICE: the BAM has no ${c} contig (reference without ALT contigs); ${g} is not callable on this reference and is skipped"
        NOT_CALLABLE+=("$g")
        CALL_BAM_GENES=$(tr ' ' '\n' <<<"$CALL_BAM_GENES" | grep -vxF "$g" | paste -sd ' ' -) || CALL_BAM_GENES="" ;;
    esac
  done
fi

echo ""
echo "Running pypgx for $(echo "$BAM_GENES" "$VCF_GENES" | wc -w | tr -d ' ') genes..."
echo "BAM-based (SV detection): ${BAM_GENES}"
echo "VCF-based: ${VCF_GENES}"
echo ""

# Run all genes in a single Docker container to avoid repeated startup overhead.
# pypgx requires a two-phase setup before calling genes:
#   1. prepare-depth-of-coverage — compute read depth from BAM for SV genes
#   2. compute-control-statistics — normalize read depth using a control gene (VDR)
# Then per-gene calling:
#   - SV genes: --depth-of-coverage + --control-statistics (no --variants to avoid
#     pseudogene-confounded VCF calls in CYP2D6/CYP2D7 region)
#   - VCF genes: --variants only
# Individual gene failures are logged but do not stop the loop.
run_in --cpus 4 --memory 8g \
  -v "${PYPGX_BUNDLE}:/tmp/pypgx-bundle:ro" -e PYPGX_BUNDLE=/tmp/pypgx-bundle \
  "${PYPGX_IMAGE}" \
  bash -c '
    SAMPLE="'"${SAMPLE}"'"
    BAM_GENES="'"${CALL_BAM_GENES}"'"
    VCF_GENES="'"${VCF_GENES}"'"
    OUTBASE="/genome/${SAMPLE}/pypgx"
    BAM="/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
    VCF="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
    DOC="${OUTBASE}/depth_of_coverage.zip"
    CTRL="${OUTBASE}/control_statistics.zip"
    FAILED=""
    SUCCEEDED=0

    # Phase 1: Prepare depth of coverage for the SV genes (one-time, from BAM)
    echo "--- Preparing depth of coverage for SV genes: ${BAM_GENES:-none} ---"
    if [ -z "$BAM_GENES" ]; then
      DOC=""
    elif ! pypgx prepare-depth-of-coverage \
      "$DOC" "$BAM" --assembly GRCh38 --genes $BAM_GENES 2>&1; then
      echo "ERROR: prepare-depth-of-coverage failed — cannot call SV genes"
      # Fall through to VCF-only genes; mark all BAM genes as failed
      for GENE in $BAM_GENES; do FAILED="${FAILED} ${GENE}"; done
      DOC=""
    fi

    # Phase 2: Compute control statistics from VDR for read-depth normalization
    if [ -n "$DOC" ]; then
      echo "--- Computing control statistics (VDR) ---"
      if ! pypgx compute-control-statistics \
        VDR "$CTRL" "$BAM" --assembly GRCh38 2>&1; then
        echo "WARNING: compute-control-statistics failed; SV calling proceeds without normalization"
        CTRL=""
      fi
    fi

    # Phase 3a: BAM-based genes — SV detection via read depth + VCF variants
    # Uses both --variants and --depth-of-coverage per upstream WGS workflow
    if [ -n "$DOC" ]; then
      for GENE in $BAM_GENES; do
        echo "--- Calling ${GENE} (BAM + VCF) ---"
        EXTRA=""
        [ -f "$CTRL" ] && EXTRA="--control-statistics $CTRL"
        pypgx run-ngs-pipeline "$GENE" "${OUTBASE}/${GENE}" \
          --variants "$VCF" \
          --depth-of-coverage "$DOC" \
          --assembly GRCh38 \
          --force \
          $EXTRA 2>&1 \
          && SUCCEEDED=$((SUCCEEDED + 1)) \
          || { echo "WARNING: ${GENE} failed"; FAILED="${FAILED} ${GENE}"; }
      done
    fi

    # Phase 3b: VCF-based genes — star alleles from variant calls only
    for GENE in $VCF_GENES; do
      echo "--- Calling ${GENE} (VCF-based) ---"
      pypgx run-ngs-pipeline "$GENE" "${OUTBASE}/${GENE}" \
        --variants "$VCF" \
        --assembly GRCh38 \
        --force 2>&1 \
        && SUCCEEDED=$((SUCCEEDED + 1)) \
        || { echo "WARNING: ${GENE} failed"; FAILED="${FAILED} ${GENE}"; }
    done

    echo ""
    echo "pypgx pipeline: ${SUCCEEDED} genes succeeded"
    if [ -n "$FAILED" ]; then
      echo "Failed genes:${FAILED}"
    fi
    # Exit non-zero if ALL genes failed
    [ "$SUCCEEDED" -gt 0 ] || exit 1
  '

echo ""
echo "Extracting results and building summary..."

# Consolidate per-gene results into a summary TSV (bin/pgx_outside_calls.py
# pypgx-summary, the same code as the PYPGX module). It is written as .partial
# and renamed only after the depth check's verdict is applied below, so a run
# that stops in between never leaves an unchecked CYP2D6 call for step 36.
SUMMARY="${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
NC_ARG=()
for g in ${NOT_CALLABLE[@]+"${NOT_CALLABLE[@]}"}; do NC_ARG+=(--not-callable "$g"); done
run_in --cpus 2 --memory 4g \
  -v "${PYPGX_BUNDLE}:/tmp/pypgx-bundle:ro" -e PYPGX_BUNDLE=/tmp/pypgx-bundle \
  -v "${PGP_ROOT}/bin:/pgp-bin:ro" \
  "${PYPGX_IMAGE}" \
  python3 /pgp-bin/pgx_outside_calls.py pypgx-summary \
    --outbase "$(cpath "$OUTPUT_DIR")" --bam-genes "$BAM_GENES" --vcf-genes "$VCF_GENES" \
    --vcf "$(cpath "$VCF")" --gq-floor "$GQ_FLOOR" ${NC_ARG[@]+"${NC_ARG[@]}"} \
    --out "$(cpath "${SUMMARY}.partial")" 2>&1

# A CYP2D6 call from multi-mapped depth is not a call: the row says so; the
# call itself stays in CYP2D6/results.zip.
if [ "$DEPTH_STATUS" != ok ] && [ -f "${SUMMARY}.partial" ]; then
  awk -F'\t' -v OFS='\t' '$1 == "CYP2D6" {$2 = "Indeterminate"; $3 = "Indeterminate (CYP2D6 depth check)"; $4 = "N/A (CYP2D6 depth check)"; $6 = "no call (CYP2D6 depth check)"} {print}' \
    "${SUMMARY}.partial" > "${SUMMARY}.tmp"
  mv "${SUMMARY}.tmp" "${SUMMARY}.partial"
  MSG=$(awk -F'\t' '$1 == "message" {print $2}' "$CHECK" 2>/dev/null || true)
  echo "WARNING: ${MSG:-the CYP2D6 depth check wrote no result}"
  echo "  The CYP2D6 row of ${SUMMARY} says Indeterminate; see ${CHECK}."
fi
mv "${SUMMARY}.partial" "$SUMMARY"
echo "Summary written: ${SUMMARY}"

# The PharmCAT comparison is written by step 27 (CPIC lookup), which runs after
# both PharmCAT (step 7) and this step; here it could read a missing or
# previous-run PharmCAT report, because run-all.sh starts steps 7 and 32 together.
# Step 27 removes and rewrites the file each time it runs, so this step leaves it
# alone: re-running step 32 on its own keeps the comparison the reports read.

# Print summary
echo ""
echo "============================================"
echo "  pypgx complete: ${SAMPLE}"
echo "============================================"
echo "Results:    ${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
echo "Comparison with PharmCAT: run step 27, which writes ${OUTPUT_DIR}/${SAMPLE}_pharmcat_comparison.tsv"
echo ""
if [ -f "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv" ]; then
  echo "Summary:"
  column -t -s $'\t' "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv" 2>/dev/null || cat "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
fi
