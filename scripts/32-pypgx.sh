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

# CYP2D6 depth check: mosdepth over CYP2D6 and its flanks, all reads and
# MAPQ >= 1, judged by bin/cyp2d6_depth_check.py (the regions come from it too).
DEPTH_DIR="${OUTPUT_DIR}/cyp2d6_depth"
CHECK="${OUTPUT_DIR}/${SAMPLE}_cyp2d6_depth_check.tsv"
mkdir -p "$DEPTH_DIR"
rm -f "$CHECK"
echo "CYP2D6 depth check..."
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py bed > "${DEPTH_DIR}/regions.bed"
for Q in 0 1; do
  run_in --cpus 2 --memory 2g "${MOSDEPTH_IMAGE}" \
    mosdepth -n -c chr22 -t 2 -Q "$Q" -b "$(cpath "${DEPTH_DIR}/regions.bed")" \
      "$(cpath "${DEPTH_DIR}/q${Q}")" "$(cpath "$BAM")"
done
run_in -v "${PGP_ROOT}/bin:/pgp-bin:ro" "${PYTHON_IMAGE}" \
  python3 /pgp-bin/cyp2d6_depth_check.py check \
    --all "$(cpath "${DEPTH_DIR}/q0.regions.bed.gz")" \
    --mapq1 "$(cpath "${DEPTH_DIR}/q1.regions.bed.gz")" \
    --out "$(cpath "$CHECK")"
# A check that wrote nothing counts as failed: CYP2D6 is then Indeterminate.
DEPTH_STATUS=$(awk -F'\t' '$1 == "status" {print $2}' "$CHECK" 2>/dev/null || true)

# Curated gene list: CPIC Level A/B + key genes PharmCAT misses
# BAM-based (structural variation): CYP2D6, CYP2A6, GSTM1, GSTT1
# VCF-based (additional coverage): CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP3A4, CYP3A5,
#   CYP4F2, DPYD, TPMT, NUDT15, UGT1A1, SLCO1B1, VKORC1, NAT2, COMT, MTHFR, ABCB1, G6PD, IFNL3
BAM_GENES="CYP2D6 CYP2A6 GSTM1 GSTT1"
VCF_GENES="CYP1A2 CYP2B6 CYP2C9 CYP2C19 CYP3A4 CYP3A5 CYP4F2 DPYD TPMT NUDT15 UGT1A1 SLCO1B1 VKORC1 NAT2 COMT MTHFR ABCB1 G6PD IFNL3"

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
    BAM_GENES="'"${BAM_GENES}"'"
    VCF_GENES="'"${VCF_GENES}"'"
    OUTBASE="/genome/${SAMPLE}/pypgx"
    BAM="/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"
    VCF="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
    DOC="${OUTBASE}/depth_of_coverage.zip"
    CTRL="${OUTBASE}/control_statistics.zip"
    FAILED=""
    SUCCEEDED=0

    # GSTT1 lies on chr22_KI270879v1_alt in GRCh38. A BAM aligned to a reference
    # without ALT contigs has no such contig, and depth preparation then fails for
    # every SV gene. Leave GSTT1 out in that case and say so.
    if ! python3 -c "import pysam, sys; sys.exit(0 if \"chr22_KI270879v1_alt\" in pysam.AlignmentFile(sys.argv[1]).references else 1)" "$BAM"; then
      echo "NOTICE: the BAM has no chr22_KI270879v1_alt contig (reference without ALT contigs); GSTT1 cannot be called from depth and is skipped"
      BAM_GENES=$(echo "$BAM_GENES" | tr " " "\n" | grep -vx GSTT1 | tr "\n" " ")
      FAILED="${FAILED} GSTT1"
    fi

    # Phase 1: Prepare depth of coverage for the SV genes (one-time, from BAM)
    echo "--- Preparing depth of coverage for SV genes: ${BAM_GENES} ---"
    if ! pypgx prepare-depth-of-coverage \
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

# Consolidate per-gene results into a summary TSV. It is written as .partial
# and renamed only after the depth check's verdict is applied below, so a run
# that stops in between never leaves an unchecked CYP2D6 call for step 36.
run_in --cpus 2 --memory 4g \
  -v "${PYPGX_BUNDLE}:/tmp/pypgx-bundle:ro" -e PYPGX_BUNDLE=/tmp/pypgx-bundle \
  "${PYPGX_IMAGE}" \
  python3 -c "
import os, sys, csv, subprocess

sample = '${SAMPLE}'
outbase = f'/genome/{sample}/pypgx'
bam_genes = '${BAM_GENES}'.split()
vcf_genes = '${VCF_GENES}'.split()
all_genes = bam_genes + vcf_genes

summary_path = f'{outbase}/{sample}_pypgx_summary.tsv'
rows = []

for gene in all_genes:
    results_zip = f'{outbase}/{gene}/results.zip'
    if not os.path.isfile(results_zip):
        rows.append([gene, 'FAILED', 'N/A', 'N/A', 'BAM' if gene in bam_genes else 'VCF'])
        continue

    diplotype = 'N/A'
    phenotype = 'N/A'
    cnv = 'N/A'
    source = 'BAM' if gene in bam_genes else 'VCF'
    # pypgx print-data results.zip outputs a TSV whose first column is the sample
    # and whose named columns include Genotype, Phenotype and CNV (the copy-number
    # call pypgx made from read depth, e.g. Normal or WholeDel1).
    out = subprocess.run(['pypgx', 'print-data', results_zip], capture_output=True, text=True)
    if out.returncode != 0:
        print(f'WARNING: pypgx print-data failed for {gene}: {out.stderr.strip()}', file=sys.stderr)
    else:
        lines = out.stdout.rstrip().split('\n')
        if len(lines) >= 2:
            headers = lines[0].split('\t')
            values = lines[1].split('\t')
            def col(name):
                if name in headers and headers.index(name) < len(values):
                    return values[headers.index(name)]
                return None
            diplotype = col('Genotype') or 'N/A'
            phenotype = col('Phenotype') or 'N/A'
            # Copy number is called from read depth, so it is meaningful only for the
            # BAM-based genes. Report pypgx's own value, never a guess from allele names.
            if source == 'BAM':
                cnv = col('CNV')
                if cnv is None:
                    print(f'WARNING: {gene} results have no CNV column', file=sys.stderr)
                    cnv = 'N/A'
    rows.append([gene, diplotype, phenotype, cnv, source])

with open(summary_path + '.partial', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['Gene', 'Diplotype', 'Phenotype', 'CNV_call', 'Source'])
    w.writerows(rows)

print(f'Genes called: {sum(1 for r in rows if r[1] != \"FAILED\")}/{len(rows)}')
" 2>&1

# A CYP2D6 call from multi-mapped depth is not a call: the row says so; the
# call itself stays in CYP2D6/results.zip.
SUMMARY="${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
if [ "$DEPTH_STATUS" != ok ] && [ -f "${SUMMARY}.partial" ]; then
  awk -F'\t' -v OFS='\t' '$1 == "CYP2D6" {$2 = "Indeterminate"; $3 = "Indeterminate (CYP2D6 depth check)"} {print}' \
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
