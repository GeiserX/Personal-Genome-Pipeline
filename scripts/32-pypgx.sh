#!/usr/bin/env bash
# pypgx — Comprehensive pharmacogenomic star allele calling with SV detection
# Input: BAM + VCF from alignment/variant calling steps
# Output: Per-gene star allele calls, consolidated summary TSV, PharmCAT comparison TSV
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

# Consolidate per-gene results into a summary TSV
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

with open(summary_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['Gene', 'Diplotype', 'Phenotype', 'CNV_call', 'Source'])
    w.writerows(rows)

print(f'Summary written: {summary_path}')
print(f'Genes called: {sum(1 for r in rows if r[1] != \"FAILED\")}/{len(rows)}')
" 2>&1

# Cross-reference with PharmCAT if output exists (newest report wins)
PHARMCAT_JSON=""
for DIR in "${GENOME_DIR}/${SAMPLE}/pharmcat" "${GENOME_DIR}/${SAMPLE}/vcf"; do
  [ -d "$DIR" ] || continue
  CANDIDATE=$(find "$DIR" -maxdepth 1 \( -name "*.report.json" -o -name "*_pharmcat.json" \) -print0 2>/dev/null \
    | xargs -0 ls -t 2>/dev/null | head -1)
  if [ -n "$CANDIDATE" ]; then
    PHARMCAT_JSON="$CANDIDATE"
    break
  fi
done

if [ -n "$PHARMCAT_JSON" ]; then
  echo ""
  echo "PharmCAT output found, generating comparison..."

  run_in --cpus 2 --memory 4g \
    "${PYTHON_IMAGE}" \
    python3 -c "
import json, csv, os, re, sys

sample = '${SAMPLE}'
outbase = f'/genome/{sample}/pypgx'
comparison_path = f'{outbase}/{sample}_pharmcat_comparison.tsv'

# Load pypgx summary
pypgx_data = {}
summary_path = f'{outbase}/{sample}_pypgx_summary.tsv'
if os.path.isfile(summary_path):
    with open(summary_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pypgx_data[row['Gene']] = row['Diplotype']

# Load PharmCAT results. PharmCAT 3.x writes 'genes' either flat ({gene -> data})
# or nested ({source -> {gene -> data}}); 2.x used a list. Same logic as
# scripts/27-cpic-lookup.sh. A report that cannot be read, or that yields no gene,
# is an error: an empty comparison would show 0 conflicts.
pharmcat_path = '$(echo "$PHARMCAT_JSON" | sed "s|${GENOME_DIR}|/genome|")'
with open(pharmcat_path) as f:
    data = json.load(f)

def parse_gene(g):
    if not isinstance(g, dict):
        return None
    dips = g.get('sourceDiplotypes') or g.get('recommendationDiplotypes') or []
    if not dips:
        return None
    dip = dips[0]
    a1 = (dip.get('allele1') or {}).get('name', '?')
    a2 = (dip.get('allele2') or {}).get('name', '?')
    return dip.get('label') or f'{a1}/{a2}'

pharmcat_data = {}
genes = data.get('genes')
if isinstance(genes, dict):
    for key, val in genes.items():
        if not isinstance(val, dict):
            continue
        if 'sourceDiplotypes' in val or 'recommendationDiplotypes' in val:
            d = parse_gene(val)                       # flat: key is the gene
            if d and key not in pharmcat_data:
                pharmcat_data[key] = d
        else:
            for gene_name, g in val.items():          # nested: key is the source
                d = parse_gene(g)
                if d and gene_name not in pharmcat_data:
                    pharmcat_data[gene_name] = d
elif isinstance(genes, list):
    for entry in genes:
        gene = entry.get('geneSymbol', entry.get('gene', ''))
        d = parse_gene(entry)
        if gene and d and gene not in pharmcat_data:
            pharmcat_data[gene] = d

if not pharmcat_data:
    print(f'ERROR: parsed 0 genes from PharmCAT report {pharmcat_path}; refusing to write an empty comparison', file=sys.stderr)
    sys.exit(1)
print(f'PharmCAT genes parsed: {len(pharmcat_data)}')

# PharmCAT names VKORC1 alleles 'rs9923231 variant (T)' where pypgx writes
# 'rs9923231', and either tool may put the two alleles in either order. Compare
# the sorted allele names without that suffix; the TSV keeps the raw strings.
def norm(diplotype):
    return sorted(re.sub(r' (variant|reference) \([ACGT]+\)', '', a).strip() for a in diplotype.split('/'))

# Build comparison for overlapping genes
all_genes = sorted(set(list(pypgx_data.keys()) + list(pharmcat_data.keys())))

with open(comparison_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['Gene', 'PharmCAT_diplotype', 'pypgx_diplotype', 'Match', 'Called_by'])
    matches = 0
    mismatches = 0
    for gene in all_genes:
        pc = pharmcat_data.get(gene, 'Not called')
        pg = pypgx_data.get(gene, 'Not called')
        # PharmCAT writes Unknown/Unknown when it could not call a gene; pypgx
        # writes FAILED. Neither is a call, so neither can conflict.
        pc_called = pc != 'Not called' and any(x != 'Unknown' for x in pc.split('/'))
        pg_called = pg not in ('Not called', 'FAILED')
        if not pc_called and not pg_called:
            continue
        if not pc_called:
            match = called_by = 'pypgx only'
            mismatches += 1
        elif not pg_called:
            match = called_by = 'PharmCAT only'
            mismatches += 1
        elif norm(pc) == norm(pg):
            match, called_by = 'Yes', 'both'
            matches += 1
        else:
            match, called_by = 'No', 'both'
            mismatches += 1
        w.writerow([gene, pc, pg, match, called_by])

print(f'Comparison written: {comparison_path}')
print(f'Concordant: {matches}, Discordant/partial: {mismatches}')
" 2>&1
else
  echo ""
  echo "NOTE: No PharmCAT output found. Run step 7 first if you want a comparison."
  echo "  Expected in: ${GENOME_DIR}/${SAMPLE}/pharmcat/ or ${GENOME_DIR}/${SAMPLE}/vcf/"
fi

# Print summary
echo ""
echo "============================================"
echo "  pypgx complete: ${SAMPLE}"
echo "============================================"
echo "Results:    ${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
if [ -n "$PHARMCAT_JSON" ]; then
  echo "Comparison: ${OUTPUT_DIR}/${SAMPLE}_pharmcat_comparison.tsv"
fi
echo ""
if [ -f "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv" ]; then
  echo "Summary:"
  column -t -s $'\t' "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv" 2>/dev/null || cat "${OUTPUT_DIR}/${SAMPLE}_pypgx_summary.tsv"
fi
