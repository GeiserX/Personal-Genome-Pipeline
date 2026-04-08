#!/usr/bin/env bash
# pypgx — Comprehensive pharmacogenomic star allele calling with SV detection
# Input: BAM + VCF from alignment/variant calling steps
# Output: Per-gene star allele calls, consolidated summary TSV, PharmCAT comparison TSV
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=../versions.env
. "${SCRIPT_DIR}/../versions.env"

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

# Validate sample name to prevent shell injection in bash -c / python3 -c strings
if [[ "$SAMPLE" =~ [^a-zA-Z0-9._-] ]]; then
  echo "ERROR: Sample name contains invalid characters. Use only a-z, A-Z, 0-9, ., _, -" >&2
  exit 1
fi
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
  echo "  cd ${GENOME_DIR}/reference && git clone --branch 0.26.0 --depth 1 https://github.com/sbslee/pypgx-bundle.git" >&2
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
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  -v "${PYPGX_BUNDLE}:/root/pypgx-bundle:ro" \
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

    # Phase 1: Prepare depth of coverage for all SV genes (one-time, from BAM)
    echo "--- Preparing depth of coverage for SV genes ---"
    if ! pypgx prepare-depth-of-coverage \
      "$DOC" "$BAM" --assembly GRCh38 2>&1; then
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

    # Phase 3a: BAM-based genes — SV detection via read depth
    # Omits --variants to avoid pseudogene-confounded VCF calls (CYP2D6/CYP2D7)
    if [ -n "$DOC" ]; then
      for GENE in $BAM_GENES; do
        echo "--- Calling ${GENE} (read-depth SV detection) ---"
        EXTRA=""
        [ -f "$CTRL" ] && EXTRA="--control-statistics $CTRL"
        pypgx run-ngs-pipeline "$GENE" "${OUTBASE}/${GENE}" \
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
docker run --rm --user root \
  --cpus 2 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  -v "${PYPGX_BUNDLE}:/root/pypgx-bundle:ro" \
  "${PYPGX_IMAGE}" \
  python3 -c "
import os, sys, zipfile, csv, io

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
    try:
        import subprocess
        # pypgx print-data results.zip outputs a TSV with columns:
        # [sample] Genotype Phenotype Haplotype1 Haplotype2 ...
        out = subprocess.run(
            ['pypgx', 'print-data', results_zip],
            capture_output=True, text=True
        )
        if out.returncode == 0:
            lines = out.stdout.rstrip().split('\n')
            if len(lines) >= 2:
                headers = lines[0].split('\t')
                values = lines[1].split('\t')
                if 'Genotype' in headers:
                    idx = headers.index('Genotype')
                    if idx < len(values):
                        diplotype = values[idx]
                if 'Phenotype' in headers:
                    idx = headers.index('Phenotype')
                    if idx < len(values):
                        phenotype = values[idx]
    except Exception as e:
        print(f'WARNING: Error extracting {gene}: {e}', file=sys.stderr)

    # SV detection only meaningful for BAM-based genes (read-depth analysis)
    # '*5' = gene deletion in CYP2D6 context; 'x2'/'x3' = duplications
    source = 'BAM' if gene in bam_genes else 'VCF'
    if gene in bam_genes:
        sv_detected = 'Yes' if any(x in (diplotype or '') for x in ['DEL', 'DUP', 'x2', 'x3', '*5']) else 'No'
    else:
        sv_detected = 'N/A'
    rows.append([gene, diplotype, phenotype, sv_detected, source])

with open(summary_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Gene', 'Diplotype', 'Phenotype', 'SV_detected', 'Source'])
    w.writerows(rows)

print(f'Summary written: {summary_path}')
print(f'Genes called: {sum(1 for r in rows if r[1] != \"FAILED\")}/{len(rows)}')
" 2>&1

# Cross-reference with PharmCAT if output exists
PHARMCAT_JSON=""
for DIR in "${GENOME_DIR}/${SAMPLE}/pharmcat" "${GENOME_DIR}/${SAMPLE}/vcf"; do
  for FILE in "${DIR}"/*.report.json "${DIR}"/*_pharmcat.json; do
    if [ -f "$FILE" ] 2>/dev/null; then
      PHARMCAT_JSON="$FILE"
      break 2
    fi
  done
done

if [ -n "$PHARMCAT_JSON" ]; then
  echo ""
  echo "PharmCAT output found, generating comparison..."

  docker run --rm --user root \
    --cpus 2 --memory 4g \
    -v "${GENOME_DIR}:/genome" \
    "${PYTHON_IMAGE}" \
    python3 -c "
import json, csv, os, sys

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

# Load PharmCAT results
pharmcat_data = {}
pharmcat_path = '$(echo "$PHARMCAT_JSON" | sed "s|${GENOME_DIR}|/genome|")'
try:
    with open(pharmcat_path) as f:
        data = json.load(f)

    if 'genes' in data and isinstance(data['genes'], dict):
        for source, gene_dict in data['genes'].items():
            if not isinstance(gene_dict, dict):
                continue
            for gene_name, g in gene_dict.items():
                dips = g.get('sourceDiplotypes', [])
                if not dips:
                    continue
                dip = dips[0]
                a1_obj = dip.get('allele1')
                a2_obj = dip.get('allele2')
                a1 = a1_obj.get('name', '?') if a1_obj else '?'
                a2 = a2_obj.get('name', '?') if a2_obj else '?'
                pharmcat_data[gene_name] = f'{a1}/{a2}'
            break
    elif 'genes' in data and isinstance(data['genes'], list):
        for entry in data['genes']:
            gene = entry.get('geneSymbol', entry.get('gene', ''))
            src_dips = entry.get('sourceDiplotypes', [])
            if src_dips:
                pharmcat_data[gene] = src_dips[0].get('label', 'N/A')
except Exception as e:
    print(f'WARNING: Could not parse PharmCAT JSON: {e}', file=sys.stderr)

# Build comparison for overlapping genes
all_genes = sorted(set(list(pypgx_data.keys()) + list(pharmcat_data.keys())))

with open(comparison_path, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Gene', 'PharmCAT_diplotype', 'pypgx_diplotype', 'Match'])
    matches = 0
    mismatches = 0
    for gene in all_genes:
        pc = pharmcat_data.get(gene, 'Not called')
        pg = pypgx_data.get(gene, 'Not called')
        if pc == 'Not called' and pg == 'Not called':
            continue
        if pc == pg:
            match = 'Yes'
            matches += 1
        elif pc == 'Not called':
            match = 'pypgx only'
            mismatches += 1
        elif pg == 'Not called' or pg == 'FAILED':
            match = 'PharmCAT only'
            mismatches += 1
        else:
            match = 'No'
            mismatches += 1
        w.writerow([gene, pc, pg, match])

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
