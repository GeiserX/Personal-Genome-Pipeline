#!/usr/bin/env bash
# 26-ancestry.sh — [EXPERIMENTAL] Ancestry-informative PCA from WGS VCF
# Usage: ./scripts/26-ancestry.sh <sample_name>
#
# EXPERIMENTAL: This runs PCA on a SINGLE sample intersected with 1000 Genomes
# reference sites. Without joint PCA against the actual 1000G genotypes, the
# resulting principal components are NOT directly interpretable for ancestry.
# The PCs only become meaningful when projected alongside reference samples.
#
# This step is included as a starting point. For proper ancestry estimation,
# you need multi-sample PCA with a reference panel (not implemented here).
#
# Requires: VCF from step 3, reference panel (downloaded automatically)
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}

VCF="${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
OUTDIR="${GENOME_DIR}/${SAMPLE}/ancestry"
REFDIR="${GENOME_DIR}/ancestry_ref"
mkdir -p "$OUTDIR" "$REFDIR"

if [ ! -f "$VCF" ]; then
  echo "ERROR: VCF not found: ${VCF}"
  echo "  Run step 3 (DeepVariant) first."
  exit 1
fi

echo "============================================"
echo "  Step 26: Ancestry PCA"
echo "  Tool: plink2"
echo "  Sample: ${SAMPLE}"
echo "  Input:  ${VCF}"
echo "  Output: ${OUTDIR}/"
echo "============================================"
echo ""

# Step 1: Download 1000 Genomes reference sites (if not present)
# Using a pre-computed set of common, LD-pruned SNPs from 1000G Phase 3
KG_SITES="${REFDIR}/1kg_common_snps.vcf.gz"
KG_POPS="${REFDIR}/1kg_populations.tsv"

if [ ! -f "$KG_SITES" ]; then
  echo "[1/5] Downloading 1000 Genomes reference SNPs..."
  echo "  This is a one-time download (~1 GB)."

  RAW_VCF="${REFDIR}/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz"
  # Use a lock file to prevent parallel runs from downloading simultaneously
  LOCKFILE="${REFDIR}/.download_lock"
  if [ -f "$LOCKFILE" ]; then
    echo "  Another download in progress (lock file exists). Waiting..."
    while [ -f "$LOCKFILE" ]; do sleep 5; done
  fi

  if [ ! -f "$KG_SITES" ]; then
    touch "$LOCKFILE"
    trap 'rm -f "$LOCKFILE"' EXIT

    if [ ! -f "$RAW_VCF" ]; then
      wget -q -O "$RAW_VCF" \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz" 2>/dev/null || {
        echo "WARNING: Could not download 1000G reference sites."
        echo "  The 1000 Genomes FTP may be temporarily unavailable."
        echo "  Try again later, or download manually from:"
        echo "  https://www.internationalgenome.org/data-portal/data-collection/30x-grch38"
        rm -f "$LOCKFILE"
        exit 1
      }
    fi

    # Index the raw file first (required by bcftools --regions)
    if [ ! -f "${RAW_VCF}.tbi" ]; then
      docker run --rm --user root \
        -v "${GENOME_DIR}:/genome" \
        staphb/bcftools:1.21 \
        bcftools index -t "/genome/ancestry_ref/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz"
    fi

    # Extract common biallelic SNPs (MAF > 5%, autosomal only)
    echo "  Filtering to common biallelic autosomal SNPs..."
    docker run --rm --user root \
      --cpus 4 --memory 4g \
      -v "${GENOME_DIR}:/genome" \
      staphb/bcftools:1.21 \
      bash -c "
        bcftools view -m2 -M2 -v snps \
          -i 'AF>=0.05 && AF<=0.95' \
          --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
          /genome/ancestry_ref/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz \
          -Oz -o /genome/ancestry_ref/1kg_common_snps.vcf.gz &&
        bcftools index -t /genome/ancestry_ref/1kg_common_snps.vcf.gz
      "
    rm -f "$LOCKFILE"
    echo "  [OK] Reference SNPs prepared."
  fi
else
  echo "[1/5] 1000G reference SNPs already downloaded."
fi

# Step 2: Download population labels
if [ ! -f "$KG_POPS" ]; then
  echo "[2/5] Downloading population labels..."
  wget -q -O "${REFDIR}/integrated_call_samples_v3.20130502.ALL.panel" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" 2>/dev/null || true

  # Create simple population mapping
  if [ -f "${REFDIR}/integrated_call_samples_v3.20130502.ALL.panel" ]; then
    awk -F'\t' 'NR>1 {print $1"\t"$3}' "${REFDIR}/integrated_call_samples_v3.20130502.ALL.panel" > "$KG_POPS"
  else
    echo "WARNING: Could not download population labels. PCA will run without labels."
    touch "$KG_POPS"
  fi
else
  echo "[2/5] Population labels already present."
fi

# Step 3: Extract overlapping SNPs between your sample and reference
echo "[3/5] Finding shared SNPs between your sample and reference panel..."
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "
    bcftools isec -n=2 -w1 \
      /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
      /genome/ancestry_ref/1kg_common_snps.vcf.gz \
      -Oz -o /genome/${SAMPLE}/ancestry/${SAMPLE}_shared.vcf.gz &&
    bcftools index -t /genome/${SAMPLE}/ancestry/${SAMPLE}_shared.vcf.gz
  "

SHARED_COUNT=$(docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools view -H "/genome/${SAMPLE}/ancestry/${SAMPLE}_shared.vcf.gz" 2>/dev/null | wc -l || echo 0)
echo "  Shared SNPs: ${SHARED_COUNT}"

if [ "$SHARED_COUNT" -lt 1000 ]; then
  echo "WARNING: Very few shared SNPs (${SHARED_COUNT}). PCA results may be unreliable."
  echo "  This usually means your VCF uses different variant IDs or genome build."
fi

# Step 4: LD pruning (requires >=50 samples — skip for single-sample input)
echo "[4/5] LD pruning..."
PRUNED_COUNT=0
if docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  pgscatalog/plink2:2.00a5.10 \
  plink2 \
    --vcf "/genome/${SAMPLE}/ancestry/${SAMPLE}_shared.vcf.gz" \
    --indep-pairwise 50 5 0.2 \
    --out "/genome/${SAMPLE}/ancestry/${SAMPLE}_ld" \
    --threads 4 \
    --memory 6000 \
    --chr 1-22 \
    --allow-extra-chr \
    --output-chr chrM 2>&1; then
  PRUNED_COUNT=$(wc -l < "${OUTDIR}/${SAMPLE}_ld.prune.in" 2>/dev/null || echo 0)
  echo "  LD-pruned SNPs: ${PRUNED_COUNT}"
  EXTRACT_ARGS=(--extract "/genome/${SAMPLE}/ancestry/${SAMPLE}_ld.prune.in")
else
  echo "  LD pruning skipped (requires >=50 samples). Using all shared SNPs."
  EXTRACT_ARGS=()
  PRUNED_COUNT="${SHARED_COUNT}"
fi

# Step 5: PCA
# Note: PCA on a single sample is mathematically degenerate (0 PCs from 1 sample).
# plink2 --pca requires >=2 samples. For proper ancestry estimation, you need
# joint PCA with a multi-sample reference panel (e.g., 1000G genotypes).
# We still record the shared SNP and LD-pruned counts as useful QC metrics.
echo "[5/5] Attempting PCA (requires >=2 samples)..."
if docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v "${GENOME_DIR}:/genome" \
  pgscatalog/plink2:2.00a5.10 \
  plink2 \
    --vcf "/genome/${SAMPLE}/ancestry/${SAMPLE}_shared.vcf.gz" \
    "${EXTRACT_ARGS[@]}" \
    --pca 10 \
    --out "/genome/${SAMPLE}/ancestry/${SAMPLE}_pca" \
    --threads 4 \
    --memory 6000 \
    --chr 1-22 \
    --allow-extra-chr \
    --output-chr chrM 2>&1; then
  echo "  PCA completed."
else
  echo "  PCA could not run (expected for single-sample input)."
  echo "  Single-sample PCA is mathematically undefined (need >=2 samples)."
  echo "  The shared SNP counts and LD-pruned variants above are still useful."
fi

echo ""
echo "============================================"
echo "  Ancestry PCA complete: ${SAMPLE}"
echo ""
echo "  Shared SNPs with 1000G: ${SHARED_COUNT:-N/A}"
echo "  LD-pruned SNPs: ${PRUNED_COUNT:-N/A}"
if [ -f "${OUTDIR}/${SAMPLE}_pca.eigenvec" ]; then
  echo "  PCA results: ${OUTDIR}/${SAMPLE}_pca.eigenvec"
  echo "  Eigenvalues: ${OUTDIR}/${SAMPLE}_pca.eigenval"
else
  echo "  PCA: Not available (single-sample limitation)"
fi
echo "============================================"
echo ""

# Display PC1 and PC2 values
if [ -f "${OUTDIR}/${SAMPLE}_pca.eigenvec" ]; then
  echo "  Principal Components (first 5):"
  echo "  ────────────────────────────────"
  head -2 "${OUTDIR}/${SAMPLE}_pca.eigenvec" | column -t
  echo ""
  echo "  To interpret these results, you need to project your PCs onto a"
  echo "  1000 Genomes reference PCA plot. General guidelines:"
  echo ""
  echo "  - PC1 primarily separates African vs non-African ancestry"
  echo "  - PC2 separates European vs East Asian ancestry"
  echo "  - PC3 separates South Asian ancestry"
  echo "  - PC4+ capture finer-grained population structure"
  echo ""
  echo "  For single-sample PCA (without a reference panel), the absolute"
  echo "  PC values are not directly interpretable. They become meaningful"
  echo "  when projected alongside the 1000G reference samples."
fi

echo ""
echo "NOTE: Single-sample PCA has limited power. For robust ancestry"
echo "estimation, you need multi-sample PCA with a reference panel."
echo "See docs/26-ancestry.md for detailed interpretation."
