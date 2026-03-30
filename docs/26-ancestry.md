# Step 26: Ancestry Estimation via PCA

## What This Does

Runs principal component analysis (PCA) on your sample using common SNPs shared with the 1000 Genomes Project reference panel. Because this is single-sample PCA (not a joint projection with the 1000G cohort), the resulting PC values capture your genome's internal variance structure but are **not directly comparable** to published population cluster plots. See [Single-sample limitation](#single-sample-limitation) below for details.

## Why

Knowing your genetic ancestry is useful for two practical reasons:

1. **PRS interpretation**: Polygenic risk scores (step 25) are ancestry-dependent. Knowing where you fall on the ancestry spectrum helps contextualize your scores.
2. **Variant filtering**: Some variants are common in one population but rare in another. Ancestry helps distinguish benign population-specific variants from truly rare findings.

Note: This step produces single-sample PCA, which is a starting point but cannot place you on a population map without joint analysis. See [Interpreting Results](#interpreting-results) for details.

## Tool

- **plink2** for PCA computation and LD pruning
- **bcftools** for variant intersection and filtering
- **1000 Genomes Project** Phase 3 as the reference panel

## Docker Images

```
pgscatalog/plink2:2.00a5.10
staphb/bcftools:1.21
```

## Input

- VCF from DeepVariant (step 3): `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz`

## Command

```bash
./scripts/26-ancestry.sh your_name
```

## What the Script Does Internally

1. **Downloads 1000 Genomes reference SNPs** (one-time, ~100 MB): fetches the GRCh38 biallelic SNV sites file and filters to common autosomal SNPs (MAF 5-95%)
2. **Downloads population labels**: maps each 1000G sample to its super-population (AFR, AMR, EAS, EUR, SAS)
3. **Intersects your VCF with the reference**: finds SNPs present in both your sample and the 1000G panel using `bcftools isec`
4. **LD prunes**: removes correlated SNPs (window 50, step 5, r-squared threshold 0.2) to avoid redundant signal. PCA requires independent markers.
5. **Runs PCA**: computes the first 10 principal components from the LD-pruned SNP set

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_pca.eigenvec` | Principal component values (10 PCs per sample) |
| `${SAMPLE}_pca.eigenval` | Eigenvalues showing variance explained by each PC |
| `${SAMPLE}_shared.vcf.gz` | SNPs shared between your sample and 1000G |
| `${SAMPLE}_ld.prune.in` | SNPs retained after LD pruning |
| `${SAMPLE}_ld.prune.out` | SNPs removed by LD pruning |

All output is written to `${GENOME_DIR}/${SAMPLE}/ancestry/`. Reference data is cached in `${GENOME_DIR}/ancestry_ref/`.

## Runtime

~15-30 minutes (dominated by the initial 1000G download on first run; subsequent runs are faster).

## Interpreting Results

The `eigenvec` file contains your sample's coordinates on 10 principal components. The `eigenval` file shows how much variance each PC explains.

### Single-sample limitation

This script runs PCA on **your sample alone**, not jointly with the 1000G reference panel. This is a fundamental limitation: in population-structure PCA (Price et al. 2006), the PC axes are defined by the variance across many individuals. With a single sample, the axes instead capture internal genotype variance (e.g., heterozygosity patterns), which does not map onto population-level structure.

The PC values from this step are **not comparable** to published 1000G PCA plots, where PC1 separates African from non-African ancestry and PC2 separates European from East Asian. Those axis interpretations require joint PCA across a multi-population cohort.

To properly place yourself on a population map, you would need to:

1. Download the full 1000G genotype data (~30-50 GB)
2. Merge your sample with the 1000G samples
3. Run joint PCA on the combined dataset
4. Plot your sample against the 1000G population clusters

This pipeline does not perform joint PCA. The single-sample output is included as a starting point for users who want to extend it with their own reference panel.

## Limitations

- Single-sample PCA cannot produce population percentages (e.g., "85% European, 15% other"). That requires admixture analysis tools like ADMIXTURE or RFMix with a reference panel.
- The 1000G panel does not represent all global populations equally. Fine-grained ancestry (e.g., distinguishing Spanish from Italian) requires specialized reference panels.
- Low variant overlap between your VCF and the reference panel weakens results. The script warns if fewer than 1,000 shared SNPs are found.
- The reference sites download URL from the 1000 Genomes FTP may occasionally be unavailable.

## Notes

- Reference data (1000G SNPs and population labels) is downloaded once and cached in `${GENOME_DIR}/ancestry_ref/`. Delete this directory to force re-download.
- LD pruning parameters (window=50, step=5, r2=0.2) are standard for ancestry PCA.
- 10 PCs are computed by default. For single-sample PCA this is more than sufficient; additional PCs would not add interpretable signal without a reference cohort.
- For a more complete ancestry analysis, consider uploading your VCF to tools like [Gnomix](https://github.com/AI-sandbox/gnomix) or using the PLINK `--admixture` approach.

## Links

- [1000 Genomes Project](https://www.internationalgenome.org/)
- [plink2 PCA documentation](https://www.cog-genomics.org/plink/2.0/strat)
- [1000G data portal (GRCh38)](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)
- [Price et al. 2006 (PCA for population structure)](https://doi.org/10.1038/ng1847)
