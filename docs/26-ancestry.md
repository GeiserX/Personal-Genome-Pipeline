# Step 26: Ancestry SNP Intersection [EXPERIMENTAL]

## What This Does

Intersects your sample's common SNPs with the 1000 Genomes Project reference panel and runs LD pruning. The script attempts single-sample PCA with plink2, but **single-sample PCA is mathematically degenerate** — PCA defines axes from variance across a cohort (Price et al. 2006), so one sample cannot produce interpretable principal components. The output is best understood as a prepared SNP set for users who want to extend it with joint multi-sample PCA against a reference panel.

## Why

The intermediate outputs (shared SNPs, LD-pruned variant set) are useful for two reasons:

1. **PRS interpretation**: Polygenic risk scores (step 25) are ancestry-dependent. The ancestry SNP set helps identify which population reference to use.
2. **Variant filtering**: Population-specific variant frequencies help distinguish benign variants from truly rare findings.

**This step does NOT produce a usable ancestry estimate.** For ancestry analysis from WGS, you need joint PCA or admixture analysis against a multi-population reference panel (not implemented here). For a quick ancestry check, consumer services (23andMe, AncestryDNA) or tools like [Gnomix](https://github.com/AI-sandbox/gnomix) with a reference cohort are more appropriate.

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

`run-all.sh` does not run this step unless you ask for it, because on one sample it downloads about 900 MB to produce a count:

```bash
ANCESTRY=true ./scripts/run-all.sh your_name <male|female>
```

## What the Script Does Internally

1. **Downloads 1000 Genomes reference SNPs** (one-time, ~900 MB): fetches the GRCh38 biallelic SNV sites file and filters to common autosomal SNPs (MAF 5-95%)
2. **Downloads population labels**: maps each 1000G sample to its super-population (AFR, AMR, EAS, EUR, SAS)
3. **Intersects your VCF with the reference**: finds SNPs present in both your sample and the 1000G panel using `bcftools isec`, and prints how many there are
4. **Tries LD pruning** (window 50, step 5, r-squared threshold 0.2). plink2 needs at least 50 samples for this, so on one sample it fails and the script carries on with all shared SNPs.
5. **Tries PCA**. plink2 needs at least 2 samples, so on one sample it fails too and the script says so.

On a single sample the only result is the shared SNP set and its count.

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_shared.vcf.gz` (+ `.tbi`) | SNPs shared between your sample and 1000G |

The `${SAMPLE}_ld.prune.in`/`.prune.out` and `${SAMPLE}_pca.eigenvec`/`.eigenval` files are only written when plink2 gets enough samples, which never happens with one sample.

All output is written to `${GENOME_DIR}/${SAMPLE}/ancestry/`. Reference data is cached in `${GENOME_DIR}/ancestry_ref/`.

## Runtime

~15-30 minutes (dominated by the initial 1000G download on first run; subsequent runs are faster).

## Interpreting Results

The step reports how many of your SNPs are also common SNPs in the 1000G panel. A count below 1,000 usually means a different genome build or a VCF with few variants. It does not tell you anything about your ancestry by itself.

### Single-sample limitation

This script can only attempt PCA on **your sample alone**, not jointly with the 1000G reference panel, and plink2 refuses it. Even if it ran, this is a fundamental limitation: in population-structure PCA (Price et al. 2006), the PC axes are defined by the variance across many individuals. With a single sample, the axes instead capture internal genotype variance (e.g., heterozygosity patterns), which does not map onto population-level structure.

Single-sample PC values would **not be comparable** to published 1000G PCA plots, where PC1 separates African from non-African ancestry and PC2 separates European from East Asian. Those axis interpretations require joint PCA across a multi-population cohort.

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
- The script asks plink2 for 10 PCs. With one sample plink2 computes none; with a reference cohort merged in, 10 is the usual number.
- For a more complete ancestry analysis, consider running a local-ancestry tool such as [Gnomix](https://github.com/AI-sandbox/gnomix) on your VCF (it runs on your own machine with a reference panel), or ADMIXTURE.

## Links

- [1000 Genomes Project](https://www.internationalgenome.org/)
- [plink2 PCA documentation](https://www.cog-genomics.org/plink/2.0/strat)
- [1000G data portal (GRCh38)](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)
- [Price et al. 2006 (PCA for population structure)](https://doi.org/10.1038/ng1847)
