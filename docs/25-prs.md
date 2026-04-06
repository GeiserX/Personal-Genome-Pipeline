# Step 25: Polygenic Risk Scores (PRS)

## What This Does

Calculates polygenic risk scores for 10 common conditions using validated scoring files from the PGS Catalog and plink2. Each PRS aggregates the tiny effects of hundreds to millions of genetic variants into a single number representing your relative genetic predisposition for a trait or disease.

## Why

Most common diseases (heart disease, diabetes, cancer) are not caused by a single gene. They result from the combined effect of many variants, each contributing a small amount of risk. A PRS sums these contributions using weights derived from large genome-wide association studies (GWAS). While no single variant is predictive on its own, the aggregate score can be informative.

## Tool

- **plink2** (Chang et al., GigaScience 2015) -- the standard tool for large-scale genomic computation
- **PGS Catalog** -- curated repository of published polygenic scoring files

## Docker Image

```
pgscatalog/plink2:2.00a5.10
```

## Input

- VCF from DeepVariant (step 3): `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz`

## Command

```bash
./scripts/25-prs.sh your_name
```

## Conditions Scored

| Condition | PGS ID | Source |
|---|---|---|
| Coronary artery disease | PGS000018 | Khera et al. 2018 |
| Type 2 diabetes | PGS000014 | Mahajan et al. 2018 |
| Breast cancer | PGS000004 | Mavaddat et al. 2019 |
| Prostate cancer | PGS000662 | Conti et al. 2021 |
| Atrial fibrillation | PGS000016 | Khera et al. 2018 |
| Alzheimer's disease | PGS000334 | De Rojas et al. 2021 |
| Body mass index | PGS000027 | Khera et al. 2019 |
| Schizophrenia | PGS000738 | PGC 2022 |
| Inflammatory bowel disease | PGS000020 | Khera et al. 2018 |
| Colorectal cancer | PGS000055 | Huyghe et al. 2019 |

## What the Script Does Internally

1. Downloads GRCh38-harmonized scoring files from the PGS Catalog FTP (one-time, cached in `${GENOME_DIR}/prs_scores/`)
2. Converts your VCF to plink2 binary format (pgen/pvar/psam), restricting to autosomes (chr1-22) and assigning variant IDs in `chr:pos` format (matching PGS Catalog convention)
3. For each scoring file, reformats the PGS Catalog columns (chromosome, position, effect allele, weight) into plink2's `--score` input format, deduplicating entries with the same variant ID and allele
4. Runs `plink2 --score` for each condition, producing a `.sscore` file with the aggregate score and the number of variants matched
5. Collects all results into a summary TSV

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_prs_summary.tsv` | Tab-delimited summary: condition, PGS ID, score, variants used, variants total |
| `${PGS_ID}.sscore` | Raw plink2 score output per condition |
| `${PGS_ID}_formatted.tsv` | Reformatted scoring file used for each calculation |
| `${SAMPLE}.pgen/.pvar/.psam` | plink2 binary genotype files (intermediate) |

All output is written to `${GENOME_DIR}/${SAMPLE}/prs/`.

## Runtime

~20-40 minutes total (dominated by VCF-to-plink conversion and scoring across all 10 conditions).

## Interpreting Results

The summary TSV contains a raw score for each condition. Here is what the columns mean:

- **Score**: Weighted sum of risk alleles you carry. Higher = more genetic predisposition.
- **Variants_Used**: How many scoring variants matched your VCF.
- **Variants_Total**: Total variants in the scoring file.

### What these scores are NOT

- They are NOT percentiles. A raw score of 0.5 does not mean 50th percentile.
- They are NOT probabilities. A high score does not mean you will develop the condition.
- They are NOT comparable across conditions. A score of 10 for CAD and 10 for T2D mean entirely different things.
- They are NOT stable across arbitrary pipeline changes. If you change the PGS file version, genome build harmonization, or variant matching rules, you need to recompute and reinterpret the score.

### How to make them meaningful

Raw PRS become useful only when compared against a population distribution. To convert your score into a percentile, you need a reference panel of thousands of individuals with scores computed using the same scoring file. The PGS Catalog provides some population-level statistics, but full percentile calculation requires a reference cohort (not included in this pipeline).

Comparing two people is only defensible when both were scored with the same PGS ID, the same scoring file version, the same genome build conventions, and the same preprocessing. Even then, treat the comparison as directional rather than clinically calibrated unless you also have a matched reference distribution.

**Do not convert raw scores to percentiles using generic SD thresholds.** The mapping between a raw score and a population percentile depends on the score distribution in a matched reference cohort (same ancestry, same scoring file, same preprocessing). Without that cohort, statements like "top 16%" or "top 2.5%" are not grounded. See the [PGS Catalog Calculator interpretation guide](https://pgsc-calc.readthedocs.io/) and the ACMG points-to-consider for PRS reporting.

### Variant matching

Check the `Variants_Used / Variants_Total` ratio. If fewer than 50% of scoring variants matched, the score is less reliable. Low matching rates usually indicate:
- The scoring file was built on array data with different variant coverage than WGS
- Variant ID format mismatches between your VCF and the scoring file

## Limitations

- PRS were predominantly developed in European-ancestry populations. They are less accurate for other ancestries.
- A PRS captures only the genetic component. Lifestyle, environment, and family history are often more predictive.
- Sex-specific conditions (breast cancer, prostate cancer) should be interpreted accordingly.
- Scoring file availability and quality vary. Some PGS IDs may fail to download if the PGS Catalog FTP is unavailable.
- No mean imputation is used (`no-mean-imputation` flag), so missing variants reduce the score proportionally rather than being imputed to population averages.

## Notes

- Scoring files are downloaded once and cached in `${GENOME_DIR}/prs_scores/`. Delete this directory to force re-download.
- The script prefers GRCh38-harmonized scoring files. If unavailable, it falls back to the original (which may be on GRCh37 and produce poor variant matching).
- You can add more PGS IDs by editing the `PGS_IDS` associative array in the script. Browse available scores at [pgscatalog.org](https://www.pgscatalog.org/).

## Maintenance

- Recheck the PGS Catalog against its latest release page at least quarterly before treating this step as "current."
- A scoring file update is a **result-changing event**. If the harmonized file version/date changes, rerun step 25 and treat the output as a new baseline.
- If you publish or compare PRS results over time, keep the `PGS ID`, the harmonized scoring file version/date, and the pipeline commit together so score changes remain auditable.

## Links

- [PGS Catalog](https://www.pgscatalog.org/)
- [plink2 documentation](https://www.cog-genomics.org/plink/2.0/)
- [PGS Catalog scoring file format](https://www.pgscatalog.org/downloads/#scoring_files)
- [Khera et al. 2018 (multi-trait PRS)](https://doi.org/10.1038/s41588-018-0183-z)
