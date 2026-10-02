# Step 25: Polygenic Risk Scores (PRS)

## What This Does

Calculates polygenic risk scores for 9 common conditions using validated scoring files from the PGS Catalog and plink2. Each PRS aggregates the tiny effects of hundreds to millions of genetic variants into a single number representing your relative genetic predisposition for a trait or disease.

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

Each label is the `trait_reported` value the [PGS Catalog REST API](https://www.pgscatalog.org/rest/) returns for that score, and the script prints the same label.

| Condition (PGS Catalog trait) | PGS ID | Variants | Publication |
|---|---|---|---|
| Coronary artery disease | [PGS000018](https://www.pgscatalog.org/score/PGS000018/) | 1,745,179 | Inouye et al. 2018, J Am Coll Cardiol |
| Type 2 diabetes (T2D) | [PGS000014](https://www.pgscatalog.org/score/PGS000014/) | 6,917,436 | Khera et al. 2018, Nat Genet |
| Breast cancer | [PGS000004](https://www.pgscatalog.org/score/PGS000004/) | 313 | Mavaddat et al. 2018, Am J Hum Genet |
| Prostate cancer | [PGS000662](https://www.pgscatalog.org/score/PGS000662/) | 269 | Conti et al. 2021, Nat Genet |
| Atrial fibrillation | [PGS000016](https://www.pgscatalog.org/score/PGS000016/) | 6,730,541 | Khera et al. 2018, Nat Genet |
| Late-onset Alzheimer’s disease | [PGS000334](https://www.pgscatalog.org/score/PGS000334/) | 22 | Zhang et al. 2020, Nat Commun |
| Body mass index (BMI) | [PGS000027](https://www.pgscatalog.org/score/PGS000027/) | 2,100,302 | Khera et al. 2019, Cell |
| Inflammatory bowel disease | [PGS000017](https://www.pgscatalog.org/score/PGS000017/) | 6,907,112 | Khera et al. 2018, Nat Genet |
| Colorectal cancer | [PGS000055](https://www.pgscatalog.org/score/PGS000055/) | 76 | Schmit et al. 2019, J Natl Cancer Inst |

There is no schizophrenia score yet. An earlier version listed PGS000738 as schizophrenia, but that score is for vitiligo; a schizophrenia row comes back only once a score is chosen from the catalog and checked against the API.

To check a label before adding a score:

```bash
curl -s https://www.pgscatalog.org/rest/score/PGS000017 | jq -r '.trait_reported, .variants_number'
```

## What the Script Does Internally

1. Downloads the GRCh38-harmonized scoring file of each score from the PGS Catalog FTP (one-time, cached in `${GENOME_DIR}/prs_scores/`). The download goes to a `.part` file first, and the file is kept only if its `#HmPOS_build` header says `GRCh38`. If the download fails or the build is anything else, the step stops with an error. There is no fallback to the author-reported file, which is often GRCh37 or rsID-only and would score the wrong positions without any visible sign.
2. Converts your VCF to plink2 binary format (pgen/pvar/psam), restricting to autosomes (chr1-22) and assigning variant IDs in `chr:pos` format (matching PGS Catalog convention)
3. For each scoring file, reformats the harmonized PGS Catalog columns (`hm_chr`, `hm_pos`, effect allele, weight) into plink2's `--score` input format, deduplicating entries with the same variant ID and allele. Rows the catalog could not map to GRCh38 have no `hm_pos` and are dropped
4. Deletes any `.sscore` left by an earlier run, then runs `plink2 --score ... cols=+scoresums` for each condition. A plink2 failure stops the step. The one exception is a score with no variant at all in your VCF, which is reported as `NA` with 0 matched
5. Collects all results into a summary TSV

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_prs_summary.tsv` | Tab-delimited summary: `Condition`, `PGS_ID`, `Score_SUM`, `Variants_Matched`, `Variants_Total` |
| `${PGS_ID}.sscore` | Raw plink2 score output per condition |
| `${PGS_ID}_formatted.tsv` | Reformatted scoring file used for each calculation |
| `${SAMPLE}.pgen/.pvar/.psam` | plink2 binary genotype files (intermediate) |

All output is written to `${GENOME_DIR}/${SAMPLE}/prs/`.

## Runtime

~20-40 minutes total (dominated by VCF-to-plink conversion and scoring across all 9 conditions).

## Interpreting Results

The summary TSV contains a raw score for each condition. Here is what the columns mean:

- **Score_SUM**: Weighted sum of the effect alleles you carry (plink2's `SCORE1_SUM` column). Higher = more genetic predisposition.
- **Variants_Matched**: How many scoring variants were found in your VCF (plink2's `ALLELE_CT` divided by 2).
- **Variants_Total**: Variants in the scoring file with a GRCh38 position.

### The score is biased until the pipeline keeps hom-ref sites

The VCF from step 3 lists only sites where you differ from the reference. A scoring variant whose effect allele is the reference allele is therefore missing from the VCF when you carry two copies of it, and it adds nothing to your sum. So the sum misses the weight of every reference-allele effect allele you carry on two copies, and `Variants_Matched` counts only the sites present in the VCF. The step prints this line under every score:

```
hom-ref sites are absent from this VCF, so the score is biased; not comparable to published distributions
```

Scoring from a gVCF, which records hom-ref sites, removes the bias.

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

Check the `Variants_Matched / Variants_Total` ratio. If fewer than 50% of scoring variants matched, the score is less reliable. Low matching rates usually indicate:
- The scoring file was built on array data with different variant coverage than WGS
- Variant ID format mismatches between your VCF and the scoring file

## Limitations

- PRS were predominantly developed in European-ancestry populations. They are less accurate for other ancestries.
- A PRS captures only the genetic component. Lifestyle, environment, and family history are often more predictive.
- Sex-specific conditions (breast cancer, prostate cancer) should be interpreted accordingly.
- Scoring file availability and quality vary. If the PGS Catalog FTP is unavailable the step stops; rerun it later (files already downloaded stay cached).
- No mean imputation is used (`no-mean-imputation` flag), so missing variants reduce the score proportionally rather than being imputed to population averages.

## Notes

- Scoring files are downloaded once and cached in `${GENOME_DIR}/prs_scores/`. Delete this directory to force re-download.
- The script uses only GRCh38-harmonized scoring files. A cached file without a `#HmPOS_build=GRCh38` header (for example one written by an older version of this script) is deleted and downloaded again.
- You can add more scores by adding a `"<PGS ID>|<trait_reported>"` line to the `PGS_SCORES` list in the script, with the label copied from the API. Browse available scores at [pgscatalog.org](https://www.pgscatalog.org/).

## Maintenance

- Recheck the PGS Catalog against its latest release page at least quarterly before treating this step as "current."
- A scoring file update is a **result-changing event**. If the harmonized file version/date changes, rerun step 25 and treat the output as a new baseline.
- If you publish or compare PRS results over time, keep the `PGS ID`, the harmonized scoring file version/date, and the pipeline commit together so score changes remain auditable.

## Links

- [PGS Catalog](https://www.pgscatalog.org/)
- [plink2 documentation](https://www.cog-genomics.org/plink/2.0/)
- [PGS Catalog scoring file format](https://www.pgscatalog.org/downloads/#scoring_files)
- [Khera et al. 2018 (multi-trait PRS)](https://doi.org/10.1038/s41588-018-0183-z)
