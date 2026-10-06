# Step 25: Polygenic Risk Scores (PRS)

## What This Does

Scores you for 9 common conditions with the PGS Catalog's own calculator, [pgsc_calc](https://github.com/PGScatalog/pgsc_calc). Each polygenic score adds up the small effects of hundreds to millions of variants into one number. With the ancestry reference panel installed, each score is also given as a **percentile**: where your score falls among the reference samples whose genetic ancestry is most similar to yours. Without the panel the step reports the raw sum and says "raw score only", because a raw sum from one person cannot be compared with anyone.

## Why

Most common diseases (heart disease, diabetes, cancer) are not caused by a single gene. They come from the combined effect of many variants, each adding a small amount of risk. A score sums those contributions with weights from large genome-wide association studies (GWAS). The sum only means something next to the sums of other people scored the same way, and the distribution of sums differs between ancestries. That is why the percentile is taken among the most similar reference group, and why pgsc_calc is used: it matches variants the way the PGS Catalog intends, and with a reference panel it reports ancestry-adjusted percentiles.

## Tool

- **pgsc_calc** (Lambert et al., Nat Genet 2024), the PGS Catalog Calculator, release `PGSC_CALC_VERSION` in `versions.env`, Apache-2.0. It is a Nextflow pipeline of its own; step 25 starts it.
- **PGS Catalog**: curated repository of published polygenic scoring files.

## Docker Image

pgsc_calc runs its steps in these images, all pinned in `versions.env`: `PGSC_UTILS_IMAGE` (matching and ancestry adjustment), `PLINK2_IMAGE` (scoring), `PGSC_FRAPOSA_IMAGE` (projection onto the panel), `PGSC_ZSTD_IMAGE`, `PGSC_PYYAML_IMAGE` and `PGSC_REPORT_IMAGE` (its HTML report). Step 25 also uses `PYTHON_IMAGE` and `BCFTOOLS_IMAGE`. [Image versions](versions.md) lists the tags.

## Input

- VCF from DeepVariant (step 3): `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz`
- gVCF from DeepVariant (step 3), when present: `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz`. The score positions are genotyped from it, so sites where you match the reference count.
- The scores of [`assets/pgs_scores.tsv`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/assets/pgs_scores.tsv), downloaded by `setup.sh` into `${GENOME_DIR}/prs_scores/`.
- Optional: the ancestry reference panel, `${GENOME_DIR}/reference/pgsc_calc/pgsc_1000G_v1.tar.zst` with its site list beside it (`scripts/setup.sh --ancestry-panel <genome_dir>`).
- Java 17+ and Nextflow on the host, as for `run-all.sh`.

## Command

```bash
./scripts/setup.sh /path/to/genome_dir                     # the scores, pgsc_calc and its plugin (once)
./scripts/setup.sh --ancestry-panel /path/to/genome_dir    # optional: percentiles (about 7 GB)
./scripts/25-prs.sh your_name
```

`ANCESTRY_PANEL=none` scores without the panel even when it is installed. `PGSC_MAX_MEMORY` (for example `12.GB`) caps the memory pgsc_calc may use; the default is three quarters of the machine's RAM.

## Conditions Scored

The list is [`assets/pgs_scores.tsv`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/assets/pgs_scores.tsv). Each label is the `trait_reported` value the [PGS Catalog REST API](https://www.pgscatalog.org/rest/) returns for that score; `scripts/ci/check-pgs-labels.sh` compares the two, so an ID cannot be printed under the wrong disease.

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

To add a score, add a line to `assets/pgs_scores.tsv` with the label the API gives:

```bash
curl -s https://www.pgscatalog.org/rest/score/PGS000017 | jq -r '.trait_reported, .variants_number'
```

Only additive scores are accepted (no `dosage_*_weight` columns, no `is_dominant` or `is_recessive` rows); step 25 stops on any other.

## What the Script Does Internally

1. **Scoring files.** Downloads the GRCh38-harmonised file of each score from the PGS Catalog FTP when it is not in `${GENOME_DIR}/prs_scores/` yet, checked against the md5 the catalog publishes beside it. A file whose `#HmPOS_build` header is not `GRCh38` is refused. There is no fallback to the author-reported file, which is often GRCh37 or rsID-only and would score the wrong positions without any visible sign.
2. **Scores as pgsc_calc reads them.** `bin/collect_summary.py prs-format` writes each file as a custom GRCh38 scoring file: `chr_name` and `chr_position` from the harmonised `hm_chr` and `hm_pos`, the effect allele, the other allele (from `other_allele`, else `hm_inferOtherAllele` when it names one allele), the weight, and the catalog's trait as its label. Rows the catalog could not place on GRCh38 are dropped. pgsc_calc then needs no network and no liftover.
3. **Genotypes.**
   - **With step 3's gVCF** (the default since step 3 writes one): every score position, and with the panel every panel SNV, is genotyped from the gVCF. `bcftools convert --gvcf2vcf` turns each reference block over a position into a 0/0 call with the reference base; a position with no coverage (`./.`) or outside every block stays missing. A 0/0 record gets as its ALT the position's first allele (a score's effect or other allele, the panel's ALT) that is not the reference, so pgsc_calc can match it. These genotypes are kept as `prs/pgsc_calc/target.vcf.gz`.
   - **Without a gVCF** (an older run): the variant-only VCF is scored, so every site where you match the reference is missing (see below).
4. **pgsc_calc.** Runs `pgsc_calc` (from `${GENOME_DIR}/tools/pgsc_calc-<release>`, which `setup.sh` or the step itself unpacks from GitHub's archive of the release, checked against `PGSC_CALC_SHA256`) with the images of `versions.env`, offline, its containers without network. With the panel it adds `--run_ancestry`. pgsc_calc matches each score's variants to your genotypes (strand flips, ambiguous A/T and C/G pairs dropped, one best match per variant), scores them with plink2 and, with the panel, projects you onto the panel's principal components and compares your score with the reference group most similar to you. A score that matches under 75% of its variants is dropped by pgsc_calc and gets no sum.
5. **Summary.** `bin/collect_summary.py prs-table` reads pgsc_calc's match summary and scores into `${SAMPLE}_prs_summary.tsv`, and with the panel writes step 26's ancestry table. pgsc_calc's work folder is deleted; its results (its own HTML report, the match log) are kept.

The Nextflow pipeline does the same with four processes: `PRS_PREPARE`, `PRS_SCORE_SITES`, `PRS` (pgsc_calc, which runs on the host because it starts its own containers) and `PRS_SUMMARY`. Pass `--ancestry_ref` for the panel and `--pgsc_calc ${GENOME_DIR}/tools/pgsc_calc-<release>` to run offline.

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_prs_summary.tsv` | `Condition`, `PGS_ID`, `Score_SUM`, `Variants_Matched`, `Variants_Total`, `Matched_Pct`, `Percentile`, `Ancestry_Group`, `Input` |
| `pgsc_calc/target.vcf.gz` | The genotypes pgsc_calc scored (from the gVCF) |
| `pgsc_calc/results/sample/score/` | pgsc_calc's scores, its HTML report `report.html`, and with the panel the ancestry-adjusted scores and the principal components |
| `pgsc_calc/results/sample/match/` | pgsc_calc's match log and summary |
| `pgsc_calc/pgsc_calc.log` | pgsc_calc's console output |
| `../ancestry/${SAMPLE}_ancestry.tsv` | With the panel: step 26's table (population, its probability, principal components) |

All output is written to `${GENOME_DIR}/${SAMPLE}/prs/`.

## Runtime

Without the panel, about 20-40 minutes: most of it is the gVCF pass over every score position (the large scores reach most of the genome) and pgsc_calc's conversion of those genotypes. The panel adds its extraction, the intersection with your genotypes, the panel's PCA and the projection; see the measured numbers below.

## The reference panel on a GitHub-hosted runner

The panel is pgsc_calc's 1000 Genomes database, `pgsc_1000G_v1.tar.zst` (`PGSC_PANEL` in `versions.env`), published by the PGS Catalog at https://ftp.ebi.ac.uk/pub/databases/spot/pgs/resources/. Measured by the E2E case `tests/e2e/prs-3-panel-measure.sh` on `ubuntu-latest` (4 CPUs, 16 GB RAM), with pgsc_calc projecting the PGS Catalog's synthetic genome-wide test target (600 samples) onto the panel and scoring PGS000018:

MEASUREMENT_TABLE

## Interpreting Results

- **Score_SUM**: weighted sum of the effect alleles you carry. Higher = more genetic predisposition, but only relative to other people scored the same way.
- **Variants_Matched**: score variants pgsc_calc matched in your genotypes.
- **Variants_Total**: variants of the scoring file with a GRCh38 position.
- **Matched_Pct**: `Variants_Matched / Variants_Total` as a percentage; the step warns below 50%, and pgsc_calc gives no sum below 75%.
- **Percentile**: with the panel, where your score falls among the reference samples of `Ancestry_Group` (pgsc_calc's empirical percentile, `percentile_MostSimilarPop`). `NA` without the panel.
- **Ancestry_Group**: the panel population whose genetic ancestry is most similar to yours (for the 1000 Genomes panel: AFR, AMR, EAS, EUR or SAS). `NA` without the panel.
- **Input**: `gvcf` when the positions were genotyped from the gVCF, `vcf` when only the variant-only VCF was there.

Both reports show the percentile with its group, or one line saying "Raw score only" when no panel was used.

### Hom-ref sites come from the gVCF

The VCF from step 3 lists only sites where you differ from the reference. A scoring variant whose effect allele is the reference allele is missing from that VCF when you carry two copies of it, and adds nothing to your sum. Step 3 also writes a gVCF, which records where you match the reference, and step 25 reads those sites from it as 0/0. With the gVCF, a site is missing only when it was not covered.

Without a gVCF (`Input` is `vcf`, a sample called by an older version), the sum misses the weight of every reference-allele effect allele you carry on two copies, and few scores reach pgsc_calc's 75% match rate. The step then prints:

```
NOTE: hom-ref sites are absent from this VCF, so the score is biased; not comparable to published distributions
```

Call the sample again with step 3 to get the gVCF.

### What these scores are NOT

- A raw sum is NOT a percentile. Without the panel there is no percentile at all.
- A percentile is NOT a probability. Being at the 90th percentile does not mean a 90% chance of the condition.
- They are NOT comparable across conditions.
- The percentile compares you with the reference samples most similar to you, not with your own family or community. For a group the panel represents poorly, the percentile is less reliable; the ancestry table says when the population match is low-confidence.
- They are NOT stable across scoring file versions. A new harmonised file is a new baseline.

## Limitations

- Most scores were developed in European-ancestry populations and predict less well for other ancestries, even with an ancestry-adjusted percentile.
- A score captures only the genetic component. Lifestyle, environment, and family history are often more predictive.
- Sex-specific conditions (breast cancer, prostate cancer) should be interpreted accordingly.
- Chip data matches far fewer score variants than a genome; most scores then fall under pgsc_calc's 75% match rate and get no sum (see [Using chip data](chip-data-guide.md)).

## Notes

- Scoring files are cached in `${GENOME_DIR}/prs_scores/`. Delete a file to download it again. A cached file without a `#HmPOS_build=GRCh38` header is deleted and downloaded again.
- pgsc_calc starts with a fresh work folder each run: it keeps converted genotypes between runs, and an old run's would be scored instead of the current VCF.

## Maintenance

- Recheck the PGS Catalog and pgsc_calc releases at least quarterly. A pgsc_calc bump means reading its `conf/modules.config` at the new release and moving the `PGSC_*_IMAGE` lines of `versions.env` with it, plus `PGSC_CALC_SHA256` and `PGSC_CALC_NF_SCHEMA`.
- A scoring file update is a **result-changing event**: rerun step 25 and treat the output as a new baseline.
- Keep the `PGS ID`, the harmonised file version, the pgsc_calc release and the pipeline commit together so score changes stay auditable.

## Links

- [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) and its [documentation](https://pgsc-calc.readthedocs.io/)
- [PGS Catalog](https://www.pgscatalog.org/)
- [PGS Catalog scoring file format](https://www.pgscatalog.org/downloads/#scoring_files)
- [Khera et al. 2018 (multi-trait PRS)](https://doi.org/10.1038/s41588-018-0183-z)
