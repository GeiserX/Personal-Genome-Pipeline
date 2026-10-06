# Step 26: Ancestry (projection onto a reference panel)

## What This Does

Places your genome among the samples of a reference panel whose populations are known, and names the panel population your genetic ancestry is most similar to. It uses pgsc_calc's ancestry projection: the panel's principal components are computed from the panel's own samples, and your genotypes are projected onto them, so one sample is enough. The same pgsc_calc run gives step 25 its percentiles.

Without the panel the step prints one line and exits 0:

```
Step 26 skipped: no ancestry reference panel at <genome_dir>/reference/pgsc_calc/pgsc_1000G_v1.tar.zst; install it with scripts/setup.sh --ancestry-panel <genome_dir>
```

## Why

1. **PRS interpretation**: a polygenic score (step 25) only means something next to people of similar genetic ancestry. The population found here is the group step 25 compares your score with.
2. **Context for other results**: population frequencies and some risk estimates differ between ancestries.

A principal component analysis of one genome alone cannot work, since the axes come from the variation between many people. Projection onto axes that a reference panel defines is the method that works on one sample.

## Tool

- **pgsc_calc** (`PGSC_CALC_VERSION` in `versions.env`), run by step 25: FRAPOSA's online augmentation, decomposition and Procrustes projection (`--projection_method oadp`), then a random forest on the first principal components to assign the most similar population.
- **Reference panel**: pgsc_calc's 1000 Genomes database `pgsc_1000G_v1` (`PGSC_PANEL` in `versions.env`), samples of the five 1000 Genomes super-populations (AFR, AMR, EAS, EUR, SAS), published by the PGS Catalog.

## Docker Images

The images of step 25, all pinned in `versions.env`: `PGSC_UTILS_IMAGE`, `PLINK2_IMAGE`, `PGSC_FRAPOSA_IMAGE`, `PGSC_ZSTD_IMAGE`, `PGSC_PYYAML_IMAGE`, `PGSC_REPORT_IMAGE`, with `PYTHON_IMAGE` and `BCFTOOLS_IMAGE`. [Image versions](versions.md) lists the tags.

## Input

- VCF from DeepVariant (step 3), and its gVCF beside it: the panel's SNVs are genotyped from the gVCF, so the sites where you match the reference count in the projection. Without a gVCF only your variant sites are projected, which weakens it.
- The panel and its site list, installed once (7.4 GB, and about 24 GB of disk while pgsc_calc runs): `./scripts/setup.sh --ancestry-panel <genome_dir>`. It is opt-in for that reason; see the measured numbers in [step 25](25-prs.md#the-reference-panel-on-a-github-hosted-runner).
- Java 17+ and Nextflow, as for step 25.

## Command

```bash
./scripts/setup.sh --ancestry-panel /path/to/genome_dir   # once
./scripts/26-ancestry.sh your_name
```

`run-all.sh` runs it after the pipeline when asked:

```bash
ANCESTRY=true ./scripts/run-all.sh your_name <male|female>
```

`ANCESTRY_PANEL=/path/to/panel.tar.zst` points at another pgsc_calc panel (its `_GRCh38_sites.tsv` must be beside it); `ANCESTRY_PANEL=none` skips the step.

## What the Script Does Internally

1. Without the panel: prints the one line above and exits 0.
2. Runs step 25, which uses the panel when it is installed: the score and panel positions are genotyped from the gVCF, and pgsc_calc runs with `--run_ancestry`. pgsc_calc intersects your genotypes with the panel, keeps unrelated panel samples and common, LD-thinned SNVs, computes the panel's principal components, projects you onto them, and assigns the population with a random forest trained on the panel's labels.
3. Prints the population and checks that the ancestry table was written.

In the Nextflow pipeline the same happens inside `prs` when `--ancestry_ref` is set; `ancestry` in `--tools` needs `prs` and `--ancestry_ref`.

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_ancestry.tsv` | `key` and `value` rows: `sample`, `reference_panel`, `population` (the most similar panel population), `population_low_confidence`, `probability_<POP>` for each panel population, and `PC1` to `PC10` |

Written to `${GENOME_DIR}/${SAMPLE}/ancestry/`. pgsc_calc's own files, including the principal components of every panel sample, are in `${GENOME_DIR}/${SAMPLE}/prs/pgsc_calc/results/sample/score/` (`sample_popsimilarity.txt.gz`).

## Runtime

The panel's extraction, QC and PCA come on top of step 25's run; see the measured numbers in [step 25](25-prs.md#the-reference-panel-on-a-github-hosted-runner).

## Interpreting Results

- **population** is the reference group whose genetic ancestry is most similar to yours. It is a statement about similarity to five continental groups of the 1000 Genomes Project, not about identity, nationality or ethnicity.
- **probability_<POP>** are the random forest's probabilities. A low top probability (`population_low_confidence` is `True`) means you sit between groups or far from all of them; the step 25 percentile is then less reliable.
- **PC1 to PC10** are your coordinates on the panel's principal components. They can be plotted against the panel samples in `sample_popsimilarity.txt.gz`.

## Limitations

- Five continental groups only. Mixed ancestry is shown as the single most similar group, with its probabilities; this step does not estimate admixture fractions or local ancestry.
- Groups the 1000 Genomes Project does not sample well are placed less reliably.
- Fine-grained ancestry (for example Spanish versus Italian) needs other panels and tools.

## Notes

- The panel is one file kept as downloaded; pgsc_calc unpacks the GRCh38 part into its work folder on each run. The site list beside it (`pgsc_1000G_v1_GRCh38_sites.tsv`) is made by `setup.sh` with plink2 from the panel's own genotypes: the biallelic autosomal SNVs with a panel frequency of 5% or more, the threshold pgsc_calc's projection applies (`maf_ref`).
- pgsc_calc's synthetic HAPNEST panel (`GRCh38_HAPNEST_reference`, 268 MB) is what the pull-request tests project onto; `ANCESTRY_PANEL_NAME=GRCh38_HAPNEST_reference ./scripts/setup.sh --ancestry-panel <genome_dir>` installs it, but its populations are simulated and say nothing about you.

## Links

- [pgsc_calc ancestry documentation](https://pgsc-calc.readthedocs.io/en/latest/explanation/geneticancestry.html)
- [1000 Genomes Project](https://www.internationalgenome.org/)
- [FRAPOSA](https://github.com/daviddaiweizhang/fraposa) and the [fork pgsc_calc runs](https://github.com/PGScatalog/fraposa_pgsc)
