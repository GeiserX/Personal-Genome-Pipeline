# Step 44: Ancient Populations (AADR PCA)

## What This Does

Places your genome on a principal component analysis (PCA) of present-day West Eurasians, the same axes onto which about 11,000 ancient West Eurasians of the Allen Ancient DNA Resource (AADR) are projected, and lists the ancient groups that sit nearest to you, with their dates and places. It also draws PC1 against PC2 with every ancient person coloured by age and you marked. Opt-in.

Closeness on a PCA is genetic similarity, not descent. A genome near an ancient group resembles it on these axes; that does not mean you descend from those people, and two very different histories can land on the same spot. The step does not fit qpAdm models of ancestry proportions: that needs chosen source and outgroup populations and f-statistics, and it is out of scope here.

## Why

Most people who look at their ancestry want to know where they fall among ancient genomes, and the AADR is the public, curated collection of those genomes (CC0). The standard way to compare a present-day genome with ancient ones is the one the ancient DNA papers use: compute the axes from present-day people only, then project the ancient people onto them by least squares (smartpca's `lsqproject`), because ancient genomes have most sites missing and would distort axes they helped build.

## Tool

- **smartpca** and **convertf** (EIGENSOFT 9.0.0): the PCA, and the conversion of the present-day people to EIGENSTRAT text. 9.0.0 is the first release that reads the AADR's TGENO genotype format.
- **bcftools** (`convert --gvcf2vcf`, `view`, `query`) reads your genotypes at the panel's sites from the gVCF.
- **`bin/ancient_report.py`** (Python standard library): chooses the panel's people, lifts the AADR sites to GRCh38, adds your column to the present-day genotypes, aligns the axes and writes the report and the plot.

## Docker Image

- `EIGENSOFT_IMAGE`, with `BCFTOOLS_IMAGE` and `PYTHON_IMAGE`, pinned in `versions.env`; [Image versions](versions.md) lists the tags.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz` with its `.tbi`: the gVCF of [step 3](03-variant-calling.md). A site where you match the reference reads as a real genotype and a site without reads stays missing; a variant-only VCF cannot tell those apart. The gVCF must hold one sample.
- The panel `setup.sh --ancient-panel` builds in `reference/ancient/aadr_<release>/` (`ANCIENT_PANEL_DIR` picks another folder inside `GENOME_DIR`, `ANCIENT_PANEL_DIR=none` skips the step). Without it the step prints one line saying so and exits 0.

In the Nextflow pipeline (`ancient` in `--tools`, with `--ancient_panel <folder>`) the processes `ANCIENT_SITES`, `ANCIENT_GENO`, `ANCIENT_PCA` and `ANCIENT_REPORT` run the same commands for every sample with a gVCF; a sample without one is skipped with a warning. `run-all.sh` lists the step only when `TOOLS` names `ancient`, and passes the panel when its `DONE` file exists.

## The Panel

`scripts/setup.sh --ancient-panel <genome_dir>` builds it once:

1. Downloads the AADR release `AADR_VERSION` (v66.p1, 1240K) from Harvard Dataverse: the `.anno`, `.ind`, `.snp` and the 7.1 GB `.geno`, each checked against the md5 Dataverse publishes (`AADR_*_MD5` in `versions.env`). UCSC's hg19ToHg38 chain is checked against `HG19_TO_HG38_CHAIN_SHA256`.
2. Chooses the people. The axes come from 41 present-day groups of Europe, the Caucasus and the Near East (shotgun diploid genomes of the 1000 Genomes Project, HGDP and SGDP, about 820 people; North African groups and groups with recent East Eurasian ancestry are left out). The projected people are the ancient ones the AADR's own assessment passes, dated, with at least 100,000 of the 1240K sites covered, between latitudes 25 and 75 and longitudes -25 and 75, one entry per individual (the one with most sites).
3. Lifts the 1240K sites from GRCh37 to GRCh38 and keeps the autosomal sites one of whose two alleles is the GRCh38 reference base; a site the chain places twice, or on another chromosome, is dropped.
4. Runs smartpca on the whole release with the present-day groups as the population list and every ancient person projected (`panel.evec`), then runs it again on the present-day people alone, as step 44 does, and stops if any of the first four axes (the ones the distances use) of the two runs correlate below 0.999 up to sign.

`ANCIENT_SOURCE_DIR=<dir>` copies the four AADR files from a folder instead of downloading them (not checked), and `ANCIENT_CHAIN=<file inside GENOME_DIR>` uses another chain file. The download is deleted once the panel is built.

## How Step 44 Works

1. bcftools reads your genotype at each panel site from the gVCF.
2. Your genotypes become one more column of the present-day panel's genotypes, counted as copies of the AADR's first allele; a site with no genotype, or with an allele the AADR does not list, is missing.
3. smartpca builds the axes from the same present-day people and projects you onto them (`lsqproject`). Each axis is then turned to the sign of the panel's run, through the present-day people both runs share, so your position and the ancient people's positions are on the same axes. The step stops if an axis disagrees (correlation below 0.999): the panel and the run do not match.
4. Each ancient group (the AADR's Group ID) gets the mean position of its people. The distance from you to that mean uses the first four PCs, each weighted by the square root of its eigenvalue, so an axis that explains more variation counts more.

## Command

```bash
./scripts/setup.sh --ancient-panel /path/to/genome_dir   # once
./scripts/44-ancient.sh your_name
```

## Output

| File | Contents |
|---|---|
| `ancient/${SAMPLE}_ancient.tsv` | the 20 nearest ancient groups: `rank`, `group`, `n` (people in the group), `distance`, `pc1`, `pc2` (the group's mean), `date_bp` (the median date, years before 1950), `dates` (the group's span, in BCE/CE), `country`, `locality` (the most common in the group) |
| `ancient/${SAMPLE}_ancient_info.tsv` | `release`, `sites_used` and `panel_sites`, the numbers of present-day and ancient people and of ancient groups, the PCs in the distance, your `pc1` and `pc2` |
| `ancient/${SAMPLE}_ancient.txt` | the same in words, with the caveat |
| `ancient/${SAMPLE}_ancient_pca.svg` | PC1 against PC2 on a dark background: present-day people in grey, ancient people coloured by age, the five nearest groups named, you as a ring |

Dates follow the AADR: years before present count from 1950, so 4250 BP is 2300 BCE, and 1950 BP is 1 BCE.

Both reports have an Ancient Populations section: the ten nearest groups, the caveat, and in the HTML report the plot.

## Runtime

- Setup: the 7.1 GB download, then convertf and two smartpca runs over the whole release. smartpca holds the packed genotypes of every person in memory, about 7 GB, so the setup container gets 12 GB.
- Step 44: minutes. bcftools reads the gVCF at the panel's sites through its index, and smartpca works on about 820 people plus you; it gets 8 GB.

## Notes

- A group with one or two people has a noisy mean; the `n` column says how many there are.
- The ranking compares projected positions only. smartpca's `shrinkmode`, which corrects the shrinking of projected samples towards the centre, is not used: it would need one double per person and site of the whole release in memory, and the people it corrects for are the ones that build the axes, not the projected ones.
- The axes describe variation within West Eurasia. A genome with most of its ancestry from elsewhere falls at the edge of the plot, and its nearest groups say little.
- Ancient genomes are pseudo-haploid at most sites and differ in coverage; `lsqproject` places them without the missing sites, but a low-coverage person still lands less precisely than a present-day genome.
- The AADR is distributed under CC0. Cite it as Mallick S, et al. The Allen Ancient DNA Resource (AADR): a curated compendium of ancient human genomes. Sci Data 11, 182 (2024), https://doi.org/10.1038/s41597-024-03031-7, and the release as https://doi.org/10.7910/DVN/FFIDCW.
- Further reading: Patterson N, Price AL, Reich D. Population structure and eigenanalysis. PLoS Genet 2, e190 (2006), https://doi.org/10.1371/journal.pgen.0020190 (smartpca); Lazaridis I, et al. Ancient human genomes suggest three ancestral populations for present-day Europeans. Nature 513, 409-413 (2014), https://doi.org/10.1038/nature13673 (the projection of ancient genomes onto present-day West Eurasian axes); EIGENSOFT, https://github.com/DReichLab/EIG. Links checked on 2026-10-25.
