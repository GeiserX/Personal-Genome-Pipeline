# Step 26b: Fine-Grained Ancestry (admixture fractions and a chromosome painting)

## What This Does

Estimates which reference groups your genome is most similar to, as fractions with a bootstrap range, at two levels: 6 continental groups and 18 regional ones. It also paints each chromosome: each segment of each of your two copies gets the continental group it is most similar to. Opt-in.

Without the panel the step prints one line and exits 0:

```
Step 26b skipped: no fine-grained ancestry panel at <genome_dir>/reference/ancestry_fine/hgdp1kgp_v2; install it with scripts/setup.sh --ancestry-fine-panel <genome_dir>
```

## Why

[Step 26](26-ancestry.md) names the most similar of five continental groups of the 1000 Genomes Project. It gives no fractions and cannot tell a mix of two groups from a person between them. The methods behind the fractions and the painting of commercial ancestry tests are open; what differs is the reference panel. This step uses an open panel of 80 populations and open tools, so you can see the result and also where it stops being supported.

## Tools

Each choice is the tool with a maintained, pinned BioContainers image that fits a GitHub runner on a reduced panel.

- **ADMIXTURE** 1.3.0 (Alexander et al., Genome Research 2009) for the fractions. The panel's allele frequencies of each group are fitted once, in supervised mode (every reference person's group is given), and your genome is projected on them (`-P`). SCOPE, the faster alternative, has no BioContainers image.
- **Beagle** 5.5 (Browning et al., AJHG 2021; GPL-3.0) phases your genotypes against the panel, since the painting needs the two copies apart.
- **RFMix** 2.03 (Maples et al., AJHG 2013) for the painting. FLARE, the other candidate, has no BioContainers image.
- **bin/ancestry_fine_report.py** writes the outputs and the bootstrap.

**Licences: ADMIXTURE and RFMix are free for academic and non-commercial use only** (RFMix prints "licensed for non-commercial academic research use only" when it starts). Check them before any other use.

ADMIXTURE's own bootstrap (`-B`) crashes with a segmentation fault when it is combined with projection (`-P`) in this build. The spread comes instead from a block bootstrap run by the step: each of 100 replicates draws 100 blocks of neighbouring SNVs with replacement and projects them again.

## The Reference Panel

gnomAD's HGDP + 1000 Genomes callset (Koenig et al., Genome Research 2024), phased with SHAPEIT5: 4,091 people of 80 populations sequenced to about 30x on GRCh38. It is part of gnomAD's genome data, which gnomAD releases under the [Creative Commons Zero dedication](https://gnomad.broadinstitute.org/policies): open, no access agreement. `setup.sh --ancestry-fine-panel` reads it from gnomAD's public bucket (`resources/hgdp_1kg/phased_haplotypes_v2`), a region at a time, and keeps:

- the biallelic SNVs with a minor allele frequency of 5% or more, one per 2 kb (`ANCESTRY_FINE_THIN`), for everyone in the panel: phasing reference and painting reference;
- the people's groups (`samples.tsv`), from gnomAD's metadata of the release;
- after LD pruning (r2 above 0.1 within 500 kb removed, among the unrelated reference people), the SNVs ADMIXTURE uses, and each group's allele frequencies at both levels;
- Beagle's GRCh38 genetic map (HapMap, lifted from GRCh37).

The bucket reports a crc32c checksum per file. `versions.env` pins the sha256 of those checksums and sizes (`ANCESTRY_FINE_PANEL_SHA256`), and of the metadata and the map, so a file gnomAD replaces stops the install rather than changing the panel silently.

### The groups

A group name is ours: a set of reference populations we put together and named, not a nationality. Admixed populations (between continental groups) are in no group and are never a reference: ACB, ASW, CLM, MXL, PEL, PUR, Hazara and Uygur. Related people are no reference either.

| Continental group | Regional groups (reference populations) |
|---|---|
| Africa | West Africa (YRI, ESN, GWD, MSL, Yoruba, Mandenka); East Africa (LWK, BantuKenya); Central African rainforest (Biaka, Mbuti); Southern Africa (San, BantuSouthAfrica) |
| West Eurasia (Europe, Middle East, North Africa) | Northwestern Europe (CEU, GBR, Orcadian, French); Northeastern Europe (FIN, Russian); Iberia (IBS, Basque); Italy (TSI, Tuscan, BergamoItalian); Sardinia (Sardinian); Caucasus (Adygei); Levant (Druze, Palestinian); Bedouin (Bedouin); North Africa (Mozabite) |
| Central and South Asia | Central and South Asia (GIH, PJL, BEB, STU, ITU, Balochi, Brahui, Makrani, Pathan, Sindhi, Burusho, Kalash) |
| East Asia | East Asia (CHB, CHS, CDX, JPT, KHV, Han, NorthernHan, Japanese, Dai, Daur, Hezhen, Lahu, Miao, Mongolian, Naxi, Oroqen, She, Tu, Tujia, Xibo, Yi, Cambodian); Siberia (Yakut) |
| Americas (Native American) | Americas (Karitiana, Maya, Pima, Surui, Colombian) |
| Oceania | Oceania (Bougainville, PapuanHighlands, PapuanSepik) |

Europe, the Middle East and North Africa form one continental group: at that level their populations are not well apart, and a painting with them split flips between them along every European or Middle Eastern genome. The regional level splits them.

`ANCESTRY_FINE_SAMPLE_MAP=<file>` installs a `samples.tsv` of your own (the same six columns) and fits the frequencies again; the panel's genotypes are kept.

### Not in this panel

The Allen Ancient DNA Resource (AADR, release v66.p1, CC0 1.0 on the [Harvard Dataverse](https://doi.org/10.7910/DVN/FFIDCW)) has present-day people of many more regions, Greek and regional Spanish among them. It is not used here: its present-day data are array genotypes on GRCh37 in the packed EIGENSTRAT format, so they need a liftover and a converter no pinned image here has, and they are unphased, so they could not feed the painting. qpAdm models of ancient sources are not part of this step either.

## Docker Images

`ADMIXTURE_IMAGE`, `BEAGLE_IMAGE`, `RFMIX_IMAGE`, with `BCFTOOLS_IMAGE`, `PYTHON_IMAGE` and, for the install, `BCFTOOLS_HTTPS_IMAGE` (the panel is read over https; `BCFTOOLS_IMAGE` has no CA certificates) and `PLINK2_IMAGE`. All are pinned in `versions.env`; [Image versions](versions.md) lists the tags. `ADMIXTURE_IMAGE`, `BEAGLE_IMAGE`, `RFMIX_IMAGE` and `BCFTOOLS_HTTPS_IMAGE` are `# optional`: `setup.sh --ancestry-fine-panel` pulls them, a plain `setup.sh` does not.

## Input

- VCF from DeepVariant (step 3) and its gVCF beside it: the panel's SNVs are genotyped from the gVCF, so the sites where you match the reference count.
- Without a gVCF, a VCF that holds reference calls (0/0), such as the output of `chip-to-vcf.sh`. A variant-only VCF stops the step: from the sites where you differ from the reference alone, the fractions would lean toward the groups where those alleles are common. The rule: at least one in five of the genotyped panel SNVs must be 0/0.
- The panel, installed once: `./scripts/setup.sh --ancestry-fine-panel <genome_dir>`.

## Command

```bash
./scripts/setup.sh --ancestry-fine-panel /path/to/genome_dir   # once
./scripts/26b-ancestry-fine.sh your_name
```

`run-all.sh` runs it in the pipeline whenever the panel is installed (`ancestry_fine` in `--tools`, `--ancestry_fine_panel` set); `ANCESTRY_FINE_DIR=none` skips it. In the Nextflow pipeline the step is `ancestry_fine` in `--tools` with `--ancestry_fine_panel <genome_dir>/reference/ancestry_fine/hgdp1kgp_v2`. Settings: `THREADS` (default 4), `ANCESTRY_FINE_BOOTSTRAP` (replicates, default 100; 0 for none), `ANCESTRY_FINE_DIR` (another panel folder inside the genome directory). `ANCESTRY_FINE_BOOTSTRAP` applies to `26b-ancestry-fine.sh` only: the pipeline, and so `run-all.sh`, always runs 100 replicates.

## What the Script Does Internally

1. Genotypes the panel's SNVs from the gVCF (`bcftools convert --gvcf2vcf`, as step 25 does) and keeps a record only with the panel's two alleles.
2. Writes a one-person ADMIXTURE input on the SNVs ADMIXTURE uses that you have a genotype for (ADMIXTURE stops on a SNV with no genotype, so a missing one is left out), projects it on each level's frequencies, then runs the bootstrap replicates.
3. Chooses the chromosomes to paint: 100 genotyped SNVs or more spanning 5 cM or more. RFMix stops with a segmentation fault on a span of a few cM; on a whole genome every autosome spans more than 40 cM.
4. Phases those chromosomes with Beagle against the panel (`impute=false`, `seed=1`).
5. Paints them with RFMix against 50 unrelated reference people per continental group at most, taken in turn from each population (`ANCESTRY_FINE_PAINT_PER_GROUP` at install), `--random-seed=1`.
6. Writes the outputs below.

The Nextflow module (`modules/local/ancestry_fine`) runs the same parts as six processes, one per image, and the e2e test checks that it writes the fractions of the bash step.

## Output

| File | Contents |
|---|---|
| `ancestry_fine/${SAMPLE}_ancestry_fine.tsv` | one row per group: `level` (`continent`, `region`, or `painting` for the share of the painted genome), `code`, `group`, `fraction`, `se` (bootstrap standard error), `low` and `high` (2.5 and 97.5 percentiles of the bootstrap) |
| `ancestry_fine/${SAMPLE}_painting.bed` | the painting: chromosome, start, end, haplotype (1 or 2), continental group code |
| `ancestry_fine/${SAMPLE}_painting.png` | chromosomes 1 to 22 to scale, two bars each, coloured by continental group; grey is not painted |
| `ancestry_fine/${SAMPLE}_ancestry_fine.txt` | the same in plain text, with what it means |
| `ancestry_fine/work/` | the genotypes, the phased haplotypes, ADMIXTURE's and RFMix's own files |

The HTML and text reports have a section with the fractions, the painted shares and the picture.

## Runtime

Measured in CI only, on the e2e job's reduced panel: the fixture's slices plus chr22:20-35 Mb, 4 CPUs.

| Part | Time | Notes |
|---|---|---|
| `setup.sh --ancestry-fine-panel` on that panel | 1 min 43 s | 154 MB on disk, 7,789 SNVs of 4,089 people, 1,003 of them for ADMIXTURE |
| The step on the fixture's slices from a gVCF (no chromosome long enough to paint) | 4 s | 2,136 panel SNVs genotyped |
| The step with the painting of chr22:20-35 Mb | 1 min 19 s | 100 bootstrap replicates, 236 painting references |
| The same through the Nextflow pipeline, with the HTML report | 1 min 38 s | six processes |

A whole genome has not been run here. Expected, scaled from those numbers: the install reads about 19 GB from gnomAD's bucket and decodes every person's genotypes, a few hours on a fast connection; the step itself is dominated by RFMix, which grows with the number of SNVs and painting references.

## Interpreting Results

- **The fractions** say which reference groups your genome is most similar to, and in what proportion. They are not percentages of your ancestors from a country.
- **The bootstrap range** is how much a fraction moves when the SNVs are resampled. A regional fraction with a wide range is not supported by the data; read it as "somewhere in this range", never as one exact percentage.
- **The two levels** are two answers to one question. When the regional fractions inside West Eurasia add up to the continental one, the regional split is how the step divides it, not a second measurement.
- **The painting** shows where along each chromosome each continental group is most similar. A short segment of another group can be a real recent ancestor or noise; long segments are the stronger evidence.

## Limitations

- A label such as Iberia is a cluster we name after the reference populations in it, not a nationality or an ethnicity.
- Continental-level similarity is well supported. Country-level similarity is supported only where the panel has that country and its neighbours are genetically apart. Below that the step says little.
- Splits within Iberia are weak for every method; the panel has IBS (Spanish) and Basque only.
- Greek and Southern Italian ancestry overlap, and the panel has no Greek population at all: a Greek genome is described with Italy, the Levant and the Caucasus.
- Groups the panel samples poorly (Native American, Oceanian, North African) are estimated less well.
- A population the panel lacks, Ashkenazi Jewish for one, gets no group of its own: its genomes are described with the nearest groups.
- The painting uses continental groups only; local ancestry within Europe is not reliable with these tools.

## Notes

- The e2e case (`tests/e2e/ancestry-fine-grained-1.sh`) checks HG002 from its gVCF, HG002's GIAB genotypes on chr22:20-35 Mb, and a planted person whose first half comes from a Yoruba and second half from a Han Chinese person (both left out of the panel): the planted person must come out near half and half with a painting that switches at the junction. A red control in the same case moves every reference person's label to the next group and checks that the HG002 checks then fail.

## Links

- [gnomAD downloads](https://gnomad.broadinstitute.org/downloads), where the HGDP + 1000 Genomes resource is listed, and its [data policy](https://gnomad.broadinstitute.org/policies)
- [ADMIXTURE](https://dalexander.github.io/admixture/)
- [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
- [RFMix](https://github.com/slowkoni/rfmix)
