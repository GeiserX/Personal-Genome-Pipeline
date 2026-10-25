# Step 37: Y-Chromosome Haplogroup

## What This Does

Assigns the Y-chromosome haplogroup of a male sample: the branch of the paternal line's family tree its Y variants place it on. Opt-in.

## Why

The mitochondrial haplogroup (step 12) traces the maternal line only. The Y chromosome passes from father to son almost unchanged, so its haplogroup traces the paternal line. Like step 12, it is ancestry, not health: no Y haplogroup is a diagnosis.

## Tool

- **Yleaf** 3.2.1 (Ralf et al., Mol Biol Evol 2018; GPL-3.0), Erasmus MC. It reads the BAM's pileup at its Y markers and predicts the haplogroup they support.

Upstream Yleaf is at 4.x; Bioconda has 3.2.1 only, and that build is the image here.

## Docker Image

- `YLEAF_IMAGE`, and `SAMTOOLS_IMAGE` for the pileup

The Yleaf image has no samtools, which Yleaf calls for a BAM. So the step runs in three parts (`bin/yleaf_run.py`): Yleaf's marker positions from its image, `samtools idxstats` and `samtools mpileup -l <positions> -AQ20q1` (Yleaf's own flags at its default quality 20) in `SAMTOOLS_IMAGE`, then Yleaf on that pileup.

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps)
- `${SAMPLE}/indexcov/indexcov-indexcov.ped` from step 16: the sex indexcov infers from the reads. The step runs only when it is male, and asks for no sex of its own. A female or undetermined sample is skipped with one line and exit 0.

In the Nextflow pipeline (`y_haplogroup` in `--tools`) the process runs for the rows whose samplesheet sex is male. `INDEXCOV` has checked that sex against the reads before any BAM step starts, and a mismatch stops the run (unless `--sex_check warn`).

## Command

```bash
./scripts/setup.sh --yleaf-data /path/to/genome_dir   # once: Yleaf's marker tables and tree (~15 MB)
./scripts/37-y-haplogroup.sh your_name
```

The Yleaf image installs Yleaf's code but not its data folder (the GRCh38 marker positions and the haplogroup tree). `setup.sh --yleaf-data` takes that folder from GitHub's archive of the same release, checked against `YLEAF_DATA_SHA256` in `versions.env`, into `reference/yleaf-<version>/data`; the pipeline takes it as `--yleaf_data`.

With `run-all.sh`: `TOOLS=...,y_haplogroup` (it is opt-in, so a default run lists it as skipped; without the data it is skipped with the reason).

Yleaf downloads the whole hg38 FASTA on its first run unless its config file names one, and the image's config is read-only. The launcher points Yleaf at the pipeline's reference before it starts (a BAM never needs the sequence itself), so nothing is downloaded and every container runs without a network.

## Output

| File | Contents |
|---|---|
| `y_haplogroup/${SAMPLE}_y_haplogroup.txt` | Yleaf's prediction: `Hg` (the haplogroup), `Hg_marker`, `Total_reads`, `Valid_markers` (Y markers with haplogroup information and enough reads), `QC-score` |
| `y_haplogroup/yleaf/` | Yleaf's working files: the markers it read with their alleles and its log |
| `y_haplogroup/positions.txt` | the marker positions the pileup was made at |
| `y_haplogroup/${SAMPLE}_chrY.bam` and `.bai` | the BAM's chrY reads only, the file the services below read; the Nextflow pipeline publishes the same file there (`Y_PILEUP`) |

The chrY BAM is written on every male run. Expect it to be about 1 to 2 percent of the analysis BAM, so a few hundred megabytes to a gigabyte on a 30x genome. That is an estimate from chrY's share of the reads; we have not measured it on a full genome. If you do not need it, delete it and its `.bai`.

`Hg` is `NA` when too few markers had reads for a call; the step then prints `Y haplogroup: insufficient markers` with the marker count, and both reports say "insufficient markers".

With a call, both reports add the links of the next section, built from the branch name (for `R-M269*(xR-L23)`, Yleaf's mark for "R-M269 but none of these sub-branches", the links use `R-M269`).

## Explore your branch further

A haplogroup is a branch of the paternal family tree. Other people's tests tell you more about that branch: how old it is, where its men live today, which ancient remains carry it. The routes below are free first. Each was checked on its own pages on 2026-10-25; prices and policies change, so read the service's page before you send anything.

**This pipeline never uploads anything.** An upload shares the file with that service, under its own terms. A BAM of chrY holds the reads of your paternal line only, but it is still your DNA, and through the paternal line it is your father's, brothers' and sons' too. It also keeps the BAM's header (the sample name, the read groups and the command lines of the tools that made it) and, for each read, where its mate aligned. Look at the header with `samtools view -H` before you send it.

| Route | Cost | What you send | What you get | Checked 2026-10-25 |
|---|---|---|---|---|
| FamilyTreeDNA Discover, `https://discover.familytreedna.com/y-dna/<branch>` | free, no account | no file: the branch name goes in the link | the branch's age (when it formed, and its most recent common ancestor), the countries its testers name, ancient DNA samples on it, a migration map, its sub-branches | the page loads its data from `https://discover.familytreedna.com/resources/y-dna/<branch>.json`, which answered for `R-M269` with no login; its search (`/HaplogroupSearch/<name>`) found `R-M269` exactly |
| YFull YTree, `https://www.yfull.com/tree/<branch>/` | free to browse | no file: the branch name goes in the link | the branch on YFull's tree, with its age estimates and the samples placed under it | `/tree/R-M269/` answers 200, a branch YFull does not have answers 404 |
| Decoding Us Navigator, `https://decoding-us.com/download` | free | nothing: the app runs on your computer | haplogroup placement and coverage metrics from your BAM or CRAM; you choose later whether to publish anonymised summaries | the site says it accepts no file uploads and that raw data never leaves your computer; the Navigator is an alpha (v0.1.0-beta.1, 2026-09-03). Its trees are built on T2T-CHM13; we have not run it on a GRCh38 BAM |
| YFull Y analysis, `https://www.yfull.com/order/` | €44 plus VAT, paid after the analysis | a link to a shared BAM or CRAM (YFull does not accept FASTQ; a VCF only for Big Y) | your terminal branch, your novel SNPs, STRs, age estimates, and comparison with YFull's samples; the mtDNA analysis is then free in the same account | YFull's pages do not say whether a BAM with chrY only is enough; we have not tried it |

The two volunteer projects older guides name:

- **James Kane's YDNA Warehouse**: `ydna-warehouse.org` now redirects to Decoding Us (the row above).
- **Alex Williamson's The Big Tree**: `thebigtree.org` had no DNS record on 2026-10-25, so we could not check whether it takes uploads. Treat it as unavailable.

Names differ between trees. Yleaf names branches the YFull way; FamilyTreeDNA keeps its own tree. If Discover shows no page for your branch, type the branch name, or the SNP after the dash, into Discover's search: it offers the nearest branch it has (a search for YFull's `A0-T` returns FamilyTreeDNA's `A-L1098`).

**Private SNPs**, the Y mutations you carry that are on no public tree, are what would place you below your branch, next to the men closest to you. This pipeline does not list them yet: that needs the callable regions of chrY and a public SNP list beyond the markers Yleaf ships. YFull's analysis reports them as novel SNPs.

## Runtime

Not measured here; the end test does not run this step. An estimate is a few minutes: Yleaf reads the pileup at about 129,000 Y marker positions.

## Notes

- Yleaf 3 places the sample on the YFull tree (v10.01) and names the haplogroup by clade and defining marker, for example `R-M269`. Its `QC-score` (0 to 1) is how consistently the markers on the path to that branch agree; the default acceptance threshold is 0.95.
- Short reads cover only the unique parts of the Y; that is enough for a haplogroup at 30x.
- On the e2e fixture (300 kb of chrY) the step runs on a copy of HG002 with a male indexcov row, because indexcov reads the fixture's slices as female.
