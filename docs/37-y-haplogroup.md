# Step 37: Y-Chromosome Haplogroup

## What This Does

Assigns the Y-chromosome haplogroup of a male sample: the branch of the paternal line's family tree its Y variants place it on. Opt-in.

## Why

The mitochondrial haplogroup (step 12) traces the maternal line only. The Y chromosome passes from father to son almost unchanged, so its haplogroup traces the paternal line. Like step 12, it is ancestry, not health: no Y haplogroup is a diagnosis.

## Tool

- **Yleaf** 3.2.1 (Ralf et al., Mol Biol Evol 2018; GPL-3.0), Erasmus MC. It reads the BAM's pileup at its Y markers and predicts the haplogroup they support.

Upstream Yleaf is at 4.x; Bioconda has 3.2.1 only, and that build is the image here.

- **The current YFull YTree** (release 14.06.0 of 22 September 2026, CC BY 4.0, [github.com/YFullTeam/YTree](https://github.com/YFullTeam/YTree)), walked by `bin/y_tree_place.py`. Yleaf 3.2.1 ships YFull's tree v10.01, so the branch it reports can sit several levels above the one the current tree resolves. The SNPs on YFull's tree have names only; their GRCh38 positions and bases come from two tables of Yleaf 4.1.6 (GPL-3.0): `ftdna_positions_hg38.txt` (made for FamilyTreeDNA's tree, the widest of its tables) first, then `new_positions.txt` (the one Yleaf 4 uses for its YFull v14 tree). SNP names are shared between the trees, so a name is looked up as it stands. Together they place 85% of the tree's 459,468 SNPs (76% from the first table alone), counted on 2026-10-25.

The step reports both results, each with the tree it comes from.

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
./scripts/setup.sh --yleaf-data /path/to/genome_dir   # once: Yleaf's marker tables and tree, and the current tree (~75 MB)
./scripts/37-y-haplogroup.sh your_name
```

The Yleaf image installs Yleaf's code but not its data folder (the GRCh38 marker positions and the haplogroup tree). `setup.sh --yleaf-data` takes that folder from GitHub's archive of the same release, checked against `YLEAF_DATA_SHA256` in `versions.env`, into `reference/yleaf-<version>/data`; the pipeline takes it as `--yleaf_data`. The same command then puts the current tree in `data/ytree`: YFull's `current_tree.json` at its release commit, and the two Yleaf 4.1.6 position tables at the commit of that tag, each checked against its sha256 in `versions.env` (`YTREE_*`). An install made before the current tree existed is completed by running the command again. If the current tree cannot be downloaded, the command warns and still exits 0, because Yleaf's own data is in place. Without `data/ytree` the step prints `Current tree: not installed` and writes Yleaf's result only.

With `run-all.sh`: `TOOLS=...,y_haplogroup` (it is opt-in, so a default run lists it as skipped; without the data it is skipped with the reason).

Yleaf downloads the whole hg38 FASTA on its first run unless its config file names one, and the image's config is read-only. The launcher points Yleaf at the pipeline's reference before it starts (a BAM never needs the sequence itself), so nothing is downloaded and every container runs without a network.

## Output

| File | Contents |
|---|---|
| `y_haplogroup/${SAMPLE}_y_haplogroup.txt` | Yleaf's prediction: `Hg` (the haplogroup), `Hg_marker`, `Total_reads`, `Valid_markers` (Y markers with haplogroup information and enough reads), `QC-score` |
| `y_haplogroup/yleaf/` | Yleaf's working files: the markers it read with their alleles and its log |
| `y_haplogroup/positions.txt` | the marker positions the pileup was made at |
| `y_haplogroup/${SAMPLE}_y_tree.txt` | the placement on the current tree: `Hg` (the terminal branch, YFull's name), `Level` (its depth below the root), `Path_derived`, `Path_ancestral`, `Path_no_call` (the SNP calls of every branch on the path), `Ancestral_in_terminal`, `SNPs_called`, `Min_depth`, and `Tree` and `Positions`, the tree release and the position tables it used |
| `y_haplogroup/${SAMPLE}_y_tree_path.tsv` | one row per branch from the root to the terminal one: `SNPs` (its defining SNPs), `On_tables` (those the tables place), `Derived`, `Ancestral`, `No_call`, and the names of its ancestral SNPs |

How the placement works. A SNP is called when at least 10 reads cover it and 90% of them carry one of its two bases, derived or ancestral (Yleaf's own defaults, `-r 10` and `-b 90`); anything else is a no-call. A branch is positive when it has a derived SNP and more derived than ancestral ones, and negative when it has an ancestral SNP and no more derived than ancestral ones. The terminal branch is the deepest positive branch with no negative branch above it. A branch with no call at all neither helps nor blocks, so a branch whose SNPs no table places does not stop the walk. SNPs in the pseudoautosomal regions are left out: the analysis set masks them on chrY.

`Ancestral_in_terminal` lists the SNPs of the terminal branch that read ancestral. The sample then shares only part of that branch's block of SNPs, which is a sign of a sub-branch the tree does not have yet; it can also be a read error, so look at the depth before reading much into it. `Hg` is `NA` when no branch has a derived call.

`Hg` is `NA` when too few markers had reads for a call; the step then prints `Y haplogroup: insufficient markers` with the marker count, and both reports say "insufficient markers".

## Runtime

Not measured on a full genome; the end test runs the step on the fixture's chrY slices only. An estimate is a few minutes: Yleaf reads the pileup at about 129,000 Y marker positions, and the current tree's pileup is at about 369,000 more. The placement reads the tree and the tables in about a second.

## Notes

- Yleaf 3 places the sample on the YFull tree v10.01 and names the haplogroup by clade and defining marker, for example `R-M269`. Its `QC-score` (0 to 1) is how consistently the markers on the path to that branch agree; the default acceptance threshold is 0.95.
- Short reads cover only the unique parts of the Y; that is enough for a haplogroup at 30x.
- The two results can differ in depth and in name: YFull renames and splits branches between releases. When they disagree on the lineage itself, not just the depth, look at the per-level table: a lineage with ancestral calls on its path is the weaker one.
- On the e2e fixture (six stretches of chrY, 1.3 Mb in all) the step runs on a copy of HG002 with a male indexcov row, because indexcov reads the fixture's slices as female.
