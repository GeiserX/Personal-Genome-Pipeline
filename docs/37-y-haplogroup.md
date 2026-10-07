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
./scripts/37-y-haplogroup.sh your_name
```

With `run-all.sh`: `TOOLS=...,y_haplogroup` (it is opt-in, so a default run lists it as skipped).

Yleaf downloads the whole hg38 FASTA on its first run unless its config file names one, and the image's config is read-only. The launcher points Yleaf at the pipeline's reference before it starts (a BAM never needs the sequence itself), so nothing is downloaded and every container runs without a network.

## Output

| File | Contents |
|---|---|
| `y_haplogroup/${SAMPLE}_y_haplogroup.txt` | Yleaf's prediction: `Hg` (the haplogroup), `Hg_marker`, `Total_reads`, `Valid_markers` (Y markers with haplogroup information and enough reads), `QC-score` |
| `y_haplogroup/yleaf/` | Yleaf's working files: the markers it read with their alleles and its log |
| `y_haplogroup/positions.txt` | the marker positions the pileup was made at |

`Hg` is `NA` when too few markers had reads for a call; the step then prints `Y haplogroup: insufficient markers` with the marker count, and both reports say "insufficient markers".

## Runtime

A few minutes: Yleaf reads the pileup at about 129,000 Y marker positions.

## Notes

- Yleaf 3 places the sample on the YFull tree (v10.01) and names the haplogroup by clade and defining marker, for example `R-M269`. Its `QC-score` (0 to 1) is how consistently the markers on the path to that branch agree; the default acceptance threshold is 0.95.
- Short reads cover only the unique parts of the Y; that is enough for a haplogroup at 30x.
- On the e2e fixture (300 kb of chrY) the step runs on a copy of HG002 with a male indexcov row, because indexcov reads the fixture's slices as female.
