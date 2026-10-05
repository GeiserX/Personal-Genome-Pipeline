# Step 34: CRAM Archive

## What This Does
Writes the sample's alignments as CRAM beside the BAM, checks that the CRAM holds exactly the same reads, and, only when you ask and only after that check passed, deletes the BAM. `--restore` writes the BAM back from the CRAM.

## Why
A 30x genome keeps 30 to 80 GB of BAM after the analysis is done. A CRAM stores each read as its difference from the reference and is about half that size or less (on the CI fixture, 68 MB of CRAM for 133 MB of BAM), with nothing lost: every read, base quality, flag and tag is still there. Most people keep their alignments to rerun a step when a tool or database improves; a CRAM keeps that option for half the disk.

The price is the reference: a CRAM can only be read with the same FASTA it was written with. Keep `reference/GRCh38_no_alt_analysis_set.fasta` (or your `REF_FASTA`) as long as you keep the CRAM. samtools checks each contig's MD5 when it reads a CRAM and stops with an error on a different file, so a wrong reference cannot give silently wrong reads.

## Tool
- **samtools** (`view -C`, `index`, `quickcheck`, `flagstat`)

## Docker Image
- `SAMTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
./scripts/34-cram-archive.sh your_name               # write and check the CRAM; the BAM is kept
./scripts/34-cram-archive.sh your_name --delete-bam  # the same, then delete the BAM and its index
./scripts/34-cram-archive.sh your_name --restore     # write the BAM back from the CRAM
```

The script runs, in the samtools container:
```bash
samtools view -C --reference "${REF_FASTA}" -o aligned/${SAMPLE}_sorted.part.cram aligned/${SAMPLE}_sorted.bam
samtools quickcheck -v aligned/${SAMPLE}_sorted.part.cram
samtools index aligned/${SAMPLE}_sorted.part.cram
samtools flagstat aligned/${SAMPLE}_sorted.bam
samtools flagstat --input-fmt-option reference="${REF_FASTA}" aligned/${SAMPLE}_sorted.part.cram
```

The check is that the CRAM passes `samtools quickcheck` (a CRAM cut short has no end-of-file block) and that the two `flagstat` outputs are equal line for line: the same number of reads, mapped, paired, properly paired, duplicates and secondary and supplementary alignments. Only then does the CRAM lose its `.part` name. When the check fails, the CRAM is removed, the BAM is kept and the step exits non-zero, also with `--delete-bam`. An older CRAM is always written again rather than trusted, because the BAM may have been aligned again since.

`--restore` does the same in reverse: the BAM is written under a `.part` name, checked against the CRAM, and moved into place. It refuses to write over a BAM that exists.

In Nextflow, `cram_archive` in `--tools` runs the same conversion and check (`CRAM_ARCHIVE`) for every BAM of the run and publishes `<outdir>/<sample>/aligned/<sample>_sorted.cram`. It never deletes anything: delete the BAM yourself once the CRAM is there, or run this script on it.

## Output Files
| File | Description |
|---|---|
| `aligned/${SAMPLE}_sorted.cram` | The alignments |
| `aligned/${SAMPLE}_sorted.cram.crai` | Its index |
| `aligned/${SAMPLE}_sorted.cram.flagstat`, `aligned/${SAMPLE}_sorted.bam.flagstat` | The two `flagstat` outputs the check compared |

## Rerunning steps from a CRAM

| How you run | From a CRAM |
|---|---|
| Nextflow | Every step. Put the CRAM in the samplesheet's `cram` and `crai` columns (with the `vcf` columns, or a `sex` column to call it again) and pass the reference it was written with as `--reference`. `CRAM_TO_BAM` writes it out once as a BAM in the work directory, checked against the CRAM, and every BAM step reads that. The BAM needs the disk space of the original until the run's work directory is removed |
| The bash steps that read alignments | None directly: steps 03, 04, 04b, 08, 09, 10, 15, 16, 16b, 18, 19, 20, 21, 28 (its samtools flagstat), 29, 32 and 33, and the alternative callers 03a to 03d and 04a, read `aligned/${SAMPLE}_sorted.bam` (03e and 04c read `aligned_longread/`, which this step does not archive). Run `./scripts/34-cram-archive.sh your_name --restore` first, then the step |
| VCF-only steps | Every one: steps 06, 07, 11, 12, 13, 14, 17, 22, 23, 25, 26, 27, 30 and 31 read the VCF, which this step never touches |

## Runtime
Writing a CRAM reads the whole BAM and compresses it again; the check reads both files once more. On CI's fixture slice it takes seconds; a whole genome has not been timed here.

## Notes
- Keep the reference FASTA (and its `.fai`) with the CRAM. Another build of GRCh38, even one with the same contig names, cannot decode it.
- The CRAM is written with samtools' default CRAM version and compression, which every current htslib, Picard and GATK read.
- A CRAM made elsewhere (nf-core/sarek writes them) can go straight into the Nextflow samplesheet; for the bash steps, `--restore` needs it at `aligned/${SAMPLE}_sorted.cram` with its `.crai`.
