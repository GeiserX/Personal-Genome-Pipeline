# Step 16b: mosdepth Coverage Statistics

Fast per-base and per-region depth statistics from BAM. Produces coverage distributions, threshold reports, and a genome-wide summary.

---

## What It Does

Reads the BAM file and computes:
1. **Mean coverage** per chromosome and genome-wide
2. **Coverage distribution** — cumulative fraction of bases at each depth level
3. **Region coverage** — mean depth in 500bp windows across the genome
4. **Threshold report** — fraction of bases at 1x, 5x, 10x, 15x, 20x, 30x, 50x depth
5. **Gene callability** — for each gene of `assets/callability_genes.tsv` on chr1-22, X and Y, the bases with 10 or more reads of mapping quality 20 or more (a second mosdepth pass)

## Why

- **Validates sequencing depth** — confirms your WGS actually has 30X (or whatever depth was ordered)
- **Identifies low-coverage regions** — gaps in coverage cause missed variant calls
- **Detects sex chromosome anomalies** — chrX depth relative to autosomes reveals XX vs XY
- **Says which hard genes this genome could not read** — in a gene inside a segmental duplication, a false duplication of GRCh38 or on an ALT contig, short reads often get mapping quality 0, so no caller sees them and a "no finding" there means nothing. The report lists those genes (see [Gene callability](#gene-callability))
- **More precise than indexcov** — indexcov (step 16) reads only the BAM index for a quick estimate; mosdepth reads actual alignments

## Tool

**mosdepth** v0.3.14 — fast BAM/CRAM depth calculation.

- Paper: Pedersen & Quinlan, Bioinformatics 2018 (doi:10.1093/bioinformatics/btx699)
- Source: [github.com/brentp/mosdepth](https://github.com/brentp/mosdepth)

## Docker Image

- `MOSDEPTH_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command

```bash
export GENOME_DIR=/path/to/data
./scripts/16b-mosdepth.sh <sample_name>
```

## Output

| File | Description |
|---|---|
| `mosdepth/<sample>.mosdepth.summary.txt` | Per-chromosome mean coverage + genome total |
| `mosdepth/<sample>.mosdepth.global.dist.txt` | Cumulative coverage distribution (for plotting) |
| `mosdepth/<sample>.regions.bed.gz` | Mean depth per 500bp window |
| `mosdepth/<sample>.thresholds.bed.gz` | Fraction of bases at 1x/5x/10x/15x/20x/30x/50x |
| `mosdepth/<sample>_callability.bed.gz` | Per gene of `assets/callability_genes.tsv` with a position: chrom, start, end, gene, and the bases with 10 or more reads of MAPQ 20 or more (mosdepth's thresholds table, column `10X`) |

The Nextflow MOSDEPTH task writes the same files into `<outdir>/<sample>/coverage/`.

## Interpreting Results

### Summary file

```
chrom   length      bases          mean    min  max
chr1    248956422   7394254892     29.70   0    312
chr2    242193529   7185433827     29.67   0    290
...
total   3088286401  91784821430    29.72   0    312
```

- **mean ~30**: You have 30X whole-genome coverage (good for germline calling)
- **mean < 15**: Variant calling accuracy degrades significantly
- **chrX mean ~15 (male) or ~30 (female)**: Expected sex-linked depth

### Threshold report

Shows what fraction of each region is covered at key depths:
- **>95% at 10x**: Adequate for most germline variant callers
- **>90% at 20x**: Good quality WGS
- **>80% at 30x**: High-quality WGS

### Gene callability

`assets/callability_genes.tsv` lists 1,572 genes where short reads are often not placed with confidence. `scripts/ci/build-callability.sh` writes it from pinned public files, and a CI case rebuilds it and fails on any difference:

| Reason | Genes | Source |
|---|---|---|
| `CMRG` | 273 | GIAB's challenging medically relevant genes, CMRG v1.00 (Wagner et al., Nat Biotechnol 2022, PMID 35132260) |
| `segdup` | 1,276 | protein-coding genes with half or more of the gene body in a segmental duplication (GIAB genome stratifications v3.6, GRCh38) |
| `false_duplication` | 22 | genes over a false duplication of GRCh38, either copy: KCNE1, CBS and CRYAA among them (Behera et al., Genome Biol 2023, PMID 36810122; the same GIAB stratifications) |
| `collapsed` | 50 | genes over a collapsed duplication: two copies in people, one in GRCh38 (GIAB stratifications) |
| `ALT-only` | 19 | protein-coding genes only on an ALT scaffold of GRCh38, so the default no-ALT reference cannot hold them: GSTT1, HLA-DRB3, HLA-DRB4 and others (GENCODE 50, NCBI's GRCh38 assembly report) |

A gene can have several reasons. Gene positions are GENCODE 50 gene spans. The `owner` column names the step that calls the gene with a method of its own: `parascopy` (step 35, the genes of `assets/paralog_loci.tsv`), `pypgx` (step 32, the bam genes of `assets/pgx_genes.tsv`) or `cn_loci:<locus>` (the depth engine, `bin/cn_loci.tsv`); `.` for none.

A gene counts as **read** when 90% or more of its bases have 10 or more reads of MAPQ 20 or more. Both reports have a section "Genes this genome could not read reliably" in the Quality chapter. It lists every gene with no owner that was not read, with the share of its bases that reached that depth, and every ALT-only gene with no owner as "not on this reference". Genes with an owner that were not read are listed apart: their own step's section says what it called. On a partial BAM (the CI fixture, an exome) most genes are not read, because they have no reads at all.

`CALLABILITY_GENES` names another table with the same columns (the CI case uses a two-gene one). A table from an older catalog, one that lacks a gene of the current catalog, makes the section unreadable: run step 16b again. When only `<sample>_callability.bed.gz` is missing, step 16b runs the gene pass alone.

The CMRG gene list and the stratification BEDs are only inputs of the build script; a run needs only the committed table.

## Runtime

| Dataset | Threads | Time | Memory |
|---|---|---|---|
| 30-40x WGS, measured in three runs ([runtime](hardware-requirements.md#runtime-per-step)) | 4 | 2-3 min | 1.9 GB |
| The gene callability pass on the same BAM (estimate: a second read of the BAM) | 4 | about as long again | < 2 GB |
| chr22 BAM (estimate) | 2 | < 30 sec | < 1 GB |

## MultiQC Integration

mosdepth output is automatically detected by MultiQC. The `mosdepth.global.dist.txt` and `mosdepth.summary.txt` files are consumed to generate coverage distribution plots and summary statistics in the aggregated report.

## Notes

- `--fast-mode` skips per-base output (saves ~10 GB of disk for 30X WGS) while keeping all distributions and thresholds
- Supports `ALIGN_DIR` variable to use alternative alignments: `ALIGN_DIR=aligned_bwamem2 ./scripts/16b-mosdepth.sh sample`
- Supports `CAPTURE_BED` for WES on-target coverage: `CAPTURE_BED=${GENOME_DIR}/captures/my_panel.bed ./scripts/16b-mosdepth.sh sample`. **The BED file must be inside `GENOME_DIR`** — the Docker container only mounts `GENOME_DIR` as `/genome/`, so paths outside it are invisible to the container
- mosdepth uses 1 main thread + N decompression threads, so `--threads 4` actually uses 5 threads total
