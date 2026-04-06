# Step 16b: mosdepth Coverage Statistics

Fast per-base and per-region depth statistics from BAM. Produces coverage distributions, threshold reports, and a genome-wide summary.

---

## What It Does

Reads the BAM file and computes:
1. **Mean coverage** per chromosome and genome-wide
2. **Coverage distribution** — cumulative fraction of bases at each depth level
3. **Region coverage** — mean depth in 500bp windows across the genome
4. **Threshold report** — fraction of bases at 1x, 5x, 10x, 15x, 20x, 30x, 50x depth

## Why

- **Validates sequencing depth** — confirms your WGS actually has 30X (or whatever depth was ordered)
- **Identifies low-coverage regions** — gaps in coverage cause missed variant calls
- **Detects sex chromosome anomalies** — chrX depth relative to autosomes reveals XX vs XY
- **More precise than indexcov** — indexcov (step 16) reads only the BAM index for a quick estimate; mosdepth reads actual alignments

## Tool

**mosdepth** v0.3.13 — fast BAM/CRAM depth calculation.

- Paper: Pedersen & Quinlan, Bioinformatics 2018 (doi:10.1093/bioinformatics/btx699)
- Source: [github.com/brentp/mosdepth](https://github.com/brentp/mosdepth)

## Docker Image

```
quay.io/biocontainers/mosdepth:0.3.13--hba6dcaf_0
```

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

## Runtime

| Dataset | Threads | Time | Memory |
|---|---|---|---|
| 30X WGS (~100 GB BAM) | 4 | ~5-10 min | < 2 GB |
| chr22 BAM | 2 | < 30 sec | < 1 GB |

## MultiQC Integration

mosdepth output is automatically detected by MultiQC. The `mosdepth.global.dist.txt` and `mosdepth.summary.txt` files are consumed to generate coverage distribution plots and summary statistics in the aggregated report.

## Notes

- `--fast-mode` skips per-base output (saves ~10 GB of disk for 30X WGS) while keeping all distributions and thresholds
- Supports `ALIGN_DIR` variable to use alternative alignments: `ALIGN_DIR=aligned_bwamem2 ./scripts/16b-mosdepth.sh sample`
- Supports `CAPTURE_BED` for WES on-target coverage: `CAPTURE_BED=${GENOME_DIR}/captures/my_panel.bed ./scripts/16b-mosdepth.sh sample`. **The BED file must be inside `GENOME_DIR`** — the Docker container only mounts `GENOME_DIR` as `/genome/`, so paths outside it are invisible to the container
- mosdepth uses 1 main thread + N decompression threads, so `--threads 4` actually uses 5 threads total
