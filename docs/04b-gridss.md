# Step 4b: GRIDSS Structural Variant Calling

Assembly-based structural variant caller that identifies breakpoints through local de novo assembly. Excels at complex rearrangements, small SVs, and events with microhomology that Manta and Delly miss.

---

## What It Does

1. **Extracts SV-supporting reads** — soft-clipped, discordant, and split reads
2. **Performs local assembly** — builds contigs across breakpoints using de Bruijn graphs
3. **Calls breakpoints** — identifies both sides of each structural variant
4. **Scores variants** — probabilistic model incorporating read pair, split read, and assembly evidence

## Why

GRIDSS complements Manta and Delly in the SV consensus pipeline (step 22):
- **Complex rearrangements** — GRIDSS assembles across breakpoints rather than relying on read pair signals alone
- **Small SVs (50-300bp)** — too small for discordant read pair detection but too large for standard SNP callers
- **Breakpoint precision** — assembly provides single-nucleotide resolution at breakpoints
- **Microhomology detection** — identifies shared sequence at breakpoint junctions

## Tool

**GRIDSS** v2.13.2 — Genomic Rearrangement IDentification Software Suite.

- Paper: Cameron et al., Genome Research 2017 (doi:10.1101/gr.222109.117)
- Source: [github.com/PapenfussLab/gridss](https://github.com/PapenfussLab/gridss)

## Docker Image

- `GRIDSS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

Image size: ~1.5 GB (includes Java 11, R, bwa, samtools, all dependencies).

## Command

```bash
export GENOME_DIR=/path/to/data
./scripts/04b-gridss.sh <sample_name>
```

`run-all.sh` runs GRIDSS only when asked, because it needs a classic BWA index that `setup.sh` does not build and about 32 GB of RAM:

```bash
GRIDSS=true ./scripts/run-all.sh <sample_name> <male|female>
```

With `GRIDSS=true` and no classic BWA index, the step fails and the run's final table lists it as failed.

Before it starts, the step reads Docker's memory limit (`docker info`). Below 32 GB it prints `SKIPPED` with the amount and exits without running GRIDSS, which otherwise fails late without a clear message. `GRIDSS_MIN_MEM_GB` changes the threshold if you want to try with less. `THREADS` (default 8) sets GRIDSS's `-t` and the container's CPUs.

## Prerequisites

GRIDSS requires a **classic BWA index** (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) alongside the reference FASTA. **BWA-MEM2 index files (`.bwt.2bit.64`) are NOT compatible** — GRIDSS bundles classic `bwa` internally for its read realignment step. Generate the classic index if you don't have one:

```bash
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
source versions.env   # from the repository root: BWA_IMAGE, the image step 04b uses for this
docker run --rm -v "${GENOME_DIR}:/genome" \
  "${BWA_IMAGE}" \
  bwa index "/genome/${REF_FASTA}"
```

This takes ~1 hour and produces 5 index files (~5 GB total). Only needed once.

## Output

| File | Description |
|---|---|
| `sv_gridss/<sample>_gridss.vcf.gz` | SV calls in BND notation |
| `sv_gridss/<sample>_assembly.bam` | Assembly contigs (intermediate, can be deleted) |
| `sv_gridss/work/` | GRIDSS's working directory (`--workingdir`), about 50 GB during the run; removed once the VCF is written, kept when GRIDSS fails so you can look at it |

### BND notation

GRIDSS reports all variants as **BND** (breakend) records, not DEL/DUP/INV/INS. Each structural variant produces two VCF records (one for each breakpoint). Nothing here converts them to DEL, DUP or INV, so step 22 leaves GRIDSS out of the SV consensus (see Notes).

### Quality filtering

| Confidence | Filter |
|---|---|
| High | `QUAL >= 1000` and `AS > 0` and `RAS > 0` |
| Medium | `QUAL >= 500` and (`AS > 0` or `RAS > 0`) |

Note: GRIDSS QUAL scores are uncorrected for multiple testing and tend to be overestimated.

## Resource Requirements

| Resource | Requirement |
|---|---|
| Memory | 32 GB (28 GB JVM heap + OS overhead) |
| CPU | 8 threads |
| Disk | ~50 GB intermediate files in `sv_gridss/work/`, removed after a successful run |
| Runtime | 4-8 hours for 30X WGS |

GRIDSS is the heaviest tool in the pipeline. With `GRIDSS=true` it runs in parallel with the other heavy steps (VEP, CNVpytor, Delly).

## Runtime

| Dataset | Threads | Time |
|---|---|---|
| 30X WGS | 8 | 4-8 hours |
| chr22 BAM | 4 | ~10-20 min |

## Known Issues

- **BWA index required** — GRIDSS performs its own BWA realignment of soft-clipped reads. Without a BWA index, it fails
- **Intermediate files in /tmp** — GRIDSS uses Java's temp directory for sorts. If `/tmp` is small, set `TMP_DIR` in the container
- **Memory management** — specifying JVM heap between 32-48 GB is counterproductive (crosses Java's compressed oops threshold). Use either <=31 GB or >=49 GB
- **ENCODE blacklist** — the script auto-downloads the hg38 ENCODE blacklist (ENCFF356LFX.bed, 910 regions) to filter known problematic regions

## Notes

- GRIDSS output is not part of the SV consensus (step 22): GRIDSS reports every event as a pair of breakends (`SVTYPE=BND`), which never match the DEL, DUP and INV records of the other callers
- The assembly BAM (`_assembly.bam`) can be deleted after the VCF is generated to save disk space
- Supports `ALIGN_DIR` variable: `ALIGN_DIR=aligned_bwamem2 ./scripts/04b-gridss.sh sample`
- GRIDSS does not have ARM64 support — runs under emulation on Apple Silicon (very slow)
