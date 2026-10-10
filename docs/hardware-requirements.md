# Hardware and Storage Requirements

Everything you need to know about disk space, RAM, CPU, and runtime before starting. This page is the one place for the download sizes, the per-step runtimes and the totals; the other pages link here instead of repeating them.

## TL;DR

- **1 sample:** 500 GB free disk from a BAM, about 700 GB from FASTQ ([disk](#total-disk-requirements)), 16 GB RAM ([RAM](#ram-requirements)), 4+ CPU cores
- **2 samples:** 1 TB free disk, 32 GB RAM, 8+ CPU cores (recommended)
- **First-time setup downloads:** ~73 GB for the default run, ~248 GB with the optional annotation databases ([table](#shared-reference-data-one-time))
- **Time per sample:** plan for more than a day for a default `run-all.sh` from a BAM on 8 CPUs, and longer from FASTQ. DeepVariant is most of it: about 13.5 h in one observed ~30x run on 8 CPUs ([runtime](#runtime-per-step))
- **GRIDSS (step 4b, opt-in):** needs a 32 GB container on its own, so 32 GB of RAM or more

---

## Disk Space Breakdown

### Per-Sample Storage

| Data | Size | When Created | Can Delete After? |
|---|---|---|---|
| Raw FASTQ (gzipped) | 60-90 GB | You bring this | Keep (original data) |
| Sorted BAM + index | 80-120 GB | Step 2 (alignment) | After all BAM-dependent steps complete |
| VCF + index | 80-200 MB | Step 3 (variant calling) | Keep (needed by many steps) |
| Manta SV VCF | 1-5 MB | Step 4 | Keep |
| AnnotSV TSV | 25-35 MB | Step 5 | Keep |
| ClinVar hits | <1 MB | Step 6 | Keep |
| PharmCAT report | 1-5 MB | Step 7 | Keep |
| ExpansionHunter output | <1 MB | Step 9 | Keep |
| TelomereHunter output | 50-200 MB | Step 10 | Keep |
| VEP annotated VCF (`.vcf.gz`) | less than its 2-5 GB uncompressed size (the bgzipped size is not measured here) | Step 13 | Keep (comprehensive annotation) |
| CPSR report + data | 50-200 MB | Step 17 | Keep |
| CNVpytor .pytor file + calls | 5-15 GB | Step 18 | .pytor file can be deleted |
| Delly BCF + VCF | 5-20 MB | Step 19 | Keep VCF, delete BCF |
| Mito analysis output | 50-200 MB | Step 20 | Keep |
| **Subtotal per sample** | **150-250 GB** | | |
| Trimmed FASTQ pair, in `<sample>/nextflow/work` | 60-90 GB | `run-all.sh` from FASTQ: fastp | After a good run |
| Alignment intermediate (`.sam.gz`), in `work/` | 119 GB in one 30x run | `run-all.sh` from FASTQ: minimap2 | After a good run |
| Second copy of the BAM, in `work/` | 80-120 GB | `run-all.sh` from FASTQ: duplicate marking | After a good run |
| **Peak per sample from FASTQ, until `work/` is deleted** | **~410-580 GB** | | |

A run from FASTQ needs more space than a run from a BAM. Until you delete `<sample>/nextflow/work`, it holds the trimmed reads (60-90 GB, [step 1b](01b-fastp-qc.md)), the alignment intermediate (119 GB in one 30x run, larger than the gzipped FASTQ) and a second copy of the BAM (80-120 GB), because the pipeline publishes the BAM by copying it (`publish_dir_mode = 'copy'` in `nextflow.config`). Once the run succeeded and the results look right, delete `<sample>/nextflow/work`; the next run then starts from scratch.

### Shared Reference Data (One-Time)

One row per download, in GB as `wget` and `du -h` count them (1 GB = 2^30 bytes); the download sizes were read from each server on 2026-10-02 (the two step 17 rows on 2026-10-03). Each total is the sum of the rows above it; the `setup.sh` line is the sum of the rows it names. [Reference setup](00-reference-setup.md) has the commands, under a heading per database with the same size.

**Default run** (`setup.sh` downloads the first three rows and the Docker images; the VEP caches and the PCGR bundle are for steps 13 and 17, which a default run skips when they are missing):

| Resource | Download | On disk | Used by |
|---|---|---|---|
| GRCh38 no-ALT analysis set, FASTA + `.fai` (unpacked by `setup.sh`) | ~0.8 GB | ~3 GB | every step |
| ClinVar VCF + chr-renamed and pathogenic-only copies | ~0.2 GB | ~0.4 GB | step 6 |
| AnnotSV annotations | ~5 GB | ~20 GB | step 5 |
| VEP cache, release 116 | ~26 GB | ~30 GB | step 13 |
| PCGR/CPSR ref data bundle (20260620) | ~7 GB | ~7 GB | step 17 |
| VEP cache, release 115 | ~24 GB | ~28 GB | step 17 (CPSR's own VEP) |
| Docker images | ~10-15 GB | ~10-15 GB | every step |
| **Total, default run** | **~73-78 GB** | **~98-103 GB** | |
| Of which `setup.sh` downloads (first three rows + Docker images) | ~16-21 GB | ~33-38 GB | |

While a VEP cache installs, the tarball and the unpacked tree are on disk together until the tarball is deleted (`install_vep_cache` in `scripts/lib/common.sh`): about 56 GB free for release 116 and 52 GB for release 115.

**Optional** (only for the step named):

| Resource | Download | Used by |
|---|---|---|
| T1K HLA index (built from IPD-IMGT/HLA) | ~0.05 GB | step 8, built on first run |
| CNVpytor GC/mask files | ~0.09 GB | step 18 |
| pypgx bundle | ~0.4 GB | step 32 |
| Somatic resources (gnomAD AF-only VCF, panel of normals) | ~3 GB | step 29 |

### Annotation Databases (Optional, for Steps 30-31)

These databases enable deeper pathogenicity scoring via vcfanno (step 30) and variant prioritization via slivar (step 31). All are optional — the pipeline detects which are present and skips missing tracks.

| Resource | Download | License |
|---|---|---|
| CADD v1.7 whole-genome SNVs | ~81.5 GB | Non-commercial |
| CADD v1.7 gnomAD indels | ~1.2 GB | Non-commercial |
| SpliceAI masked SNV scores | ~27 GB | Academic and not-for-profit only |
| SpliceAI masked indel scores | ~64 GB | Academic and not-for-profit only |
| REVEL v1.3 | ~0.6 GB | Free for research |
| AlphaMissense | ~0.6 GB | CC BY-NC-SA 4.0 |
| gnomAD v4.1 gene constraint | ~0.09 GB | ODC-ODbL |
| **Total, annotation databases** | **~175 GB** | |

### Total Shared Data

| Scenario | Download Size |
|---|---|
| Default run | ~73-78 GB |
| Default run + annotation databases | ~248-253 GB |

### Total Disk Requirements

| Scenario | Minimum Free Space |
|---|---|
| 1 sample, core steps only (2-3-6-7) with the step scripts | 200 GB |
| 1 sample from a BAM, full pipeline | 500 GB |
| 1 sample from FASTQ with `run-all.sh`, until `<sample>/nextflow/work` is deleted | 700 GB (the peak above plus ~100 GB of default reference data) |
| 2 samples, full pipeline, `work/` deleted after each FASTQ run | 1 TB |
| 2 samples + keeping intermediates | 1.5 TB |

> **Tip:** After the pipeline completes, the single largest file is the BAM (80-120 GB per sample, the figure every page uses). If you're done with all BAM-dependent steps (4, 9, 10, 15, 16, 18, 19, 20), you can keep it as a CRAM with [step 34](34-cram-archive.md) (about half the size), or delete the BAM entirely if you keep the FASTQ (you can always re-align).

---

## RAM Requirements

Every step runs in a Docker container with a hard `--memory` limit. The two entry points set it differently. A step script passes its own `--cpus` and `--memory`. The pipeline gives each task the CPUs and memory of its label in `conf/base.config`, doubles the memory on the one retry after an out-of-memory exit, and caps both at `--max_cpus` and `--max_memory`.

| Step | Script: CPUs / memory | Pipeline: CPUs / memory |
|---|---|---|
| 2 (minimap2 alignment) | `THREADS` (8) / 32 GB | 8 / 32 GB; duplicate marking 4 / 8 GB |
| 3 (DeepVariant) | `THREADS` (8) / `DV_MEM` (32 GB) | 8 / 32 GB |
| 4 (Manta) | `THREADS` (8) / 16 GB | 8 / 32 GB |
| 6 (ClinVar screen) | 2 / 2 GB | 2 / 4 GB |
| 7 (PharmCAT) | up to 2 / 4 GB | 2 / 4 GB |
| 9 (ExpansionHunter) | `THREADS` (4) / 4 GB | 4 / 8 GB |
| 10 (TelomereHunter) | `THREADS` (4) / 4 GB | 4 / 8 GB |
| 13 (VEP) | `THREADS` (8) / 2 GB per thread, at least 8 GB | 8 / 32 GB |
| 17 (CPSR) | 4 / 8 GB | 4 / 8 GB |
| 18 (CNVpytor) | 4 / 8 GB | 8 / 32 GB |
| 19 (Delly) | 4 / 8 GB | 4 / 8 GB |

These are limits, not measured peaks. minimap2 peaked at 10 GB on the test reference (1.8 Gb, 57% of GRCh38; `scripts/02-alignment.sh`), so [step 2](02-alignment.md) plans for about 20 GB on GRCh38; the full-genome peak has not been measured.

**Minimum system RAM:** 16 GB. `run-all.sh` passes the machine's RAM as `--max_memory`, so on a 16 GB machine the 32 GB tasks run capped at about 15 GB. A task that needs more is killed (exit 137); its one retry asks for double, and the cap cuts that back to the same 15 GB. With the step scripts, lower DeepVariant's shards and memory with `THREADS` and `DV_MEM`, for example `THREADS=4 DV_MEM=12g ./scripts/03-deepvariant.sh <sample> [male|female]`. The other scripts set a fixed `--memory`; edit it in the script to change it.
**Recommended:** 32 GB or more, so the 8-CPU tasks get close to the 32 GB they ask for.
**GRIDSS (step 4b, opt-in):** its container takes 32 GB (a 28 GB Java heap plus overhead), so it needs a machine with 32 GB or more even when nothing else runs.

---

## CPU Requirements

Each pipeline task asks for the CPUs of its label in `conf/base.config`: 1, 2, 4 or 8. The 8-CPU tasks are minimap2 (index and alignment), DeepVariant, Manta, CNVpytor and VEP. `--max_cpus` and `--max_memory` cap each task; `run-all.sh` passes the machine's CPU count (or `THREADS`) and its RAM. They do not limit how much runs at once: Nextflow starts tasks until their requests fill the machine's CPUs and RAM. On a Mac with Docker Desktop, that is the Mac's CPUs and RAM, not the Docker VM's. `MAX_JOBS` is no longer read.

**DeepVariant's shards.** `run-all.sh` runs DeepVariant with 8 shards whatever the core count. More cores help the steps that run beside it. To give DeepVariant more, run `scripts/03-deepvariant.sh` with `THREADS=N` (its `--cpus` and `--num_shards`), or pass `-c` with a `withName: 'DEEPVARIANT'` block:

```groovy
// dv16.config, used as: ./scripts/run-all.sh <sample> <sex> -c dv16.config
process {
    withName: 'DEEPVARIANT' {
        cpus = 16
    }
}
```

`--max_cpus` still caps it. Its memory at 16 shards has not been measured; the request stays at 32 GB. The shard count is part of DeepVariant's command, so changing it makes `-resume` run DeepVariant again. `make_examples` runs one process per shard; upstream says `call_variants` on CPU scales sub-linearly ([DeepVariant details](https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-details.md#call_variants)).

- **Minimum:** 4 cores. DeepVariant then runs 4 shards and, extrapolated from the 8-CPU run, likely takes more than a day on its own (not measured).
- **Recommended:** 8 to 16 cores. On 8, each 8-CPU task runs alone ([what runs at once](#what-runs-at-once)); on 16, the BAM steps run beside DeepVariant. A 16-core run has not been measured end to end.

### A shared host

Each step has a hard memory limit. On the `run-all.sh` path, CPU is a Docker share (Nextflow passes `--cpu-shares`, 1024 per requested CPU), not a cap. An idle machine lends a task every core, and on a busy one the pipeline's tasks outweigh a default container (1024). Only the single-step scripts use hard `--cpus` caps. The repo sets no low-priority option.

To let the other services on the machine win, give the pipeline fewer shares and a smaller total in a config file:

```groovy
// shared-host.config, used as:
//   THREADS=6 ./scripts/run-all.sh <sample> <sex> --max_memory 24.GB -c shared-host.config
// This value replaces the docker profile's container options, so it repeats --network none.
process.containerOptions = '--network none --cpu-shares 256'
// Everything running at once stays within this. Keep it at least --max_cpus and
// --max_memory, or Nextflow refuses the tasks that ask for more.
executor {
    cpus   = 6
    memory = 24.GB
}
```

Check the shares on a running task: `docker ps` lists the pipeline's containers as `nxf-...`, and `docker inspect -f '{{.HostConfig.CpuShares}}' <container>` should print 256. The same `executor` block is the way to keep Docker Desktop's VM from being overbooked: set it to the CPUs and memory the VM has.

### Runtime per step

In one observed ~30x run on an 8-CPU, 24 GB budget, DeepVariant 1.10 took about 13.5 hours: `make_examples` about 10.5 h for 5.8 million candidates, `call_variants` about 3 h, `postprocess_variants` under an hour. On a busy host `make_examples` ran at about half that speed. Upstream reports 1 h 9 min on a 96-vCPU cloud machine ([DeepVariant r1.10 metrics](https://github.com/google/deepvariant/blob/r1.10/docs/metrics.md)), about 110 vCPU-hours, which agrees. Plan for more than a day for a full run from a BAM on 8 CPUs, and longer from FASTQ. We have not measured a 16-core run end to end.

The other rows are estimates from the step pages for a 16-core / 32 GB desktop with each script's default CPU limit, not measurements. Step pages link here instead of giving their own figure.

| Step | Runtime | Notes |
|---|---|---|
| 1 ORA to FASTQ | ~30 min | only for Illumina ORA input |
| 1b fastp | ~10-20 min | |
| 2 Alignment (minimap2) | ~1-2 h | plus ~30 min once to build the `.mmi` index |
| 3 DeepVariant | ~13.5 h on 8 CPUs | observed once, ~30x, DeepVariant 1.10; 1 h 9 min on 96 vCPUs upstream |
| 4 Manta | ~20 min | |
| 5 AnnotSV | ~10 min | |
| 6 ClinVar screen | ~5 min | |
| 7 PharmCAT | ~10 min | |
| 8 HLA typing (T1K) | ~30 min | plus well under a minute once to build the index (17 s in [CI run 37980248545](https://github.com/GeiserX/Personal-Genome-Pipeline/actions/runs/37980248545)) |
| 9 ExpansionHunter | ~15 min | |
| 9b Stranger | < 1 min | |
| 10 TelomereHunter | ~1 h | |
| 11 ROH | ~5 min | |
| 12 Mito haplogroup | ~1 min | |
| 13 VEP | ~2-4 h | |
| 14 Imputation prep (opt-in) | ~10 min | |
| 15 duphold | ~20 min | |
| 16 indexcov | ~5 s | |
| 16b mosdepth | ~5-10 min | |
| 17 CPSR | ~30-60 min | |
| 18 CNVpytor | ~1-3 h | |
| 19 Delly | ~2-4 h | |
| 20 Mito variants (Mutect2) | ~15-30 min | |
| 21 Cyrius | ~5-15 min | opt-in; `setup.sh --cyrius` installs it once (about a minute) |
| 22 SV consensus | ~5-15 min | |
| 23 Clinical filter | ~5-10 min | |
| 24 HTML report | ~1-3 min | |
| 25 PRS | ~20-40 min | |
| 26 Ancestry (opt-in) | ~15-30 min | most of it the first-run download |
| 27 CPIC lookup | ~1-2 min | |
| 28 MultiQC | < 1 min | |
| 29 Somatic (opt-in) | ~2-6 h | under 5 min on a few genes with `INTERVALS` |
| 30 vcfanno | ~5-15 min | |
| 31 slivar | ~5-10 min | |
| 32 pypgx | ~20-40 min | |
| 3d Octopus (alternative) | ~2-4 h | |
| 4b GRIDSS (opt-in) | ~4-8 h | |

| Run | Runtime |
|---|---|
| DeepVariant alone, 8 CPUs | about 13.5 h (one observed ~30x run) |
| DeepVariant alone, 96 vCPUs | 1 h 9 min (upstream r1.10 metrics) |
| DeepVariant alone, 4 CPUs | likely more than a day (extrapolated, not measured) |
| Default `run-all.sh` from a BAM, 8 CPUs | plan for more than a day |
| Default `run-all.sh` from FASTQ, 8 CPUs | longer than from a BAM |
| Default `run-all.sh`, 16 CPUs | not measured end to end |

### What runs at once

From FASTQ, the pipeline trims, aligns and marks duplicates first. indexcov then checks the sex from the BAM index, and from that moment DeepVariant and every step that reads the BAM are ready together. The steps that read the VCF wait for DeepVariant.

```
FASTQ ──> fastp ──> minimap2 ──> duplicate marking ──> BAM ──> indexcov (sex check)
                                                                   │
        ┌──────────────────────────────────────────────────────────┤
        ▼                                                          ▼
DeepVariant (8 CPUs)                      BAM steps: Manta, Delly, CNVpytor, ExpansionHunter,
        │                                 TelomereHunter, HLA, mito variants, mosdepth, pypgx
        ▼                                                          │
VCF steps: ClinVar, PharmCAT, VEP, CPSR,                           ▼
ROH, PRS, vcfanno, slivar                 SV chain: duphold, AnnotSV, SURVIVOR merge
        │                                                          │
        └─────────────────> HTML report, MultiQC <─────────────────┘
```

Nextflow starts a task only when its CPUs and memory fit in what is left of the machine's. On a host with 8 CPUs, every 8-CPU task (minimap2, DeepVariant, Manta, CNVpytor, VEP) runs alone, and the light SV chain steps wait for it: about 18 h behind DeepVariant in the observed run, after which they finished in seconds. On 16 CPUs, DeepVariant leaves 8 free for the BAM steps.

---

## Internet Bandwidth

### One-Time Downloads

About 73 GB for a default run and 248 GB with the annotation databases; the [table above](#shared-reference-data-one-time) has every download. The VEP caches come from slow servers, so use `wget -c` to resume.

### Ongoing Downloads

- **ClinVar updates:** ~200 MB/month (optional but recommended for latest pathogenic variant classifications)
- **Docker image updates:** Variable (only when you want newer tool versions)

> **Network during a run:** after setup, a few scripts still fetch public files (PGS scoring files, and the data of the opt-in steps they are asked for); no step container has network. None of these calls sends sample data. [Why run locally?](why-local.md#network-calls-during-a-run) lists each one.

---

## Storage Tips

### Save Disk Space

1. **Keep the BAM as a CRAM** after all BAM-dependent steps complete: `./scripts/34-cram-archive.sh <sample>` writes the CRAM, checks it against the BAM, and with `--delete-bam` then deletes the BAM ([step 34](34-cram-archive.md)). The CRAM is about half the size of the BAM, and it can only be read with the same reference FASTA, so keep that file.

2. **Delete intermediate files:**
   - `<sample>/nextflow/work` once a `run-all.sh` run succeeded and the results look right (the next run then starts from scratch)
   - CNVpytor `.pytor` files (5-15 GB each)
   - Delly `.bcf` files (after converting to VCF)

3. **Keep your original FASTQ or ORA files.** A BAM from a default, trimmed run does not hold your original reads: fastp cut bases and dropped reads before alignment, so FASTQ extracted from that BAM is only what fastp kept. Delete the BAM rather than the FASTQ if you need the space; you can re-align.

### Storage Medium Recommendations

| Medium | Suitable For | Notes |
|---|---|---|
| NVMe SSD | Active analysis | Fastest. 10-50x faster than HDD for random reads. |
| SATA SSD | Active analysis | Good performance. Adequate for all pipeline steps. |
| HDD (7200 RPM) | Storage / archive | Adequate for sequential I/O (alignment, VEP). Random access steps (DeepVariant) will be slower. |
| Network storage (NFS/SMB) | Archive only | Too slow for active analysis. Use for long-term storage after pipeline completes. |
| USB external drive | Emergency only | Severely bottlenecks I/O-intensive steps. |

---

## Cloud Cost Comparison

If you don't have suitable hardware, cloud instances work well:

| Provider | Instance | vCPUs | RAM | Cost/hr |
|---|---|---|---|---|
| AWS | c5.4xlarge | 16 | 32 GB | ~$0.68 |
| GCP | n2-standard-16 | 16 | 64 GB | ~$0.78 |
| Azure | Standard_D16s_v5 | 16 | 64 GB | ~$0.77 |
| Hetzner | CCX33 | 8 | 32 GB | ~$0.18 |

Multiply the hourly price by the run time: plan for more than a day on 8 vCPUs ([runtime](#runtime-per-step)). On 16 vCPUs `run-all.sh` still gives DeepVariant 8 shards, and a 16-vCPU run has not been measured end to end. Add ~$0.10/GB/month for persistent disk storage. A 500 GB disk costs ~$50/month.

> **Tip:** Spot or preemptible instances cost less, but a preemption stops the run. `-resume` reuses every task that finished and starts the running one over, and DeepVariant is a single task of about 13.5 h on 8 CPUs (one observed run).
