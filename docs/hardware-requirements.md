# Hardware and Storage Requirements

Everything you need to know about disk space, RAM, CPU, and runtime before starting. This page is the one place for the download sizes, the per-step runtimes and the totals; the other pages link here instead of repeating them.

## TL;DR

- **1 sample:** 500 GB free disk from a BAM, about 700 GB from FASTQ ([disk](#total-disk-requirements)), 16 GB RAM ([RAM](#ram-requirements)), 4+ CPU cores
- **2 samples:** 1 TB free disk, 32 GB RAM, 8+ CPU cores (recommended)
- **First-time setup downloads:** ~73 GB for the default run, ~248 GB with the optional annotation databases ([table](#shared-reference-data-one-time))
- **Time per sample:** a default `run-all.sh` from FASTQ took about 20-21 h on a 16-thread server, measured on 30-40x genomes; on 8 CPUs plan for more than a day (not measured). DeepVariant is most of it: 8 h 50 min to 13 h 30 min ([runtime](#runtime-per-step))
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
| Second copy of the BAM, in `work/` (none when published by hard link) | 80-120 GB | `run-all.sh` from FASTQ: alignment | After a good run |
| **Peak per sample from FASTQ, until `work/` is deleted (estimate)** | **~290-460 GB** | | |

A run from FASTQ needs more space than a run from a BAM. Until you delete `<sample>/nextflow/work`, it holds the trimmed reads (60-90 GB, [step 1b](01b-fastp-qc.md)) and a second copy of the BAM (80-120 GB) when the BAM is published by copying (`publish_dir_mode = 'copy'`, the default in `nextflow.config`). The alignment writes no SAM file: minimap2 is piped into `samtools fixmate`, `sort` and `markdup` in one task. The peak is an estimate summed from the rows above, not a measurement: the subtotal (150-250 GB) plus the trimmed reads (60-90 GB) plus the second BAM (80-120 GB) is 290-460 GB. It leaves out the temporary files `samtools sort` writes in the task's work directory while it aligns; they are removed when the task ends, and their size has not been measured here. `run-all.sh` passes `--publish_dir_mode link`, so the BAM is stored once, when a probe link from `work/` to the output directory works and it runs as root or Linux's `fs.protected_hardlinks` is 0. It chooses nothing when you pass `--publish_dir_mode` (your mode is used) or `-w`, `-work-dir`, `-c` or `-config`, or set `NXF_WORK` (the copy default stays, unless your config sets another mode). Once the run succeeded and the results look right, delete `<sample>/nextflow/work`; the next run then starts from scratch.

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
| Fine-grained ancestry panel (gnomAD HGDP + 1000 Genomes, phased) | ~19 GB read, a few GB kept | step 26b, `setup.sh --ancestry-fine-panel` |
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
| 1 sample from FASTQ with `run-all.sh`, until `<sample>/nextflow/work` is deleted | 700 GB (the estimated peak above plus ~100 GB of default reference data is ~560 GB; the rest is room for the sort's temporary files, not measured) |
| 2 samples, full pipeline, `work/` deleted after each FASTQ run | 1 TB |
| 2 samples + keeping intermediates | 1.5 TB |

> **Tip:** After the pipeline completes, the single largest file is the BAM (80-120 GB per sample, the figure every page uses). If you're done with all BAM-dependent steps (4, 9, 10, 15, 16, 18, 19, 20), you can keep it as a CRAM with [step 34](34-cram-archive.md) (about half the size), or delete the BAM entirely if you keep the FASTQ (you can always re-align).

---

## RAM Requirements

Every step runs in a Docker container with a hard `--memory` limit. The two entry points set it differently. A step script passes its own `--cpus` and `--memory`. The pipeline gives each task the CPUs and memory that `conf/base.config` gives its label, or its name where a step's request follows its script, doubles the memory on the one retry after an out-of-memory exit, and caps both at `--max_cpus` and `--max_memory`.

| Step | Script: CPUs / memory | Pipeline: CPUs / memory |
|---|---|---|
| 2 (minimap2 alignment) | `THREADS` (8) / 32 GB | 8 / 32 GB (minimap2 piped into samtools, one task) |
| 3 (DeepVariant) | `THREADS` (8) / `DV_MEM` (32 GB) | 8 / 32 GB |
| 4 (Manta) | `THREADS` (8) / 16 GB | 8 / 16 GB |
| 6 (ClinVar screen) | 2 / 2 GB | 2 / 4 GB |
| 7 (PharmCAT) | up to 2 / 4 GB | 2 / 4 GB |
| 9 (ExpansionHunter) | `THREADS` (4) / 4 GB | 4 / 8 GB |
| 10 (TelomereHunter) | `THREADS` (4) / 4 GB | 4 / 8 GB |
| 13 (VEP) | `THREADS` (8) / 2 GB per thread, at least 8 GB | 8 / 16 GB |
| 17 (CPSR) | 4 / 8 GB | 4 / 8 GB |
| 18 (CNVpytor) | 4 / 8 GB | 4 / 8 GB |
| 19 (Delly) | 4 / 8 GB | 4 / 8 GB |
| 25 (PRS): scoring files written as pgsc_calc reads them | no limit | `PRS_PREPARE` 1 / 16 GB |
| 25 (PRS): score positions genotyped from the gVCF | 2 / 8 GB | `PRS_SCORE_SITES` 2 / 4 GB |
| 25 (PRS): pgsc_calc, each of its tasks | `PGSC_MAX_MEMORY` (three quarters of the RAM) | `PRS` 4 / 24 GB, given to pgsc_calc as `--max_memory` |

These are limits, not measured peaks. The peaks measured for each task in the end test are in the [runtime table](#runtime-per-step): on 30-40x genomes the alignment task (minimap2 and samtools together) peaked at 18.3-18.4 GB and DeepVariant at 15.9-16.5 GB. For step 25, on the nine default scores, `PRS_PREPARE` peaked at 11.5 GB and pgsc_calc's first task, run on one thread, needed about 14 to 18 GB on a GitHub runner without a genome, and 12.7 to 13.5 GiB on the three end-test genomes ([step 25, Memory](25-prs.md#memory)). Your own run records pgsc_calc's task peaks in `prs/pgsc_calc/pipeline_info/execution_trace.txt`; the `PRS` row of the pipeline's trace is only the launcher.

**Minimum system RAM:** 16 GB. `run-all.sh` passes the machine's RAM as `--max_memory`, so on a 16 GB machine the 32 GB tasks run capped at about 15 GB. A task that needs more is killed (exit 137); its one retry asks for double, and the cap cuts that back to the same 15 GB. With the step scripts, lower DeepVariant's shards and memory with `THREADS` and `DV_MEM`, for example `THREADS=4 DV_MEM=12g ./scripts/03-deepvariant.sh <sample> [male|female]`. The other scripts set a fixed `--memory`; edit it in the script to change it. Step 25's pgsc_calc is the exception that does not fit: its first task needs about 14 to 18 GB on the default scores and was killed under a 14.1 GB limit when the host has no swap, so a 16 GB machine needs swap for it ([step 25, Memory](25-prs.md#memory)).
**Recommended:** 32 GB or more, so minimap2 and DeepVariant get close to the 32 GB they ask for.
**GRIDSS (step 4b, opt-in):** its container takes 32 GB (a 28 GB Java heap plus overhead), so it needs a machine with 32 GB or more even when nothing else runs.

---

## CPU Requirements

Each pipeline task asks for the CPUs that `conf/base.config` gives its label or its name: 1, 2, 4 or 8. The 8-CPU tasks are minimap2 (the index, and the alignment with its duplicate marking in the same task), fastp, DeepVariant, Manta and VEP. `--max_cpus` and `--max_memory` cap each task; `run-all.sh` passes the machine's CPU count (or `THREADS`) and its RAM. They do not limit how much runs at once: Nextflow starts tasks until their requests fill the machine's CPUs and RAM. On a Mac with Docker Desktop, that is the Mac's CPUs and RAM, not the Docker VM's. `MAX_JOBS` is no longer read.

**DeepVariant's shards.** `run-all.sh` runs DeepVariant with at most 8 shards, however many cores the machine has, and fewer when `--max_cpus` (`THREADS`) is lower. More cores help the steps that run beside it. To give DeepVariant more, run `scripts/03-deepvariant.sh` with `THREADS=N` (its `--cpus` and `--num_shards`), or pass `-c` with a `withName: 'DEEPVARIANT'` block:

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
- **Recommended:** 8 to 16 cores. On 8, each 8-CPU task runs alone ([what runs at once](#what-runs-at-once)); on 16, the BAM steps run beside DeepVariant. On a 16-thread server a default run from FASTQ took about 20-21 h ([runtime](#runtime-per-step)).

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

Check the shares on a running task: `docker ps` lists the pipeline's containers as `nxf-...`, and `docker inspect -f '{{.HostConfig.CpuShares}}' <container>` should print 256. The same `executor` block is the way to keep Docker Desktop's VM from being overbooked: set it to the CPUs and memory the VM has, and pass `THREADS` and `--max_memory` no larger than the VM's (e.g. `THREADS=6 ./scripts/run-all.sh <sample> <sex> --max_memory 14.GB`). On a Mac, `run-all.sh` otherwise passes the Mac's own CPU count and RAM, and Nextflow refuses a task that asks for more than the `executor` block allows.

### Runtime per step

The figures below are measured, from our end test: three real 30-40x human whole genomes, run with `run-all.sh` in October 2026 on a 16-thread, 62 GB desktop-class server with the data on an SSD, while the server ran its other services. DeepVariant ran on the CPU. Each task had its default request from `conf/base.config`; on this path a task's CPUs are a Docker share, not a cap, so a task could use more cores than it asked for when the host had them free. fastp and the alignment have two runs, from the two samples that started from FASTQ in these runs (the third sample's alignment ran on an earlier version of the pipeline that split it into two tasks, and is left out); the other steps have one run per sample, and a step that ran again on a sample counts every completed run. Each figure is the range across those runs, taken from Nextflow's trace: the task's wall time, its `peak_rss` and its average `%cpu` (in cores). `peak_rss` adds up the resident memory of the task's processes, shared and file-mapped pages included, so it can be more than the container was charged: CPSR showed up to 15.6 GB under its 8 GB limit and was not killed.

DeepVariant 1.10 with 8 shards took 8 h 50 min to 13 h 30 min. In the 13.5 h run, on an 8-CPU, 24 GB budget, `make_examples` took about 10.5 h, `call_variants` about 3 h and `postprocess_variants` under an hour; on a busy host `make_examples` ran at about half the speed it reached on a quiet one. Upstream reports 1 h 9 min on a 96-vCPU cloud machine ([DeepVariant r1.10 metrics](https://github.com/google/deepvariant/blob/r1.10/docs/metrics.md)), about 110 vCPU-hours, which agrees.

A step the end test did not run says "not measured"; where it shows a time, that is the step page's estimate. Step pages link here instead of giving their own figure.

| Step | Wall time | Peak memory | CPU (cores) | Notes |
|---|---|---|---|---|
| 1 ORA to FASTQ | not measured (estimate ~30 min) | | | only for Illumina ORA input |
| 1b fastp | 15-20 min | 1.4 GB | 7.1-7.7 | |
| 2 minimap2 index | 1-3 min | 11.4 GB | 1.0-1.6 | once per reference |
| 2 Alignment (minimap2, fixmate, sort, markdup) | 3 h 51 min - 5 h 21 min | 18.3-18.4 GB | 7.4-8.1 | one task |
| 3 DeepVariant | 8 h 50 min - 13 h 30 min | 15.9-16.5 GB | 7.3-8.9 | 8 shards; 1 h 9 min on 96 vCPUs upstream |
| 4 Manta | 22-41 min | 1.0-1.2 GB | 7.8-7.9 | |
| 5 AnnotSV | 1-1.5 min | 1.3-1.7 GB | 0.8-1.0 | |
| 6 ClinVar screen | 30 s - 1 min | 42 MB | 0.8-1.0 | |
| 7 PharmCAT | 30 s - 1 min | 0.6-0.8 GB | | the VCF preprocessor 25-48 s, then PharmCAT 6-16 s |
| 8 HLA typing (T1K) | 35 min - 1 h 5 min | 2.8-3.0 GB | 1.1-1.6 | plus 12-26 s once to build the index |
| 9 ExpansionHunter | 2-6 s | 0.2 GB | 1.1-2.2 | |
| 9b Stranger | about 1 s | | | |
| 10 TelomereHunter | 1 h 11 min - 1 h 32 min | 0.1 GB | 1.0 | |
| 11 ROH | 10-18 s | 53 MB | 0.9-1.0 | |
| 12 Mito haplogroup | 4-22 s | 0.2 GB | | the chrM extraction, haplogrep3 and haplocheck together |
| 13 VEP | 4 h 3 min - 4 h 37 min | 2.6-3.0 GB | 5.7-6.1 | PASS records only, 8 forks |
| 14 Imputation prep (opt-in) | not measured (estimate ~10 min) | | | |
| 15 duphold | 3-4 min | 1.0 GB | 3.2-4.1 | plus about a second for its filter |
| 16 indexcov | 1-3 s | | | |
| 16b mosdepth | 2-3 min | 1.9 GB | 2.7-3.5 | |
| 17 CPSR | 17-36 min | 12.3-15.6 GB (`peak_rss`, see above) | 1.8-2.3 | |
| 18 CNVpytor | 14-20 min | 1.4-2.7 GB | 2.3-3.7 | |
| 19 Delly | 21 min - 3 h 19 min | 1.2-3.1 GB | 1.0-1.1 | five runs: 21 and 25 min on two samples, 1 h 35 min to 3 h 19 min in three runs on the third |
| 20 Mito variants (Mutect2) | 1-2 min | 1.3-2.5 GB | 2.7-3.2 | |
| 21 Cyrius (opt-in) | not measured (estimate ~5-15 min) | | | `setup.sh --cyrius` installs it once |
| 22 SV consensus | 2-10 s | | | SURVIVOR's preparation, merge and sort together |
| 23 Clinical filter | 2.5-6.5 min | under 50 MB | 0.9-1.2 | |
| 24 HTML report | not measured (estimate ~1-3 min) | | | `run-all.sh` runs it after the pipeline, outside the trace |
| 25 PRS | 54 min - 1 h 17 min per sample, plus 6-11 min once per run | see notes | | genotyping the score positions 21-31 min (1.1 GB), then pgsc_calc 33-51 min with 24 GB granted ([its tasks' memory](25-prs.md#memory)); preparing the scoring files 6-11 min, 11.3-12 GB, once per run |
| 26 Ancestry (opt-in) | not measured on a genome | | | [step 25's panel numbers](25-prs.md#the-reference-panel-on-a-github-hosted-runner) were measured on a test target |
| 26b Fine-grained ancestry (opt-in) | not measured on a genome | | | measured in CI on a reduced panel: 1 min 19 s for the fractions and one painted chromosome ([step 26b](26b-ancestry-fine.md#runtime)) |
| 27 CPIC lookup | 3-5 s | | | |
| 28 MultiQC | about 3 s | 0.1-0.2 GB | | |
| 29 Somatic (opt-in) | not measured (estimate ~2-6 h) | | | under 5 min on a few genes with `INTERVALS` (estimate) |
| 30 vcfanno | not measured (estimate ~5-15 min, plus 2-3 min with CADD) | | | not in the end test |
| 31 slivar | 4-10 min | under 50 MB | 0.9-1.2 | almost all of it the prioritising pass |
| 32 pypgx | 5.5-10 min | 3.9-4.1 GB | 2.0-2.4 | plus 5-10 s for the CYP2D6 depth check |
| 36 PGx consensus | 1-7 s | | | |
| 3d Octopus (alternative) | not measured (estimate ~2-4 h) | | | |
| 4b GRIDSS (opt-in) | not measured (estimate ~4-8 h) | | | |

| Run | Runtime |
|---|---|
| Default `run-all.sh` from FASTQ, 16 threads | 20 h 14 min and 20 h 52 min, measured (two runs; one stopped at PRS and finished it with a 39-minute `-resume`) |
| Default `run-all.sh` from a BAM, 16 threads | not measured end to end; the two FASTQ runs less their fastp and alignment come to about 14.5-17 h |
| `-resume` that runs pgsc_calc again, and what follows it | 39-52 min, measured (two runs) |
| DeepVariant alone, 8 shards | 8 h 50 min - 13 h 30 min, measured (three runs) |
| DeepVariant alone, 96 vCPUs | 1 h 9 min (upstream r1.10 metrics) |
| DeepVariant alone, 4 CPUs | likely more than a day (extrapolated, not measured) |
| Default `run-all.sh`, 8 CPUs | not measured; plan for more than a day, because each 8-CPU task then runs alone |

### What runs at once

From FASTQ, the pipeline trims, aligns and marks duplicates first. indexcov then checks the sex from the BAM index, and from that moment DeepVariant and every step that reads only the BAM are ready together. The steps that read the VCF, pypgx included, wait for DeepVariant.

```
FASTQ ──> fastp ──> minimap2 + markdup (one task)  ──> BAM ──> indexcov (sex check)
                                                                   │
        ┌──────────────────────────────────────────────────────────┤
        ▼                                                          ▼
DeepVariant (8 CPUs)                      BAM steps: Manta, Delly, CNVpytor, ExpansionHunter,
        │                                 TelomereHunter, HLA, mito variants, mosdepth
        ▼                                                          │
VCF steps: ClinVar, PharmCAT, pypgx (also                          ▼
reads the BAM), VEP, CPSR, ROH, PRS,      SV chain: duphold, AnnotSV, SURVIVOR merge
vcfanno, slivar
        │                                                          │
        └─────────────────> HTML report, MultiQC <─────────────────┘
```

Nextflow starts a task only when its CPUs and memory fit in what is left of the machine's. On a host with 8 CPUs, every 8-CPU task (fastp, the alignment with its duplicate marking, DeepVariant, Manta, VEP) runs alone, and the light SV chain steps wait for it: about 18 h behind DeepVariant in the observed run, after which they finished in seconds. On 16 CPUs, DeepVariant leaves 8 free for the BAM steps.

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

Multiply the hourly price by the run time: plan for more than a day on 8 vCPUs (not measured). On 16 vCPUs `run-all.sh` still gives DeepVariant 8 shards; a default run from FASTQ took about 20-21 h on a 16-thread server ([runtime](#runtime-per-step)), and a cloud vCPU may be slower or faster than that server's threads. Add ~$0.10/GB/month for persistent disk storage. A 500 GB disk costs ~$50/month.

> **Tip:** Spot or preemptible instances cost less, but a preemption stops the run. `-resume` reuses every task that finished and starts the running one over, and DeepVariant is a single task of 8 h 50 min to 13 h 30 min with 8 shards (three measured runs).
