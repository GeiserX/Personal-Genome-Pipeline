# Hardware and Storage Requirements

Everything you need to know about disk space, RAM, CPU, and runtime before starting. This page is the one place for the download sizes, the per-step runtimes and the totals; the other pages link here instead of repeating them.

## TL;DR

- **1 sample:** 500 GB free disk, 16 GB RAM (steps run a few at a time, see [RAM](#ram-requirements)), 4+ CPU cores
- **2 samples:** 1 TB free disk, 32 GB RAM, 8+ CPU cores (recommended)
- **First-time setup downloads:** ~75 GB for the default run, ~250 GB with the optional annotation databases ([table](#shared-reference-data-one-time))
- **Total time per sample:** 6-12 hours for a default `run-all.sh` on a 16-core desktop ([per step](#runtime-per-step))
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
| VEP annotated VCF | 2-5 GB | Step 13 | Keep (comprehensive annotation) |
| CPSR report + data | 50-200 MB | Step 17 | Keep |
| CNVpytor .pytor file + calls | 5-15 GB | Step 18 | .pytor file can be deleted |
| Delly BCF + VCF | 5-20 MB | Step 19 | Keep VCF, delete BCF |
| Mito analysis output | 50-200 MB | Step 20 | Keep |
| **Subtotal per sample** | **150-250 GB** | | |

### Shared Reference Data (One-Time)

One row per download, in GB as `wget` and `du -h` count them (1 GB = 2^30 bytes); the download sizes were read from each server on 2026-10-02. Each total is the sum of the rows above it. [Reference setup](00-reference-setup.md) has the commands, under a heading per database with the same size.

**Default run** (`setup.sh` downloads the first three rows and the Docker images; the VEP caches and the PCGR bundle are for steps 13 and 17, which a default run skips when they are missing):

| Resource | Download | On disk | Used by |
|---|---|---|---|
| GRCh38 FASTA + `.fai` | ~3 GB | ~3 GB | every step |
| ClinVar VCF + chr-renamed and pathogenic-only copies | ~0.2 GB | ~0.4 GB | step 6 |
| AnnotSV annotations | ~5 GB | ~20 GB | step 5 |
| VEP cache, release 116 | ~26 GB | ~30 GB | step 13 |
| PCGR/CPSR ref data bundle (20250314) | ~5 GB | ~5 GB | step 17 |
| VEP cache, release 113 | ~23 GB | ~27 GB | step 17 (CPSR's own VEP) |
| Docker images | ~10-15 GB | ~10-15 GB | every step |
| **Total, default run** | **~72-77 GB** | **~95-100 GB** | |

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
| Default run | ~72-77 GB |
| Default run + annotation databases | ~250 GB |

### Total Disk Requirements

| Scenario | Minimum Free Space |
|---|---|
| 1 sample, core steps only (2-3-6-7) | 200 GB |
| 1 sample, full pipeline | 500 GB |
| 2 samples, full pipeline | 1 TB |
| 2 samples + keeping intermediates | 1.5 TB |

> **Tip:** After the pipeline completes, the single largest file is the BAM (80-120 GB per sample, the figure every page uses). If you're done with all BAM-dependent steps (4, 9, 10, 15, 16, 18, 19, 20), you can convert to CRAM to save 40-60% space, or delete the BAM entirely if you keep the FASTQ (you can always re-align).

---

## RAM Requirements

Each pipeline step runs in a Docker container with a `--memory` limit. Here's what each step actually needs:

| Step | Memory Limit | Peak Usage | Notes |
|---|---|---|---|
| 2 (minimap2 alignment) | 16 GB | 6-10 GB | minimap2 is RAM-efficient |
| 3 (DeepVariant) | 32 GB | 8-20 GB | Scales with `--cpus` |
| 4 (Manta) | 8 GB | 4-6 GB | Moderate |
| 6 (ClinVar screen) | 4 GB | 1-2 GB | Light |
| 7 (PharmCAT) | 4 GB | 2-3 GB | Light |
| 9 (ExpansionHunter) | 8 GB | 4-6 GB | Moderate |
| 10 (TelomereHunter) | 8 GB | 4-6 GB | Moderate |
| 13 (VEP) | 16 GB | 4-8 GB | Cache loaded into memory |
| 17 (CPSR) | 8 GB | 4-6 GB | Moderate |
| 18 (CNVpytor) | 8 GB | 4-6 GB | .pytor (HDF5) file can be large |
| 19 (Delly) | 8 GB | 4-6 GB | Moderate |

**Minimum system RAM:** 16 GB. Every default step fits in it except possibly DeepVariant, whose peak at its 8 shards can pass 16 GB (its container may use up to 32 GB); on a 16 GB machine it can be killed for lack of memory, and lowering `--num_shards` in `scripts/03-deepvariant.sh` lowers the peak. Also, `run-all.sh` starts several containers at once (up to `MAX_JOBS`, half the CPU count with a minimum of 4) and counts CPUs, not memory. On a 16 GB machine run it with `MAX_JOBS=2`.
**Recommended:** 32 GB (run multiple steps in parallel)
**Ideal:** 64 GB (run everything in parallel)
**GRIDSS (step 4b, opt-in):** its container takes 32 GB (a 28 GB Java heap plus overhead), so it needs a machine with 32 GB or more even when nothing else runs.

> **Reducing memory limits:** If you have less RAM, edit the `--memory` flag in each script. Most steps will work with less -- they'll just be slower or may fail on edge cases. DeepVariant is the most memory-hungry.

---

## CPU Requirements

All scripts use `--cpus` to limit Docker container CPU usage. More cores = faster, but with diminishing returns above 16 cores for most tools.

| Step | Default --cpus | Scales Linearly? | Notes |
|---|---|---|---|
| 2 (minimap2) | 8 | Yes, up to ~16 | I/O bound above 16 cores |
| 3 (DeepVariant) | 8 | Yes, up to ~32 | Most CPU-intensive step |
| 4 (Manta) | 8 | Yes | Already very fast |
| 13 (VEP) | 8 | Yes (--fork) | Can use all available cores |
| 18 (CNVpytor) | 4 | Yes (-j flag) | Multi-threaded via -j |
| 19 (Delly) | 4 | Limited | Per-chromosome parallelism |

**Minimum:** 4 cores (very slow but works)
**Recommended:** 16 cores (good balance of speed and availability)
**No benefit beyond:** ~32 cores for any single step

### Runtime per step

On a 16-core / 32 GB desktop, with each script's default CPU limit. These are estimates from the step pages, except the rows marked measured. Step pages link here instead of giving their own figure.

| Step | Runtime | Notes |
|---|---|---|
| 1 ORA to FASTQ | ~30 min | only for Illumina ORA input |
| 1b fastp | ~10-20 min | |
| 2 Alignment (minimap2) | ~1-2 h | plus ~30 min once to build the `.mmi` index |
| 3 DeepVariant | ~3-5 h | measured with DeepVariant 1.6.0 on 8 threads ([benchmarking](benchmarking.md#runtime-full-genome-30x-wgs)) |
| 4 Manta | ~20 min | |
| 5 AnnotSV | ~10 min | |
| 6 ClinVar screen | ~5 min | |
| 7 PharmCAT | ~10 min | |
| 8 HLA typing (T1K) | ~30 min | plus ~35 min once to build the index |
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
| 21 Cyrius | ~5-15 min | includes the pip install |
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
| Minimum useful run (steps 2, 3, 6, 7) | ~4-7 h |
| Default `run-all.sh`, steps in parallel | ~6-12 h |

### Parallelization Strategy

After step 3 (variant calling) completes, many steps can run simultaneously:

```
Step 3 done ──┬──> Steps 4, 6, 7, 9, 11, 12, 16 (quick, ~1 hr total)
              ├──> Step 13 (VEP, ~2-4 hr) ──> Step 30 (vcfanno, ~15 min) ──> Step 31 (slivar)
              ├──> Step 17 (CPSR, ~30-60 min)
              ├──> Step 18 (CNVpytor, ~1-3 hr)    ← These 3 use BAM, need RAM
              ├──> Step 19 (Delly, ~2-4 hr)        ← Run 1-2 at a time
              ├──> Step 32 (pypgx, ~20-40 min)     ← Uses BAM, parallel with above
              ├──> Step 10 (TelomereHunter, ~1 hr)
              └──> Step 20 (GATK Mutect2 mito, ~15-30 min)
```

---

## Internet Bandwidth

### One-Time Downloads

About 75 GB for a default run and 250 GB with the annotation databases; the [table above](#shared-reference-data-one-time) has every download. The VEP caches come from slow servers, so use `wget -c` to resume.

### Ongoing Downloads

- **ClinVar updates:** ~200 MB/month (optional but recommended for latest pathogenic variant classifications)
- **Docker image updates:** Variable (only when you want newer tool versions)

> **Network during a run:** after setup, a few steps still fetch public files (the HLA database, PGS scoring files, Cyrius from PyPI, MultiQC's update check). None sends sample data. [Why run locally?](why-local.md#network-calls-during-a-run) lists each one.

---

## Storage Tips

### Save Disk Space

1. **Convert BAM to CRAM** after all BAM-dependent steps complete:
   ```bash
   samtools view -C -T reference.fasta input.bam > output.cram
   ```
   Saves 40-60% (30-50 GB per sample).

2. **Delete intermediate files:**
   - CNVpytor `.pytor` files (5-15 GB each)
   - Delly `.bcf` files (after converting to VCF)

3. **Compress VEP output:**
   ```bash
   bgzip sample_vep.vcf  # Compresses from ~3.5 GB to ~400 MB
   ```

4. **Delete FASTQ** if you have the BAM and don't plan to re-align. You can always re-extract FASTQ from BAM if needed.

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

| Provider | Instance | vCPUs | RAM | Cost/hr | ~Cost per Sample |
|---|---|---|---|---|---|
| AWS | c5.4xlarge | 16 | 32 GB | ~$0.68 | ~$5-8 |
| GCP | n2-standard-16 | 16 | 64 GB | ~$0.78 | ~$6-10 |
| Azure | Standard_D16s_v5 | 16 | 64 GB | ~$0.77 | ~$6-10 |
| Hetzner | CCX33 | 8 | 32 GB | ~$0.18 | ~$2-3 |

Add ~$0.10/GB/month for persistent disk storage. A 500 GB disk costs ~$50/month.

> **Tip:** Use spot/preemptible instances for 60-80% savings. The pipeline is restartable -- if your instance gets preempted, just re-run the interrupted step.
