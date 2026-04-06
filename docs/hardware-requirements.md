# Hardware and Storage Requirements

Everything you need to know about disk space, RAM, CPU, and runtime before starting.

## TL;DR

- **1 sample:** 500 GB free disk, 16 GB RAM, 4+ CPU cores
- **2 samples:** 1 TB free disk, 32 GB RAM, 8+ CPU cores (recommended)
- **First-time setup downloads:** ~70-75 GB (reference genome, databases, Docker images)
- **Total time per sample:** 6-12 hours on a 16-core desktop

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
| CNVnator ROOT file + calls | 5-15 GB | Step 18 | ROOT file can be deleted |
| Delly BCF + VCF | 5-20 MB | Step 19 | Keep VCF, delete BCF |
| Mito analysis output | 50-200 MB | Step 20 | Keep |
| **Subtotal per sample** | **150-250 GB** | | |

### Shared Reference Data (One-Time)

| Resource | Download Size | Extracted Size | Notes |
|---|---|---|---|
| GRCh38 FASTA + index | ~3.1 GB | ~3.5 GB | Core reference genome (not compressed) |
| ClinVar VCF + index | ~200 MB | ~200 MB | Updated monthly |
| VEP cache (Ensembl 112) | ~26 GB | ~30 GB | Largest single download |
| PCGR/CPSR data bundle + VEP 113 cache | ~31 GB | ~35 GB | ClinVar + gnomAD + panels |
| Docker images (all steps) | ~10-15 GB | ~10-15 GB | Cached by Docker |
| **Subtotal (shared)** | **~73 GB** | **~85 GB** | |

### Total Disk Requirements

| Scenario | Minimum Free Space |
|---|---|
| 1 sample, core steps only (2-3-6-7) | 200 GB |
| 1 sample, full pipeline | 500 GB |
| 2 samples, full pipeline | 1 TB |
| 2 samples + keeping intermediates | 1.5 TB |

> **Tip:** After the pipeline completes, the single largest file is the BAM (80-120 GB per sample). If you're done with all BAM-dependent steps (4, 9, 10, 15, 16, 18, 19, 20), you can convert to CRAM to save 40-60% space, or delete the BAM entirely if you keep the FASTQ (you can always re-align).

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
| 18 (CNVnator) | 8 GB | 4-6 GB | ROOT file can be large |
| 19 (Delly) | 8 GB | 4-6 GB | Moderate |

**Minimum system RAM:** 16 GB (run one step at a time with reduced `--memory` flags)
**Recommended:** 32 GB (run multiple steps in parallel)
**Ideal:** 64 GB (run everything in parallel)

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
| 18 (CNVnator) | 4 | Limited | Mostly single-threaded |
| 19 (Delly) | 4 | Limited | Per-chromosome parallelism |

**Minimum:** 4 cores (very slow but works)
**Recommended:** 16 cores (good balance of speed and availability)
**No benefit beyond:** ~32 cores for any single step

### Runtime Estimates

On a 16-core / 32 GB desktop (e.g., AMD Ryzen 9 5950X):

| Step Group | Steps | Runtime | Can Parallelize? |
|---|---|---|---|
| Alignment | 2 | 1-2 hours | No (one BAM per sample) |
| Variant Calling | 3 | 2-4 hours | No (needs BAM from step 2) |
| Quick Analyses | 4, 5, 6, 7, 9, 11, 12, 16 | ~1 hour total | Yes (all independent after step 3) |
| Heavy Annotation | 13, 17 | 2-5 hours total | Yes (both use VCF) |
| Optional SV Callers | 18, 19 | 2-4 hours each | Yes (both use BAM) |
| Optional Mito/Telomere | 10, 20 | 1-2 hours total | Yes (both use BAM) |
| **Total (sequential)** | All 20 | **12-20 hours** | |
| **Total (parallelized)** | All 20 | **6-10 hours** | |

### Parallelization Strategy

After step 3 (variant calling) completes, many steps can run simultaneously:

```
Step 3 done ──┬──> Steps 4, 6, 7, 9, 11, 12, 16 (quick, ~1 hr total)
              ├──> Step 13 (VEP, ~2-4 hr)
              ├──> Step 17 (CPSR, ~30-60 min)
              ├──> Step 18 (CNVnator, ~2-4 hr)    ← These 3 use BAM, need RAM
              ├──> Step 19 (Delly, ~2-4 hr)        ← Run 1-2 at a time
              ├──> Step 10 (TelomereHunter, ~1 hr)
              └──> Step 20 (GATK Mutect2 mito, ~15-30 min)
```

---

## Internet Bandwidth

### One-Time Downloads

| Resource | Size | Notes |
|---|---|---|
| GRCh38 reference | ~3.1 GB | Not compressed, fast download |
| ClinVar | ~200 MB | Fast download |
| VEP cache | ~26 GB | Slow servers, `wget -c` recommended for resume |
| PCGR + VEP 113 cache | ~31 GB | Can be slow |
| Docker images | ~10-15 GB | Pulled automatically by `docker run` |
| **Total** | **~70-75 GB** | |

### Ongoing Downloads

- **ClinVar updates:** ~200 MB/month (optional but recommended for latest pathogenic variant classifications)
- **Docker image updates:** Variable (only when you want newer tool versions)

> **Offline operation:** After the initial setup, the core pipeline runs offline. A few steps fetch small public resources on first use if not already present: step 4b downloads the ENCODE blacklist (~50 KB), and steps 25/26 download scoring files and reference panels from public FTP servers. All downloads are cached after the first run. No sample data is ever uploaded.

---

## Storage Tips

### Save Disk Space

1. **Convert BAM to CRAM** after all BAM-dependent steps complete:
   ```bash
   samtools view -C -T reference.fasta input.bam > output.cram
   ```
   Saves 40-60% (30-50 GB per sample).

2. **Delete intermediate files:**
   - CNVnator `.root` files (5-15 GB each)
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
