<p align="center">
  <img src="docs/images/banner.svg" alt="Genomics Pipeline banner" width="900"/>
</p>

<h1 align="center">Genomics Pipeline</h1>

<p align="center">
  <strong>Analyze your own whole genome sequencing (WGS) data on consumer hardware.</strong><br>
  No cloud accounts, no subscriptions, no bioinformatics degree required.
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/GeiserX/genomics-pipeline?style=flat-square" alt="License"></a>
  <a href="https://github.com/GeiserX/genomics-pipeline/actions/workflows/lint.yml"><img src="https://img.shields.io/github/actions/workflow/status/GeiserX/genomics-pipeline/lint.yml?style=flat-square&label=CI" alt="CI"></a>
  <a href="https://github.com/GeiserX/genomics-pipeline/stargazers"><img src="https://img.shields.io/github/stars/GeiserX/genomics-pipeline?style=flat-square&logo=github" alt="GitHub Stars"></a>
  <a href="https://www.docker.com/"><img src="https://img.shields.io/badge/runs%20with-Docker-0db7ed?style=flat-square&logo=docker&logoColor=white" alt="Docker"></a>
  <a href="https://geiserx.github.io/genomics-pipeline"><img src="https://img.shields.io/badge/docs-GitHub%20Pages-blue?style=flat-square&logo=github" alt="Docs"></a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/genome-GRCh38%2Fhg38-22c55e?style=flat-square" alt="GRCh38">
  <img src="https://img.shields.io/badge/analysis%20steps-27-f97316?style=flat-square" alt="27 steps">
  <img src="https://img.shields.io/badge/data%20privacy-100%25%20local-06b6d4?style=flat-square&logo=data:image/svg%2bxml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiPjxwYXRoIGQ9Ik0xMiAxYTUgNSAwIDAgMC01IDV2Mkg1djE0aDE0VjhIMTdWNmE1IDUgMCAwIDAtNS01em0tMyA1YTMgMyAwIDEgMSA2IDB2MkgzVjZ6Ii8+PC9zdmc+" alt="100% Local">
  <img src="https://img.shields.io/badge/platform-linux%20%7C%20macOS%20%7C%20WSL2-lightgrey?style=flat-square" alt="Platform">
</p>

This pipeline takes raw sequencing data (FASTQ/BAM/VCF) from any vendor and runs 27 analysis steps to produce a comprehensive genomic profile: variant calling, pharmacogenomics, structural variants, cancer predisposition screening, polygenic risk scores, ancestry estimation, telomere length, mitochondrial analysis, and more. Everything runs locally in Docker containers with resource limits so it won't crash your machine.

**Time:** 6-12 hours per sample on a 16-core desktop | **Disk:** 500 GB minimum per sample | **Cost:** Free (you just need your data)

---

## Who Is This For?

- You got WGS from **Nebula/DNA Complete, Dante Labs, Sequencing.com, Novogene**, or any other vendor and want to analyze it yourself
- You have clinical WGS data (Illumina DRAGEN, BAM+VCF from a hospital) and want deeper analysis than the lab report
- You're a biohacker, researcher, or patient advocate who wants full control over your genomic data
- You can't afford a $500/hour genetics consultant but you have a computer and curiosity

> **Not WGS?** If you have genotyping array data (23andMe, AncestryDNA, MyHeritage), this pipeline won't work directly. Those services test ~600K positions, not the full 3 billion. You can convert array data to VCF and use some downstream steps (ClinVar screening, pharmacogenomics), but alignment and variant calling require actual sequencing reads. See the [vendor guide](docs/vendor-guide.md) for details.

---

## What You Get

| Category | What It Finds | Steps |
|---|---|---|
| **Variant Calling** | SNPs, indels, structural variants, copy number variants | 3, 4, 18, 19 |
| **Clinical Screening** | Pathogenic variants, carrier status, cancer predisposition (ACMG SF v3.2) | 6, 17 |
| **Pharmacogenomics** | Drug-gene interactions (21+ genes, CYP2C19, CYP2D6, DPYD, etc.) | 7, 21, 27 |
| **Structural Variants** | Deletions, duplications, inversions, translocations (3 callers + consensus) | 4, 5, 15, 18, 19, 22 |
| **Functional Annotation** | Impact prediction for every variant (VEP: SIFT, PolyPhen, gnomAD frequencies) | 13 |
| **Repeat Expansions** | Huntington's, Fragile X, ALS, and 50+ other repeat expansion disorders | 9 |
| **Ancestry & Haplogroups** | Mitochondrial haplogroup, consanguinity check, ancestry PCA | 11, 12, 26 |
| **Telomere Length** | Biological age proxy from telomere content | 10 |
| **Mitochondrial** | Heteroplasmy detection, mitochondrial disease variants | 12, 20 |
| **Polygenic Risk** | Risk scores for 10 common conditions (CAD, T2D, cancers, etc.) | 25 |
| **Quality Control** | Coverage uniformity, sex chromosome verification, SV false positive filtering | 15, 16 |

---

## Pipeline Overview

```
FASTQ/BAM ──> Alignment ──> Sorted BAM ──┬──> DeepVariant (SNPs/indels) ──> VCF
              (minimap2)                  │        │
                                          │        ├──> ClinVar Screen
                                          │        ├──> PharmCAT (PGx) ──> CPIC Recommendations
                                          │        ├──> VEP Annotation ──> Clinical Filter
                                          │        ├──> CPSR (Cancer Predisposition)
                                          │        ├──> ExpansionHunter (STRs)
                                          │        ├──> ROH Analysis
                                          │        ├──> PRS (Polygenic Risk Scores)
                                          │        ├──> Ancestry PCA
                                          │        └──> Imputation Prep
                                          │
                                          ├──> Manta (SVs) ──┐
                                          ├──> Delly (SVs) ──┼──> SV Consensus ──> duphold ──> AnnotSV
                                          ├──> CNVnator (CNVs)┘
                                          ├──> Cyrius (CYP2D6)
                                          ├──> TelomereHunter
                                          ├──> indexcov (Coverage QC)
                                          ├──> MToolBox (Mitochondrial)
                                          └──> Haplogrep3 (mtDNA haplogroup)
                                                              └──> HTML Report (aggregates all)
```

### All Steps

| # | Step | Tool | Docker Image | Runtime | Required? |
|---|---|---|---|---|---|
| 1 | [ORA to FASTQ](docs/01-ora-to-fastq.md) | orad | `orad` binary | ~30 min | Only for Illumina ORA files |
| 2 | [Alignment](docs/02-alignment.md) | minimap2 + samtools | `staphb/samtools:1.20` | ~1-2 hr | Yes (if starting from FASTQ) |
| 3 | [Variant Calling](docs/03-variant-calling.md) | DeepVariant | `google/deepvariant:1.6.0` | ~2-4 hr | Yes |
| 4 | [Structural Variants](docs/04-structural-variants.md) | Manta | `quay.io/biocontainers/manta` | ~20 min | Recommended |
| 5 | [SV Annotation](docs/05-annotsv.md) | AnnotSV | `getwilds/annotsv:latest` | ~10 min | If step 4 run |
| 6 | [ClinVar Screen](docs/06-clinvar-screen.md) | bcftools isec | `staphb/bcftools:1.21` | ~5 min | Yes |
| 7 | [Pharmacogenomics](docs/07-pharmacogenomics.md) | PharmCAT | `pgkb/pharmcat:2.15.5` | ~10 min | Yes |
| 8 | [HLA Typing](docs/08-hla-typing.md) | T1K | `quay.io/biocontainers/t1k:1.0.9` | ~30 min | Optional |
| 9 | [STR Expansions](docs/09-str-expansions.md) | ExpansionHunter | `weisburd/expansionhunter` | ~15 min | Recommended |
| 10 | [Telomere Length](docs/10-telomere-analysis.md) | TelomereHunter | `lgalarno/telomerehunter:latest` | ~1 hr | Optional |
| 11 | [ROH Analysis](docs/11-roh-analysis.md) | bcftools roh | `staphb/bcftools:1.21` | ~5 min | Recommended |
| 12 | [Mito Haplogroup](docs/12-mito-haplogroup.md) | haplogrep3 | `genepi/haplogrep3` | ~1 min | Optional |
| 13 | [VEP Annotation](docs/13-vep-annotation.md) | VEP | `ensemblorg/ensembl-vep:release_112.0` | ~2-4 hr | Recommended |
| 14 | [Imputation Prep](docs/14-imputation-prep.md) | bcftools | `staphb/bcftools:1.21` | ~10 min | Optional |
| 15 | [SV Quality](docs/15-duphold.md) | duphold | `brentp/duphold:latest` | ~20 min | If step 4 run |
| 16 | [Coverage QC](docs/16-indexcov.md) | indexcov | `quay.io/biocontainers/goleft` | ~5 sec | Recommended |
| 17 | [Cancer Predisposition](docs/17-cpsr.md) | CPSR | `sigven/pcgr:1.4.1` | ~30-60 min | Recommended |
| 18 | [CNV Calling](docs/18-cnvnator.md) | CNVnator | `quay.io/biocontainers/cnvnator` | ~2-4 hr | Optional |
| 19 | [SV Calling (Delly)](docs/19-delly.md) | Delly | `quay.io/biocontainers/delly:1.7.3` | ~2-4 hr | Optional |
| 20 | [Mitochondrial](docs/20-mtoolbox.md) | GATK Mutect2 | `broadinstitute/gatk:4.6.1.0` | ~15-30 min | Optional |

#### Post-Processing Steps

These run after the core pipeline completes and combine outputs from earlier steps.

| # | Step | Tool | Docker Image | Runtime | Required? |
|---|---|---|---|---|---|
| 21 | [CYP2D6 Star Alleles](docs/21-cyrius.md) | Cyrius | `python:3.11-slim` | ~10 min | Optional |
| 22 | [SV Consensus Merge](docs/22-survivor-merge.md) | bcftools | `staphb/bcftools:1.21` | ~5 min | If steps 4+19 run |
| 23 | [Clinical Filter](docs/23-clinical-filter.md) | bcftools +split-vep | `staphb/bcftools:1.21` | ~5-10 min | If step 13 run |
| 24 | [HTML Report](docs/24-html-report.md) | bash + bcftools | `staphb/bcftools:1.21` | ~1-3 min | Recommended |
| 25 | [Polygenic Risk Scores](docs/25-prs.md) | plink2 | `pgscatalog/plink2:2.00a5.10` | ~30 min | Optional |
| 26 | [Ancestry PCA](docs/26-ancestry.md) | plink2 | `pgscatalog/plink2:2.00a5.10` | ~30-60 min | Optional |
| 27 | [CPIC Recommendations](docs/27-cpic-lookup.md) | Python + CPIC | `python:3.11-slim` | ~5 min | If step 7 run |

**Minimum useful run:** Steps 2, 3, 6, 7 (alignment + variant calling + ClinVar + PharmCAT) = ~4-6 hours.
**Full analysis:** All 27 steps = ~12-20 hours. Steps 4/18/19 and 10/12/20 can run in parallel.

---

## Quick Start

### Step 0: Quick Test (Optional)

Verify everything works on a small public dataset before committing to a full run:

```bash
# See docs/quick-test.md for full instructions
# Tests VCF-only steps in under 5 minutes using GIAB NA12878 data
```

### Step 0.5: Validate Your Setup

Before running on your own data, verify that all prerequisites are in place:

```bash
export GENOME_DIR=/path/to/your/data
./scripts/validate-setup.sh your_name
```

This checks Docker, disk space, reference data, Docker images, and sample files. Fix any `[FAIL]` items before proceeding.

### Path A: I Have FASTQ Files (Raw Reads)

Most common if you downloaded data from Nebula, Dante Labs, Novogene, BGI, or any sequencing provider.

```bash
# 1. Set your data directory (where your FASTQ files are)
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name

# 2. Download the GRCh38 reference genome (~3.1 GB)
mkdir -p ${GENOME_DIR}/reference
# Download Homo_sapiens_assembly38.fasta + .fai from GATK resource bundle
# See docs/00-reference-setup.md for details

# 3. Run the pipeline
./scripts/02-alignment.sh $SAMPLE        # FASTQ -> sorted BAM (~1-2 hr)
./scripts/03-deepvariant.sh $SAMPLE      # BAM -> VCF (~2-4 hr)
./scripts/06-clinvar-screen.sh $SAMPLE   # Find pathogenic variants (~5 min)
./scripts/07-pharmacogenomics.sh $SAMPLE # Drug-gene interactions (~10 min)

# 4. Optional: structural variants, annotation, etc.
./scripts/04-manta.sh $SAMPLE
./scripts/13-vep-annotation.sh $SAMPLE
./scripts/17-cpsr.sh $SAMPLE
# ... see full step table above
```

### Path B: I Have a BAM File (Aligned Reads)

Common if your lab or vendor already aligned the reads (Illumina DRAGEN output, clinical labs).

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name

# Your BAM should be at: ${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
# Skip step 2 (alignment) and start directly with variant calling:
./scripts/03-deepvariant.sh $SAMPLE
./scripts/06-clinvar-screen.sh $SAMPLE
./scripts/07-pharmacogenomics.sh $SAMPLE
```

### Path C: I Have a VCF File (Variant Calls)

If you already have variants called (from DRAGEN, GATK, or another pipeline).

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name

# Your VCF should be at: ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz
# Skip steps 2-3 and go straight to analysis:
./scripts/06-clinvar-screen.sh $SAMPLE
./scripts/07-pharmacogenomics.sh $SAMPLE
./scripts/13-vep-annotation.sh $SAMPLE
./scripts/17-cpsr.sh $SAMPLE
```

### Path D: I Have Illumina ORA Files

ORA is Illumina's proprietary compressed FASTQ format. Decompress first, then follow Path A.

```bash
./scripts/01-ora-to-fastq.sh $SAMPLE   # ORA -> FASTQ
./scripts/02-alignment.sh $SAMPLE       # FASTQ -> BAM
# ... continue as Path A
```

---

## Prerequisites

### Hardware Requirements

| Resource | Minimum | Recommended | Notes |
|---|---|---|---|
| **CPU** | 4 cores | 16+ cores | DeepVariant scales linearly with cores |
| **RAM** | 16 GB | 32 GB | Some steps need 8-16 GB; pipeline limits each container |
| **Disk** | 500 GB free | 1 TB+ | See [detailed breakdown](docs/hardware-requirements.md) |
| **Internet** | Broadband | 100+ Mbps | ~100 GB of one-time downloads (reference genome, databases, Docker images) |
| **OS** | Linux (amd64) | Ubuntu 22.04+ | macOS/ARM works but slower (see below) |

> **Disk space is the #1 surprise.** A single 30X WGS sample produces 60-90 GB of FASTQ, 30-80 GB of BAM, plus reference genomes and databases. See [docs/hardware-requirements.md](docs/hardware-requirements.md) for the full breakdown.

### Software

| Software | Version | Install |
|---|---|---|
| Docker | 20.10+ | [docs.docker.com/get-docker](https://docs.docker.com/get-docker/) |
| bash | 4.0+ | Pre-installed on Linux/macOS |
| wget or curl | Any | For downloading references |

That's it. Every analysis tool runs inside Docker -- no conda environments, no Python version conflicts, no compilation.

### Reference Data (One-Time Downloads)

| Resource | Size | Required For |
|---|---|---|
| GRCh38 reference FASTA + index | ~3.5 GB | All steps |
| ClinVar database | ~200 MB | Step 6 (ClinVar screen) |
| VEP cache | ~26 GB | Step 13 (VEP annotation) |
| PCGR/CPSR data bundle | ~21 GB | Step 17 (cancer predisposition) |
| Docker images (all steps) | ~10-15 GB | All steps |
| **Total one-time setup** | **~60-65 GB** | |

See [docs/00-reference-setup.md](docs/00-reference-setup.md) for download instructions.

---

## Platform Notes

### Linux (Recommended)
Best performance. Docker runs natively. All pipeline images are linux/amd64. No issues.

### macOS (Intel)
Works fine. Docker Desktop runs a Linux VM, so there's a ~10-20% I/O overhead on file operations. Set Docker Desktop memory to at least 16 GB (Preferences > Resources).

### macOS (Apple Silicon / M1-M4)
Works but **slower**. All bioinformatics Docker images are amd64 and run under Rosetta 2 emulation (2-5x performance penalty). DeepVariant and BWA-MEM2 are the most affected. Set Docker Desktop to use Rosetta 2 for amd64 emulation (enabled by default on newer versions).

### Windows (WSL2)
Works. Install Docker Desktop with WSL2 backend. **Critical:** Keep all genomics data on the Linux filesystem (`~/data/`, not `/mnt/c/`). Accessing Windows drives from WSL2 is 10-50x slower due to the 9P protocol. Set WSL2 memory in `%UserProfile%\.wslconfig`:
```ini
[wsl2]
memory=24GB
swap=8GB
```

### Unraid / NAS Servers
Works great for long-running analyses. Use `--cpus` and `--memory` Docker flags (already set in all scripts) to avoid starving other services. Consider running in detached mode (`-d` flag) for multi-hour steps.

---

## Data from Your Vendor

Different vendors deliver data in different formats. Here's what you need to know:

| Vendor | Format You Get | Pipeline Entry Point | Notes |
|---|---|---|---|
| **Nebula / DNA Complete** | FASTQ + VCF | Path A (FASTQ) or Path C (VCF) | Uses BGI/MGI sequencing |
| **Dante Labs** | FASTQ + BAM + VCF | Any path | Standard Illumina |
| **Sequencing.com** | FASTQ + BAM + VCF | Any path | Standard Illumina |
| **Novogene / BGI** | FASTQ | Path A | BGI read names differ from Illumina but work fine |
| **Illumina DRAGEN (clinical)** | ORA or BAM + VCF | Path D (ORA) or Path B/C | ORA needs decompression first |
| **Oxford Nanopore** | POD5/FAST5 + BAM | Not supported | Long-read data needs different tools (Clair3, Sniffles2) |
| **PacBio HiFi** | HiFi BAM | Not supported | Needs PacBio-specific pipeline (pbmm2, DeepVariant PacBio model) |
| **23andMe / Ancestry / MyHeritage** | Genotyping array TSV | Partial (VCF steps only) | Not WGS -- convert to VCF first |

See [docs/vendor-guide.md](docs/vendor-guide.md) for detailed conversion instructions for each vendor.

---

## Directory Structure

The pipeline expects this layout (created automatically by the scripts):

```
${GENOME_DIR}/
  reference/
    Homo_sapiens_assembly38.fasta      # GRCh38 reference genome
    Homo_sapiens_assembly38.fasta.fai  # FASTA index
  clinvar/
    clinvar.vcf.gz                     # ClinVar database
    clinvar.vcf.gz.tbi                 # ClinVar index
  vep_cache/                           # VEP annotation cache (~30 GB)
  pcgr_data/                           # CPSR/PCGR data bundle (~21 GB)
  ${SAMPLE}/
    fastq/                             # Raw FASTQ files (R1 + R2)
    aligned/
      ${SAMPLE}_sorted.bam             # Aligned reads
      ${SAMPLE}_sorted.bam.bai         # BAM index
    vcf/
      ${SAMPLE}.vcf.gz                 # Variant calls
      ${SAMPLE}.vcf.gz.tbi             # VCF index
    manta/                             # Structural variants (step 4)
    annotsv/                           # Annotated SVs (step 5)
    clinvar/                           # ClinVar hits (step 6)
    pharmcat/                          # Pharmacogenomics report (step 7)
    vep/                               # Functional annotation (step 13)
    cpsr/                              # Cancer predisposition (step 17)
    ...                                # Other analysis directories
```

---

## Common Issues

| Problem | Cause | Fix |
|---|---|---|
| Container exits silently | Out of memory (OOM killed) | Increase Docker memory or reduce `--memory` flag. Check `docker logs <container>`. |
| "Permission denied" writing output | Container runs as non-root | Add `--user root` to `docker run` (already done in all scripts) |
| VEP cache download fails/times out | 26 GB download over unreliable connection | Use `wget -c` (supports resume). See [docs/13-vep-annotation.md](docs/13-vep-annotation.md) |
| DeepVariant crashes on Mac | amd64 emulation + memory pressure | Reduce `--cpus` to 2 and `--memory` to 8g. Will be slow. |
| Wrong number of variants (too few) | Genome build mismatch | Ensure your BAM is aligned to GRCh38 (hg38), not hg19/GRCh37. Check with `samtools view -H your.bam \| grep SN:chr1` |
| 0-byte output files | Missing input or wrong path | Check that all input files exist. Run the script with `bash -x` for debug output. |
| "No such image" on `docker pull` | Image name/tag changed | Check the exact image name in the step's documentation. Biocontainer tags change frequently. |
| Very slow on macOS | Rosetta 2 emulation overhead | Expected. Consider running on a Linux machine or cloud instance for heavy steps. |

For detailed solutions, see:
- [docs/troubleshooting.md](docs/troubleshooting.md) — comprehensive troubleshooting guide organized by symptom
- [docs/lessons-learned.md](docs/lessons-learned.md) — every failure encountered during development
- [docs/glossary.md](docs/glossary.md) — alphabetical glossary of genomics terms

---

## FAQ

**Q: How much does WGS cost?**
$200-$1,000 depending on the vendor. Nebula/DNA Complete: $495 for 30X. Dante Labs: ~$300-600. Sequencing.com: $399. Novogene (research): ~$200-400. The pipeline itself is free.

**Q: I only have 23andMe/AncestryDNA data. Can I use this?**
Partially. Those services use genotyping arrays (~600K positions), not full sequencing (~3 billion). You can convert the raw data to VCF and run ClinVar screening (step 6) and some annotation steps, but you can't run alignment, variant calling, or structural variant analysis. For the full pipeline, you need actual WGS data.

**Q: How long does the full pipeline take?**
On a 16-core/32GB desktop: ~6-12 hours per sample for the core steps. The full 27-step pipeline takes ~12-20 hours. Many steps can run in parallel (Manta + CNVnator + Delly, or TelomereHunter + MToolBox + haplogrep3).

**Q: Can I run this on a Raspberry Pi?**
No. Most bioinformatics Docker images are amd64 only, and a Pi doesn't have enough RAM. Minimum is a desktop/server with 16 GB RAM and an x86_64 CPU.

**Q: My data is aligned to hg19/GRCh37. What do I?**
Extract FASTQ from your BAM (`samtools fastq`) and re-align to GRCh38 using step 2. LiftOver is an alternative but introduces artifacts. Re-alignment is cleaner.

**Q: I found a pathogenic variant. Should I be worried?**
Probably not. Everyone carries 2-10 pathogenic variants, almost all heterozygous (one copy) for recessive conditions. This means you're a **carrier**, not affected. Only worry if: (1) the variant is in a **dominant** gene, (2) you have **two** pathogenic variants in the same recessive gene, or (3) it is in a cancer predisposition gene (BRCA1/2, MLH1, etc.). See [interpreting-results.md](docs/interpreting-results.md) for details.

**Q: What are VUS? Should I worry about them?**
VUS (Variants of Uncertain Significance) mean there is not enough evidence to classify the variant as pathogenic or benign. The majority will eventually be reclassified as benign. They are **not actionable** — do not change your medical care based on a VUS. Check back in 1-2 years with an updated ClinVar database.

**Q: How often should I re-run the analysis?**
ClinVar and other databases are updated monthly. Re-running the ClinVar screen (step 6, ~5 minutes) and CPSR (step 17, ~30 minutes) every 6-12 months with updated databases can catch newly classified variants. The compute-heavy steps (alignment, variant calling) do not need to be re-run unless you get new sequencing data.

**Q: I ran the pipeline on two people (me and my partner). How do I compare?**
See [docs/multi-sample.md](docs/multi-sample.md) for carrier cross-screening, pharmacogenomics comparison, and family analysis.

**Q: Is this clinically validated?**
No. This is a research/educational pipeline. It uses the same tools as clinical labs (DeepVariant, VEP, ClinVar, PharmCAT) but has not been through clinical validation. Always discuss findings with a healthcare provider.

**Q: What about long-read sequencing (Nanopore, PacBio)?**
This pipeline is designed for short-read Illumina/BGI data. Long-read data requires different alignment (minimap2 with `-ax map-ont` or `pbmm2`) and variant calling tools (Clair3 for SNPs, Sniffles2 for SVs). A long-read branch may be added in the future.

**Q: My vendor's download link expired. Can I still get my data?**
It depends on the vendor. Dante Labs deletes data after 30 days. Sequencing.com archives to cold storage (1-3 day retrieval). DNA Complete requires an active subscription. Novogene/BGI keep data for ~90 days. **Always download your raw data immediately.** See [docs/vendor-guide.md](docs/vendor-guide.md) for vendor-specific details.

---

## Why Run Locally?

### Cost Comparison

| Approach | Cost | What You Get | Data Privacy |
|---|---|---|---|
| **This pipeline** | $0 (free, open source) | Full 24-step analysis | Your data never leaves your machine |
| Clinical WGS interpretation | $500-5,000 | 1-page report, selected genes only | Lab retains your data |
| Nebula/Dante report | $0-200 (included/add-on) | Web dashboard, limited depth | Data on company servers |
| 23andMe Health | $229 | ~10 health reports from array data | Data shared with research partners |
| Genetic counselor consultation | $200-500/hour | Expert interpretation of specific findings | HIPAA-protected |

**The pipeline is complementary, not a replacement.** Use it for comprehensive self-analysis, then bring specific findings to a genetic counselor or physician for clinical interpretation.

### Privacy and Security

Your genome is the most permanent piece of personal data you have. Unlike a password, you cannot change it if it leaks.

**This pipeline is 100% local:**
- No data is uploaded to any server, cloud, or API
- No internet connection is required after initial setup (reference downloads + Docker images)
- No telemetry, no analytics, no tracking
- Every tool runs in a Docker container with no network access to external services

**Recommendations for securing your data:**
- Store genomic data on an encrypted filesystem (LUKS on Linux, FileVault on macOS, BitLocker on Windows)
- Never upload raw FASTQ/BAM/VCF files to unencrypted cloud storage
- If using a NAS, enable encryption at rest
- Be cautious with VCF files in particular — they are small enough to accidentally email or upload
- Consider who has physical access to the machine where your data is stored
- If you delete your data, use `shred` (Linux) or secure erase — standard file deletion leaves data recoverable

**GDPR note:** If you are in the EU, your genomic data is classified as "special category personal data" under GDPR Article 9. Processing it locally for personal use is lawful. Sharing it with third parties (including cloud services) may require explicit consent and appropriate safeguards.

---

## Disclaimer

This pipeline is for **educational and research purposes only**. It is not a medical device and has not been clinically validated. Genomic findings should always be discussed with a qualified healthcare professional before making any medical decisions. The authors are not responsible for any actions taken based on pipeline output.

Your genome data is sensitive personal information. This pipeline runs entirely locally -- no data is uploaded anywhere. Keep your data secure.

---

## Lessons Learned

See [docs/lessons-learned.md](docs/lessons-learned.md) for every failure encountered during development and how it was resolved. This includes Docker image issues, tool-specific bugs, bcftools quirks, VEP cache problems, and general Docker tips.

---

## Learn More

- [docs/interpreting-results.md](docs/interpreting-results.md) — what your results mean, with examples
- [docs/multi-sample.md](docs/multi-sample.md) — comparing two or more genomes (partners, family)
- [docs/resources.md](docs/resources.md) — free courses, databases, and tools for learning genomics
- [docs/quick-test.md](docs/quick-test.md) — verify your setup with public test data

## Contributing

Found a bug? Have a tool suggestion? See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. The pipeline is designed to be extended -- each step is a standalone script with its own documentation.

## License

GPL-3.0 — see [LICENSE](LICENSE)
