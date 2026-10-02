# Getting started

## Who Is This For?

- You got WGS from **Nebula/DNA Complete, Dante Labs, Sequencing.com, Novogene**, or any other vendor and want to analyze it yourself
- You have clinical WGS data (Illumina DRAGEN, BAM+VCF from a hospital) and want deeper analysis than the lab report
- You're a biohacker, researcher, or patient advocate who wants full control over your genomic data
- You can't afford a $500/hour genetics consultant but you have a computer and curiosity

> **Only have 23andMe / MyHeritage / AncestryDNA?** You can still run pharmacogenomics, polygenic risk scores, ClinVar screening, and ROH analysis. See the **[chip data guide](chip-data-guide.md)** for step-by-step conversion instructions and which pipeline steps work with ~600K SNP array data.

## Prerequisites

### Hardware Requirements

| Resource | Minimum | Recommended | Notes |
|---|---|---|---|
| **CPU** | 4 cores | 16+ cores | DeepVariant scales linearly with cores |
| **RAM** | 16 GB | 32 GB | Some steps need 8-16 GB; pipeline limits each container |
| **Disk** | 500 GB free | 1 TB+ | See [detailed breakdown](hardware-requirements.md) |
| **Internet** | Broadband | 100+ Mbps | ~70-75 GB core downloads + ~175 GB optional annotation databases |
| **OS** | Linux (amd64) | Ubuntu 22.04+ | macOS/ARM works but slower (see below) |

> **Disk space is the #1 surprise.** A single 30X WGS sample produces 60-90 GB of FASTQ, 80-120 GB of BAM, plus reference genomes and databases. See [docs/hardware-requirements.md](hardware-requirements.md) for the full breakdown.

### Software

| Software | Version | Install |
|---|---|---|
| Docker | 20.10+ | [docs.docker.com/get-docker](https://docs.docker.com/get-docker/) |
| bash | 4.0+ | Pre-installed on Linux; macOS ships 3.2 — install via `brew install bash` |
| wget or curl | Any | For downloading references |
| python3 *(optional)* | 3.6+ | Used by long-read alignment (02b) for symlink resolution. Falls back to `readlink -f` on GNU/Linux if absent |

That's it. Every analysis tool runs inside Docker -- no conda environments, no Python version conflicts, no compilation.

### Reference Data (One-Time Downloads)

| Resource | Size | Required For |
|---|---|---|
| GRCh38 reference FASTA + index | ~3.5 GB | All steps |
| ClinVar database | ~200 MB | Step 6 (ClinVar screen) |
| VEP cache | ~26 GB | Step 13 (VEP annotation) |
| PCGR/CPSR data bundle + VEP 113 cache | ~31 GB | Step 17 (cancer predisposition) |
| Docker images (all steps) | ~10-15 GB | All steps |
| Annotation databases (CADD, SpliceAI, REVEL, AlphaMissense) | ~175 GB | Steps 30-31 (optional) |
| **Total one-time setup (core)** | **~70-75 GB** | |
| **Total with annotation enrichment** | **~250 GB** | |

See [docs/00-reference-setup.md](00-reference-setup.md) for download instructions.

## Platform Notes

### Linux (Recommended)
Best performance. Docker runs natively. All pipeline images are linux/amd64. No issues.

### macOS (Intel)
Works fine. Docker Desktop runs a Linux VM, so there's a ~10-20% I/O overhead on file operations. Set Docker Desktop memory to at least 16 GB (Preferences > Resources).

### macOS (Apple Silicon / M1-M4)
Works but **slower**. All bioinformatics Docker images are amd64 and run under Rosetta 2 emulation (2-5x performance penalty). DeepVariant and BWA-MEM2 are the most affected. In Docker Desktop, turn on "Use Rosetta for x86_64/amd64 emulation on Apple Silicon" if your version offers it; it needs the Apple Virtualization framework as the virtual machine manager, not Docker VMM. Without it, Docker falls back to QEMU emulation, which is slower still.

### Windows (WSL2)
Works. Install Docker Desktop with WSL2 backend. **Critical:** Keep all genomics data on the Linux filesystem (`~/data/`, not `/mnt/c/`). Accessing Windows drives from WSL2 is 10-50x slower due to the 9P protocol. Set WSL2 memory in `%UserProfile%\.wslconfig`:
```ini
[wsl2]
memory=24GB
swap=8GB
```

### Unraid / NAS Servers
Works great for long-running analyses. Use `--cpus` and `--memory` Docker flags (already set in all scripts) to avoid starving other services. Consider running in detached mode (`-d` flag) for multi-hour steps.

## Quick Start

### Step 0: Quick Test (Optional)

Verify everything works on a small public dataset before committing to a full run:

```bash
# See docs/quick-test.md for full instructions
# VCF-only steps finish in under 5 minutes — see docs/quick-test.md for full guide (~30 min including downloads)
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
./scripts/01b-fastp-qc.sh $SAMPLE        # QC + adapter trimming (~10-20 min)
./scripts/02-alignment.sh $SAMPLE        # FASTQ -> sorted BAM (~1-2 hr)
./scripts/03-deepvariant.sh $SAMPLE      # BAM -> VCF (~3-5 hr)
./scripts/06-clinvar-screen.sh $SAMPLE   # Find pathogenic variants (~5 min)
./scripts/07-pharmacogenomics.sh $SAMPLE # Drug-gene interactions (~10 min)

# 4. Optional: structural variants, annotation, etc.
./scripts/04-manta.sh $SAMPLE
./scripts/13-vep-annotation.sh $SAMPLE
./scripts/17-cpsr.sh $SAMPLE
# ... see the full step table in docs/pipeline-overview.md
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
# ORA -> FASTQ, one call per ORA file: <sample> <ora_reference_dir> <ora_file>
./scripts/01-ora-to-fastq.sh $SAMPLE /path/to/oradata /path/to/${SAMPLE}_S1_L001_R1_001.fastq.ora
./scripts/01-ora-to-fastq.sh $SAMPLE /path/to/oradata /path/to/${SAMPLE}_S1_L001_R2_001.fastq.ora
# orad keeps the ORA file name; the next steps read ${SAMPLE}_R1/_R2.fastq.gz
mv ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_S1_L001_R1_001.fastq.gz ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz
mv ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_S1_L001_R2_001.fastq.gz ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz
./scripts/01b-fastp-qc.sh $SAMPLE      # QC + adapter trimming
./scripts/02-alignment.sh $SAMPLE       # FASTQ -> BAM
# ... continue as Path A
```

### Nextflow

A Nextflow DSL2 execution path (v0.5.0) covers post-calling interpretation and clinical analysis — it accepts VCF + BAM from any upstream caller and runs the same pharmacogenomics, annotation, and clinical steps as the bash scripts. Both paths are maintained and produce biologically equivalent results (output file names and report scope may differ).

```bash
# Minimal run — default tools need no external databases
nextflow run main.nf --input samplesheet.csv --reference /path/to/GRCh38.fasta -profile docker

# Enable database-requiring tools (VEP, CPSR, ClinVar, ExpansionHunter)
nextflow run main.nf --input samplesheet.csv --reference /path/to/GRCh38.fasta \
    --tools 'pharmcat,cpic,vcfanno,roh,prs,mito_haplogroup,hla_typing,telomere_hunter,mosdepth,mito_variants,cyrius,html_report,multiqc,vep,slivar,clinical_filter,cpsr,clinvar,expansion_hunter,stranger,pypgx,ancestry' \
    --vep_cache /path/to/vep_cache \
    --pcgr_data /path/to/pcgr_data \
    --vep_cache_cpsr /path/to/vep_cache_113 \
    --clinvar /path/to/clinvar.vcf.gz \
    --clinvar_index /path/to/clinvar.vcf.gz.tbi \
    --expansion_catalog /path/to/variant_catalog.json \
    --hla_dat /path/to/hla.dat \
    --slivar_bin /path/to/slivar \
    --pypgx_bundle /path/to/pypgx-bundle \
    --ancestry_ref /path/to/1kg_common_snps.vcf.gz \
    -profile docker
```

See [docs/nextflow.md](nextflow.md) for samplesheet format, tool selection, sarek integration, and bash vs Nextflow comparison.

## Data from Your Vendor

Different vendors deliver data in different formats. Here's what you need to know:

| Vendor | Format You Get | Pipeline Entry Point | Notes |
|---|---|---|---|
| **Nebula / DNA Complete** | FASTQ + VCF | Path A (FASTQ) or Path C (VCF) | Uses BGI/MGI sequencing |
| **Dante Labs** | FASTQ + BAM + VCF | Any path | Standard Illumina |
| **Sequencing.com** | FASTQ + BAM + VCF | Any path | Standard Illumina |
| **Novogene / BGI** | FASTQ | Path A | BGI read names differ from Illumina but work fine |
| **Illumina DRAGEN (clinical)** | ORA or BAM + VCF | Path D (ORA) or Path B/C | ORA needs decompression first |
| **Oxford Nanopore** | POD5/FAST5 + BAM | [Long-read guide](long-read-guide.md) | minimap2 + Clair3 + Sniffles2 |
| **PacBio HiFi** | HiFi BAM | [Long-read guide](long-read-guide.md) | minimap2 + Clair3 + Sniffles2 |
| **23andMe / Ancestry / MyHeritage** | Genotyping array TSV | Partial (VCF steps only) | Not WGS -- convert to VCF first |

See [docs/vendor-guide.md](vendor-guide.md) for detailed conversion instructions for each vendor.

## Directory Structure

The pipeline expects this layout (created automatically by the scripts):

```
${GENOME_DIR}/
  reference/
    Homo_sapiens_assembly38.fasta      # GRCh38 reference genome
    Homo_sapiens_assembly38.fasta.fai  # FASTA index
  clinvar/
    clinvar.vcf.gz                     # ClinVar database, as downloaded
    clinvar.vcf.gz.tbi                 # ClinVar index
    clinvar_pathogenic_chr.vcf.gz      # chr-prefixed pathogenic subset that step 6 reads
    clinvar_pathogenic_chr.vcf.gz.tbi  # its index
  vep_cache/                           # VEP annotation cache (~30 GB)
  pcgr_data/                           # CPSR/PCGR data bundle (~5 GB)
  ${SAMPLE}/
    fastq/                             # Raw FASTQ files (R1 + R2)
    fastq_trimmed/                     # QC-trimmed FASTQs + fastp reports (step 1b)
    aligned/
      ${SAMPLE}_sorted.bam             # Aligned reads
      ${SAMPLE}_sorted.bam.bai         # BAM index
    vcf/
      ${SAMPLE}.vcf.gz                 # Variant calls
      ${SAMPLE}.vcf.gz.tbi             # VCF index
      ${SAMPLE}.report.html            # PharmCAT report (step 7)
    manta/                             # Structural variants (step 4)
    annotsv/                           # Annotated SVs (step 5)
    clinvar/                           # ClinVar hits (step 6)
    clinical/                          # Clinically filtered variants (step 23)
    vep/                               # Functional annotation (steps 13, 30)
    slivar/                            # Prioritized variants (step 31)
    pypgx/                             # pypgx PGx star alleles (step 32)
    cpsr/                              # Cancer predisposition (step 17)
    mito/                              # Haplogroup + mitochondrial variants (steps 12, 20)
    ...                                # Other analysis directories
```

