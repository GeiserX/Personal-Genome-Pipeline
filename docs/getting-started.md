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
| bash | 4.4+ | Pre-installed on Linux; macOS ships 3.2, install a newer one with `brew install bash` |
| Java *(for a full run)* | 17+ | Any OpenJDK build, e.g. `apt install openjdk-17-jre-headless` or [Adoptium Temurin](https://adoptium.net/) |
| Nextflow *(for a full run)* | 26.04.7 | `curl -s https://get.nextflow.io \| NXF_VER=26.04.7 bash`, see [Full run](#full-run) |
| wget or curl | Any | For downloading references |
| python3 *(optional)* | 3.6+ | Used by long-read alignment (02b) for symlink resolution. Falls back to `readlink -f` on GNU/Linux if absent |

Every analysis tool runs inside Docker: no conda environments, no Python version conflicts, no compilation. Java and Nextflow run the whole pipeline in one command (`run-all.sh`); each step also runs on its own as a script with Docker alone.

### Reference Data (One-Time Downloads)

| Resource | Size | Required For |
|---|---|---|
| GRCh38 reference FASTA + index (no-ALT analysis set) | ~0.8 GB download, ~3 GB unpacked | All steps |
| ClinVar database | ~200 MB | Step 6 (ClinVar screen) |
| VEP cache | ~26 GB | Step 13 (VEP annotation) |
| PCGR/CPSR data bundle + VEP 115 cache | ~31 GB | Step 17 (cancer predisposition) |
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

## Full run

One command runs every step of a default run for one sample. `run-all.sh` checks the setup, writes a one-row samplesheet and starts the [Nextflow pipeline](nextflow.md) (`main.nf`). Nextflow runs the steps in parallel as far as CPUs and memory allow, keeps a log for each task, and on a second run reuses every task whose inputs did not change (`-resume`), so a run that stopped half way continues where it stopped.

It needs Docker, bash 4.4 or later, Java 17 or later and Nextflow 26.04.7, the release CI validates (`NEXTFLOW_VERSION` in `versions.env`). Pin the version when you install, because the plain installer fetches the newest release:

```bash
curl -s https://get.nextflow.io | NXF_VER=26.04.7 bash
sudo mv nextflow /usr/local/bin/
```

Then, once per machine and once per sample:

```bash
export GENOME_DIR=/path/to/your/data
./scripts/setup.sh $GENOME_DIR             # reference, ClinVar, small data files and images (once)
./scripts/validate-setup.sh your_name      # what is installed and what is missing
./scripts/run-all.sh your_name male        # or female: sets the chrX/chrY ploidy and the sex check
```

Without Java or Nextflow, `run-all.sh` stops with exit 2 and prints the install line; the [paths below](#path-a-i-have-fastq-files-raw-reads) run the steps one by one instead. `./scripts/run-all.sh --help` prints the switches.

**Input.** `run-all.sh` uses the first of these that exists under `${GENOME_DIR}/<sample>/`, with its index:

1. `aligned/<sample>_sorted.bam`, with `vcf/<sample>.vcf.gz` when there is one (that VCF is then not called again);
2. `fastq/<sample>_R1.fastq.gz` and `_R2.fastq.gz`;
3. `vcf/<sample>.vcf.gz` alone.

It writes the choice to `<sample>/nextflow/samplesheet.csv` and uses that file again on the next run while every file it names exists and its BAM is the one this call names (`ALIGN_DIR` below), so a rerun finds its tasks in the cache even after the pipeline wrote a BAM and a VCF next to the FASTQ. Delete the samplesheet to choose again.

When the run starts from an existing VCF, the pipeline does not read a gVCF beside it: `vcf/<sample>.g.vcf.gz` from an earlier `03-deepvariant.sh` run stays on disk but unused, so PharmCAT and PRS read only the VCF's variant sites. Run `./scripts/07-pharmacogenomics.sh <sample>` and `./scripts/25-prs.sh <sample>` by hand to get the calls that read the gVCF. When the sample also has its indexed BAM, deleting the VCF and its index is the other way: the BAM is called again, with a gVCF the pipeline reads. With the VCF alone (no BAM, no FASTQ) that would leave no input; such a run lists the steps that read a BAM as `skipped (no BAM)`.

**Switches.** Environment variables, set before the command:

| Switch | What it does | Nextflow parameter |
|---|---|---|
| `SKIP_VALIDATION=true` | does not run `validate-setup.sh` first | |
| `THREADS=N` | at most N CPUs per task (default: the machine's CPU count) | `--max_cpus N` |
| `SKIP_TRIM=true` | aligns the raw reads, without fastp | `--skip_trim true` |
| `INTERVALS="chr20 chr22"` | the regions DeepVariant calls (whole genome by default) | `--intervals` |
| `ALIGN_DIR=dir` | takes the BAM from `<sample>/dir/` instead of `aligned/` | the samplesheet's `bam` |
| `TOOLS=a,b` | runs only these steps, by their `--tools` names; an unknown name stops it with exit 2 | `--tools` |
| `REF_FASTA=path` | another reference inside `GENOME_DIR` (see [Realigning](realignment.md)) | `--reference` |
| `EH_CATALOG=file` | another ExpansionHunter catalog; by default the one inside the ExpansionHunter image, copied once to `reference/expansionhunter_variant_catalog.json` | `--expansion_catalog` |
| `MANTA_CALL_REGIONS=file` | a bgzipped BED of the regions Manta calls | `--manta_call_regions` |
| `ANCESTRY_PANEL=file`, `PGSC_CALC_DIR=dir` | another ancestry panel (step 26 and the percentiles of step 25), another pgsc_calc checkout; both are passed when `setup.sh` installed them, and `ANCESTRY_PANEL=none` scores without the panel | `--ancestry_ref`, `--pgsc_calc` |
| `GRIDSS=true`, `IMPUTATION=true`, `SOMATIC=true` | run steps 4b, 14 and 29 as scripts after the pipeline | |
| `EXTRA_CALLERS=gatk,freebayes,strelka2,octopus` | run the alternative callers 3a to 3d after the pipeline | |
| `BENCHMARK=true` | runs `benchmark-variants.sh` after them; it needs `EXTRA_CALLERS` or a second caller's VCF from an earlier 3a-3d run, and is skipped without one | |

GRIDSS is not part of the SV consensus (step 22): it reports every event as a pair of breakends, which never match the DEL, DUP and INV records of Manta, Delly and CNVpytor.

Without `THREADS` or `--max_cpus`, `run-all.sh` passes the machine's CPU count as `--max_cpus`, and without `--max_memory` it passes the machine's RAM in whole GB (31.GB on a 32 GB Linux machine). Nextflow refuses a task that asks for more than the machine has, and the larger steps ask for 8 CPUs and 32 GB; with the caps they run with less. Options after the sex go to `nextflow run` unchanged, for example `./scripts/run-all.sh your_name male --max_memory 30.GB` to leave 2 GB to the rest of a 32 GB machine, or `--sex_check warn` to go on when the sex inferred from the BAM differs from the one given. [Nextflow Execution](nextflow.md) lists every parameter. `MAX_JOBS` is no longer read: Nextflow starts a task when its CPUs and memory fit, and `--max_cpus` and `--max_memory` cap each task.

**Data.** `run-all.sh` passes each database it finds under `GENOME_DIR`; a step without its data is listed as `skipped (data not installed: ...)`, and the run goes on:

| Step | Data under `GENOME_DIR` | Nextflow parameter |
|---|---|---|
| 6 ClinVar screen | `clinvar/clinvar_pathogenic_chr.vcf.gz` and `.tbi` | `--clinvar`, `--clinvar_index` |
| 8 HLA typing | `hla/IPD-IMGT-HLA_<release>/hla.dat` and `reference/gencode.v50.basic.genes.gtf` (`setup.sh`) | `--hla_dat`, `--hla_genes` |
| 9 ExpansionHunter, 9b Stranger | the catalog above | `--expansion_catalog` |
| 13 VEP, then 23 clinical filter, 30 vcfanno and 31 slivar | `vep_cache/homo_sapiens/116_GRCh38/` (and `annotations/gnomad_v4.1_constraint.tsv` for slivar) | `--vep_cache` (`--gnomad_constraint`) |
| 30 vcfanno | the score files of [step 30](30-vcfanno.md) under `annotations/`, with their `.tbi` | `--cadd_snv`, `--spliceai_snv`, `--revel`, `--alphamissense`, ... |
| 17 CPSR | `vep_cache/homo_sapiens/115_GRCh38/` and `pcgr_data/<bundle>/data/` | `--vep_cache_cpsr`, `--pcgr_data` |
| 18 CNVpytor | `reference/cnvpytor/gc_hg38.pytor` | `--cnvpytor_resources` |
| 5 AnnotSV | `annotsv_annotations/Annotations_Human/` | `--annotsv_annotations` |
| 32 pypgx | `reference/pypgx-bundle/` | `--pypgx_bundle` |
| 25 PRS | `prs_scores/*.txt.gz` (run `./scripts/25-prs.sh <sample>` once to download the scoring files) | `--pgs_scoring` |
| 10 TelomereHunter, 19 Delly | `reference/cytoBand.hg38.txt`, `reference/delly_human.hg38.excl.tsv` (optional) | `--cytoband`, `--delly_exclude` |

**Where things land.** Results go to `${GENOME_DIR}/<sample>/`, in the folders of the pipeline's [output structure](nextflow.md#output-structure). Most match the single scripts' folders; the pipeline writes PharmCAT to `pharmcat/`, ROH to `roh/`, depth to `coverage/` and HLA types to `hla/`, where the scripts use `vcf/`, `vcf/`, `mosdepth/` and `hla_t1k/`. The reports read both. At the end `run-all.sh` renders the full HTML report (`<sample>_report.html`, step 24), the text report (`<sample>_report.txt`) and the `summary.json` both are made from.

| Log | Where |
|---|---|
| Nextflow's own log | `<sample>/nextflow/.nextflow.log` |
| Each task's command, output and exit code | `<sample>/nextflow/work/<hash>/` (`.command.sh`, `.command.log`, `.exitcode`); the trace in `${GENOME_DIR}/pipeline_info/` gives each task's hash, the start of its folder's name |
| The script steps and the reports | `<sample>/logs/<script>.log` |
| How each step ended in this run | `<sample>/logs/run_status.tsv` (the reports mark an older result of a skipped step as stale) |

The `work/` folder holds a copy of every intermediate file. Once the run succeeded and the results look right, delete `<sample>/nextflow/work` to free the space; the next run then starts from scratch.

**After a failure** run the same command again: `-resume` reruns only the tasks that did not finish. **To rerun one step by hand**, run its script, for example `./scripts/06-clinvar-screen.sh your_name` after a ClinVar update, then `./scripts/24-html-report.sh your_name` and `./scripts/generate-report.sh your_name` to refresh the reports. A script reads the folders the scripts write; where the pipeline's folder differs (above), a script that reads another step's output may need that step run by hand first.

### Path A: I Have FASTQ Files (Raw Reads)

The paths below run the steps one by one, with Docker alone. Most common if you downloaded data from Nebula, Dante Labs, Novogene, BGI, or any sequencing provider.

```bash
# 1. Set your data directory (where your FASTQ files are)
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name

# 2. Download the GRCh38 no-ALT reference genome (~0.8 GB, ~3 GB unpacked),
#    ClinVar and the images; see docs/00-reference-setup.md for details
./scripts/setup.sh $GENOME_DIR

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
# It must be aligned to this pipeline's reference; validate-setup.sh checks
# its header against it. A BAM aligned to another reference (most vendor BAMs)
# is realigned first: docs/realignment.md.
./scripts/validate-setup.sh $SAMPLE
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

### Nextflow directly

`run-all.sh` writes the samplesheet and the parameters for one sample in `GENOME_DIR`. To run several samples at once, or with your own folders, call the pipeline yourself. A samplesheet row starts from FASTQ, a BAM, or a VCF with an optional BAM.

```bash
# Default tools only (prs and vcfanno are skipped, with a warning, until their files are set)
nextflow run main.nf --input samplesheet.csv --reference /path/to/GRCh38_no_alt_analysis_set.fasta \
    --outdir ./results -profile docker

# With the tools that need databases
nextflow run main.nf --input samplesheet.csv --reference /path/to/GRCh38_no_alt_analysis_set.fasta \
    --tools 'pharmcat,cpic,vcfanno,roh,prs,mito_haplogroup,hla_typing,telomere_hunter,mosdepth,mito_variants,html_report,multiqc,vep,slivar,clinical_filter,cpsr,clinvar,expansion_hunter,stranger,pypgx,ancestry' \
    --vep_cache /path/to/vep_cache \
    --pcgr_data /path/to/pcgr_data/20260620 \
    --vep_cache_cpsr /path/to/vep_cache \
    --clinvar /path/to/clinvar/clinvar_pathogenic_chr.vcf.gz \
    --clinvar_index /path/to/clinvar/clinvar_pathogenic_chr.vcf.gz.tbi \
    --expansion_catalog /path/to/variant_catalog.json \
    --hla_dat /path/to/hla.dat \
    --hla_genes /path/to/gencode.v50.basic.genes.gtf \
    --pypgx_bundle /path/to/pypgx-bundle \
    --ancestry_ref /path/to/reference/pgsc_calc/pgsc_1000G_v1.tar.zst \
    --pgs_scoring /path/to/prs_scores \
    --pgsc_calc /path/to/tools/pgsc_calc-<release> \
    --outdir ./results -profile docker
```

See [Nextflow Execution](nextflow.md) for the samplesheet format, tool selection, sarek output, and where the pipeline and the scripts differ.

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
    GRCh38_no_alt_analysis_set.fasta      # GRCh38 reference genome (no ALT contigs)
    GRCh38_no_alt_analysis_set.fasta.fai  # FASTA index
  clinvar/
    clinvar.vcf.gz                     # ClinVar database, as downloaded
    clinvar.vcf.gz.tbi                 # ClinVar index
    clinvar_pathogenic_chr.vcf.gz      # chr-prefixed pathogenic subset that step 6 reads
    clinvar_pathogenic_chr.vcf.gz.tbi  # its index
  vep_cache/                           # VEP annotation cache (~30 GB)
  pcgr_data/                           # CPSR/PCGR data bundle (~7 GB)
  ${SAMPLE}/
    fastq/                             # Raw FASTQ files (R1 + R2)
    nextflow/                          # run-all.sh: samplesheet.csv, .nextflow.log, work/
    logs/                              # run-all.sh: run_status.tsv, script and report logs
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

