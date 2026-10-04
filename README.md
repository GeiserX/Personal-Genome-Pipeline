<p align="center">
  <img src="docs/images/banner.svg" alt="Personal Genome Pipeline banner" width="900"/>
</p>

<h1 align="center">Personal Genome Pipeline</h1>

<p align="center">
  <strong>Analyze your own whole genome sequencing (WGS) data on consumer hardware.</strong><br>
  No cloud accounts, no subscriptions, no bioinformatics degree required.
</p>

<p align="center">
  <a href="https://github.com/GeiserX/Personal-Genome-Pipeline/releases"><img src="https://img.shields.io/github/v/release/GeiserX/Personal-Genome-Pipeline?style=flat-square" alt="Release"></a>
  <a href="https://github.com/GeiserX/Personal-Genome-Pipeline/actions/workflows/lint.yml"><img src="https://img.shields.io/github/actions/workflow/status/GeiserX/Personal-Genome-Pipeline/lint.yml?style=flat-square&label=CI" alt="CI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/github/license/GeiserX/Personal-Genome-Pipeline?style=flat-square" alt="License"></a>
  <a href="https://github.com/GeiserX/Personal-Genome-Pipeline/stargazers"><img src="https://img.shields.io/github/stars/GeiserX/Personal-Genome-Pipeline?style=flat-square&logo=github" alt="GitHub Stars"></a>
</p>

This pipeline takes raw sequencing data (FASTQ/BAM/VCF) from any vendor and runs the [analysis steps of a default run](https://geiserx.github.io/Personal-Genome-Pipeline/pipeline-overview/#what-a-default-run-covers) to produce a full genomic profile: variant calling, pharmacogenomics, structural variants, cancer predisposition screening, polygenic risk scores, ancestry estimation, telomere length, mitochondrial analysis, and more. Everything runs locally in Docker containers with resource limits so it won't crash your machine.

**Time:** 6-12 hours per sample on a 16-core desktop | **Disk:** 500 GB minimum per sample | **Cost:** Free (you just need your data)

## Features

- Takes FASTQ, BAM, VCF or Illumina ORA from any vendor (Nebula, Dante Labs, Sequencing.com, Novogene, DRAGEN), plus long reads from Nanopore and PacBio HiFi.
- Calls SNPs and indels with DeepVariant, and structural and copy number variants with Manta, Delly and CNVpytor merged into a consensus.
- Screens ClinVar and runs CPSR cancer predisposition panels, VEP annotation with CADD, SpliceAI, REVEL and AlphaMissense, and slivar prioritization.
- Pharmacogenomics with PharmCAT, pypgx (23 genes, CYP2D6 SVs), Cyrius and CPIC drug recommendations.
- Repeat expansions, HLA typing, telomere length, mitochondrial haplogroup and heteroplasmy, ROH, ancestry and polygenic risk scores.
- Every tool runs in a Docker container with CPU and memory limits, pinned by tag or digest in `versions.env`. One exception: Cyrius is installed from PyPI at run time (version and dependencies pinned). No script uploads your data; [a few steps download public files](https://geiserx.github.io/Personal-Genome-Pipeline/why-local/#network-calls-during-a-run) during a run.
- One Nextflow DSL2 pipeline from FASTQ, BAM or VCF to the report, and every step also as a bash script you can run on its own.
- Ends in an HTML report and a MultiQC summary. Alternative callers (GATK, FreeBayes, Strelka2, Octopus, BWA-MEM2, TIDDIT, GRIDSS) are there for benchmarking.

## Quick start

Start from FASTQ (other entry points in [Getting started](https://geiserx.github.io/Personal-Genome-Pipeline/getting-started/)):

```bash
export GENOME_DIR=/path/to/your/data SAMPLE=your_name
./scripts/validate-setup.sh $SAMPLE
./scripts/02-alignment.sh $SAMPLE && ./scripts/03-deepvariant.sh $SAMPLE && ./scripts/06-clinvar-screen.sh $SAMPLE && ./scripts/07-pharmacogenomics.sh $SAMPLE
```

Download the GRCh38 reference and databases first with [reference setup](https://geiserx.github.io/Personal-Genome-Pipeline/00-reference-setup/). Try a small public dataset with the [quick test](https://geiserx.github.io/Personal-Genome-Pipeline/quick-test/).

## Documentation

The full documentation is at **https://geiserx.github.io/Personal-Genome-Pipeline/**, one page per pipeline step included.

- [Getting started](https://geiserx.github.io/Personal-Genome-Pipeline/getting-started/): prerequisites, platform notes (macOS, WSL2, Unraid), the FASTQ, BAM, VCF and ORA entry paths, directory layout
- [Quick test](https://geiserx.github.io/Personal-Genome-Pipeline/quick-test/): verify the setup on public data before your own
- [Hardware and storage requirements](https://geiserx.github.io/Personal-Genome-Pipeline/hardware-requirements/): download sizes, per-step runtime, memory and disk figures
- [Reference data setup](https://geiserx.github.io/Personal-Genome-Pipeline/00-reference-setup/): the GRCh38 reference (the no-ALT analysis set) and every database
- [Realigning after a reference change](https://geiserx.github.io/Personal-Genome-Pipeline/realignment/): moving a BAM aligned to another GRCh38 file onto the pipeline's reference
- [Vendor compatibility guide](https://geiserx.github.io/Personal-Genome-Pipeline/vendor-guide/): what each provider delivers and how to get it
- [Pipeline overview](https://geiserx.github.io/Personal-Genome-Pipeline/pipeline-overview/): every step with its tool and image, which ones a default run includes, and the page for each step
- [Nextflow](https://geiserx.github.io/Personal-Genome-Pipeline/nextflow/): the workflow runner, parallel steps and resume
- [Interpreting your results](https://geiserx.github.io/Personal-Genome-Pipeline/interpreting-results/): what each report means and what to do with it
- [Multi-sample comparison](https://geiserx.github.io/Personal-Genome-Pipeline/multi-sample/): partners, siblings, parents
- [Long-read guide](https://geiserx.github.io/Personal-Genome-Pipeline/long-read-guide/): Nanopore and PacBio HiFi
- [WES guide](https://geiserx.github.io/Personal-Genome-Pipeline/wes-guide/): whole exome input
- [Chip data guide](https://geiserx.github.io/Personal-Genome-Pipeline/chip-data-guide/): 23andMe, MyHeritage and AncestryDNA files
- [Variant caller benchmarking](https://geiserx.github.io/Personal-Genome-Pipeline/benchmarking/): the alternative callers compared
- [Common issues and FAQ](https://geiserx.github.io/Personal-Genome-Pipeline/faq/) and [Troubleshooting](https://geiserx.github.io/Personal-Genome-Pipeline/troubleshooting/)
- [Lessons learned](https://geiserx.github.io/Personal-Genome-Pipeline/lessons-learned/), [Tool selection rationale](https://geiserx.github.io/Personal-Genome-Pipeline/tool-rationale/), [Glossary](https://geiserx.github.io/Personal-Genome-Pipeline/glossary/), [Resources](https://geiserx.github.io/Personal-Genome-Pipeline/resources/)
- [Why run locally?](https://geiserx.github.io/Personal-Genome-Pipeline/why-local/): cost comparison, privacy and security
- [Contributing](CONTRIBUTING.md): each step is a standalone script with its own page

## Disclaimer

This pipeline is for **educational and research purposes only**. It is not a medical device and has not been clinically validated. Genomic findings should always be discussed with a qualified healthcare professional before making any medical decisions. The authors are not responsible for any actions taken based on pipeline output.

Your genome data is sensitive personal information. This pipeline runs locally and no script uploads it. The only exception is your choice: step 14 prepares files for an imputation server, and sending them there is up to you. Keep your data secure.

## License

[GPL-3.0-or-later](LICENSE)
