# Genomics Pipeline

End-to-end pipeline for processing whole genome sequencing (WGS) data into clinically actionable reports. Designed to run on consumer hardware using Docker containers with resource limits.

## Overview

Takes raw sequencing data (ORA/FASTQ) through alignment, variant calling, and 12+ analysis panels to produce a comprehensive genomic profile. Every step runs in Docker with `--cpus` and `--memory` limits to avoid crashing the host.

## Pipeline Steps

| # | Step | Tool | Docker Image | Input | Output |
|---|---|---|---|---|---|
| 1 | [ORA to FASTQ](docs/01-ora-to-fastq.md) | orad (Illumina) | `orad` binary | `.ora` files | `.fastq.gz` |
| 2 | [Alignment](docs/02-alignment.md) | minimap2 + samtools | `staphb/samtools:1.20` | `.fastq.gz` + reference | `sorted.bam` + `.bai` |
| 3 | [Variant Calling](docs/03-variant-calling.md) | DeepVariant | `google/deepvariant:1.6.0` | `sorted.bam` | `.vcf.gz` (SNPs/indels) |
| 4 | [Structural Variants](docs/04-structural-variants.md) | Manta | `quay.io/biocontainers/manta` | `sorted.bam` | `diploidSV.vcf.gz` |
| 5 | [SV Annotation](docs/05-annotsv.md) | AnnotSV | `getwilds/annotsv:latest` | `diploidSV.vcf.gz` | `_sv_annotated.tsv` |
| 6 | [ClinVar Screen](docs/06-clinvar-screen.md) | bcftools isec | `staphb/bcftools:1.21` | `.vcf.gz` + ClinVar DB | pathogenic hits |
| 7 | [Pharmacogenomics](docs/07-pharmacogenomics.md) | PharmCAT | `pgkb/pharmcat:2.15.5` | `.vcf.gz` | PGx report HTML |
| 8 | [HLA Typing](docs/08-hla-typing.md) | T1K | `quay.io/biocontainers/t1k:1.0.9` | `sorted.bam` | HLA alleles |
| 9 | [STR Expansions](docs/09-str-expansions.md) | ExpansionHunter | `weisburd/expansionhunter` | `sorted.bam` + reference | repeat genotypes |
| 10 | [Telomere Length](docs/10-telomere-analysis.md) | TelomereHunter | `lgalarno/telomerehunter:latest` | `sorted.bam` | telomere content |
| 11 | [ROH Analysis](docs/11-roh-analysis.md) | bcftools roh | `staphb/bcftools:1.21` | `.vcf.gz` | ROH segments |
| 12 | [Mitochondrial Haplogroup](docs/12-mito-haplogroup.md) | haplogrep3 | `genepi/haplogrep3` | `.vcf.gz` (chrM) | haplogroup |
| 13 | [Functional Annotation](docs/13-vep-annotation.md) | VEP | `ensemblorg/ensembl-vep:release_112.0` | `.vcf.gz` | annotated VCF |
| 14 | [Imputation Prep](docs/14-imputation-prep.md) | bcftools | `staphb/bcftools:1.21` | `.vcf.gz` | chr-split VCFs |
| 15 | [SV Quality](docs/15-duphold.md) | duphold | `brentp/duphold:latest` | `diploidSV.vcf.gz` + BAM | quality-annotated SV VCF |
| 16 | [Coverage QC](docs/16-indexcov.md) | indexcov | `quay.io/biocontainers/goleft` | `sorted.bam` | coverage + sex check |

## Prerequisites

- Docker (amd64 host recommended; works on any Linux server with Docker)
- GRCh38 reference genome (`Homo_sapiens_assembly38.fasta` + `.fai`)
- ClinVar database (`clinvar.vcf.gz` + `.tbi`)
- 20+ cores and 32+ GB RAM recommended (each step is resource-limited)
- ~200GB disk per sample (BAM ~30-40GB, VCF ~100MB, intermediates ~50GB)

## Quick Start

```bash
# Set your data directory
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_sample
export SEX=male  # or female

# Run each step (see docs/ for full details)
./scripts/03-deepvariant.sh $SAMPLE
./scripts/04-manta.sh $SAMPLE
./scripts/05-annotsv.sh $SAMPLE
# ... etc
```

## Lessons Learned

See [docs/lessons-learned.md](docs/lessons-learned.md) for every failure encountered and how it was resolved.

## License

MIT
