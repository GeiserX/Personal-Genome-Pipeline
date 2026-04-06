# Personal Genome Pipeline Documentation

Analyze your own whole genome sequencing (WGS) data on consumer hardware. 31 analysis steps, all running locally in Docker.

[Back to GitHub repository](https://github.com/GeiserX/Personal-Genome-Pipeline)

---

## Getting Started

- [Reference Setup](00-reference-setup.md) -- download reference genome and databases
- [Hardware Requirements](hardware-requirements.md) -- what you need
- [Vendor Guide](vendor-guide.md) -- how to get your data from each provider
- [Quick Test](quick-test.md) -- verify your setup with public data
- [Validate Setup](../scripts/validate-setup.sh) -- pre-flight check script

## Pipeline Steps

### Core Analysis (Steps 1-20)

| # | Step | Doc |
|---|---|---|
| 1 | ORA to FASTQ (Illumina decompression) | [01-ora-to-fastq.md](01-ora-to-fastq.md) |
| 1b | fastp QC + Trimming | [01b-fastp-qc.md](01b-fastp-qc.md) |
| 2 | Alignment (minimap2 + samtools) | [02-alignment.md](02-alignment.md) |
| 3 | Variant Calling (DeepVariant) | [03-variant-calling.md](03-variant-calling.md) |
| 4 | Structural Variants (Manta) | [04-structural-variants.md](04-structural-variants.md) |
| 5 | SV Annotation (AnnotSV) | [05-annotsv.md](05-annotsv.md) |
| 6 | ClinVar Screen | [06-clinvar-screen.md](06-clinvar-screen.md) |
| 7 | Pharmacogenomics (PharmCAT) | [07-pharmacogenomics.md](07-pharmacogenomics.md) |
| 8 | HLA Typing (T1K) | [08-hla-typing.md](08-hla-typing.md) |
| 9 | STR Expansions (ExpansionHunter) | [09-str-expansions.md](09-str-expansions.md) |
| 10 | Telomere Length (TelomereHunter) | [10-telomere-analysis.md](10-telomere-analysis.md) |
| 11 | ROH Analysis (bcftools roh) | [11-roh-analysis.md](11-roh-analysis.md) |
| 12 | Mito Haplogroup (haplogrep3) | [12-mito-haplogroup.md](12-mito-haplogroup.md) |
| 13 | VEP Annotation | [13-vep-annotation.md](13-vep-annotation.md) |
| 14 | Imputation Prep | [14-imputation-prep.md](14-imputation-prep.md) |
| 15 | SV Quality (duphold) | [15-duphold.md](15-duphold.md) |
| 16 | Coverage QC (indexcov) | [16-indexcov.md](16-indexcov.md) |
| 16b | Coverage Statistics (mosdepth) | [16b-mosdepth.md](16b-mosdepth.md) |
| 17 | Cancer Predisposition (CPSR) | [17-cpsr.md](17-cpsr.md) |
| 18 | CNV Calling (CNVnator) | [18-cnvnator.md](18-cnvnator.md) |
| 19 | SV Calling (Delly) | [19-delly.md](19-delly.md) |
| 20 | Mitochondrial (GATK Mutect2) | [20-mtoolbox.md](20-mtoolbox.md) |

### Post-Processing (Steps 21-29)

| # | Step | Doc |
|---|---|---|
| 21 | CYP2D6 Star Alleles (Cyrius) | [21-cyrius.md](21-cyrius.md) |
| 22 | SV Consensus Merge | [22-survivor-merge.md](22-survivor-merge.md) |
| 23 | Clinical Variant Filter | [23-clinical-filter.md](23-clinical-filter.md) |
| 24 | HTML Summary Report | [24-html-report.md](24-html-report.md) |
| 25 | Polygenic Risk Scores | [25-prs.md](25-prs.md) |
| 26 | Ancestry PCA | [26-ancestry.md](26-ancestry.md) |
| 27 | CPIC Drug Recommendations | [27-cpic-lookup.md](27-cpic-lookup.md) |
| 28 | MultiQC Aggregation | [28-multiqc.md](28-multiqc.md) |
| 29 | Somatic Variants (Mutect2) | [29-mutect2-somatic.md](29-mutect2-somatic.md) |

## Guides

- [Interpreting Results](interpreting-results.md) -- what your results mean
- [Multi-Sample Analysis](multi-sample.md) -- comparing two or more genomes
- [Troubleshooting](troubleshooting.md) -- common issues and fixes
- [Lessons Learned](lessons-learned.md) -- failures encountered during development
- [Glossary](glossary.md) -- genomics terminology
- [Resources](resources.md) -- free courses, databases, and tools
- [Long-Read Guide](long-read-guide.md) -- ONT and PacBio support
- [WES Guide](wes-guide.md) -- whole exome sequencing entry path
