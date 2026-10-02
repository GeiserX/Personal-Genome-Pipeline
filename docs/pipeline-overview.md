# Pipeline overview

## What You Get

| Category | What It Finds | Steps |
|---|---|---|
| **Variant Calling** | SNPs, indels, structural variants, copy number variants | 3, 4, 4b, 18, 19 |
| **Clinical Screening** | Pathogenic variants, carrier status, cancer predisposition (CPSR panels) | 6, 17 |
| **Pharmacogenomics** | Drug-gene interactions (23+ genes, CYP2C19, CYP2D6 SV, DPYD, etc.) | 7, 21, 27, 32 |
| **Structural Variants** | Deletions, duplications, inversions, translocations (4 callers + consensus) | 4, 4b, 5, 15, 18, 19, 22 |
| **Functional Annotation** | Impact prediction for every variant (VEP + CADD, SpliceAI, REVEL, AlphaMissense) | 13, 30 |
| **Variant Prioritization** | Rare deleterious variants, compound hets, gene constraint filtering | 31 |
| **Repeat Expansions** | Huntington's, Fragile X, ALS, and 50+ other repeat expansion disorders | 9 |
| **Ancestry & Haplogroups** | Mitochondrial haplogroup, consanguinity check, ancestry SNP intersection | 11, 12, 26 |
| **Telomere Length** | Relative telomere content estimation from WGS reads | 10 |
| **Mitochondrial** | Heteroplasmy detection, mitochondrial disease variants | 12, 20 |
| **Polygenic Risk** | Risk scores for 10 common conditions (CAD, T2D, cancers, etc.) | 25 |
| **Quality Control** | Adapter trimming, coverage statistics, aggregated QC report, sex check, SV filtering | 1b, 15, 16, 16b, 28 |

## Pipeline Overview

```mermaid
graph LR
    FASTQ["FASTQ"] --> fastp["fastp<br/><small>QC + trim</small>"]
    fastp --> align["minimap2<br/><small>Alignment</small>"]
    align --> BAM["Sorted BAM"]
    BAM --> DV["DeepVariant<br/><small>SNPs + indels</small>"]
    DV --> VCF["VCF"]

    %% VCF-based analyses
    VCF --> clinvar["ClinVar Screen"]
    VCF --> pharmcat["PharmCAT<br/><small>PGx</small>"]
    pharmcat --> cpic["CPIC<br/><small>Drug recs</small>"]
    VCF --> vep["VEP<br/><small>Annotation</small>"]
    vep --> vcfanno["vcfanno<br/><small>CADD / SpliceAI<br/>REVEL / AlphaMissense</small>"]
    vcfanno --> slivar["slivar<br/><small>Prioritization</small>"]
    vcfanno --> clinical["Clinical Filter"]
    VCF --> cpsr["CPSR<br/><small>Cancer predisposition</small>"]
    VCF --> roh["ROH Analysis"]
    VCF --> prs["PRS<br/><small>Polygenic risk</small>"]
    VCF --> ancestry["Ancestry SNPs"]

    %% BAM-based analyses
    BAM --> manta["Manta<br/><small>SVs</small>"]
    BAM --> delly["Delly<br/><small>SVs</small>"]
    BAM --> cnvpytor["CNVpytor<br/><small>CNVs</small>"]
    manta --> duphold["duphold"]
    duphold --> annotsv["AnnotSV"]
    manta --> consensus["SV Consensus"]
    delly --> consensus
    cnvpytor --> consensus

    BAM --> eh["ExpansionHunter<br/><small>STRs</small>"]
    BAM --> pypgx["pypgx<br/><small>23-gene PGx<br/>+ CYP2D6 SV</small>"]
    BAM --> cyrius["Cyrius<br/><small>CYP2D6</small>"]
    BAM --> telomere["TelomereHunter"]
    BAM --> coverage["mosdepth<br/>+ indexcov"]
    BAM --> mito["Mutect2<br/><small>Mitochondrial</small>"]
    BAM --> haplo["Haplogrep3<br/><small>mtDNA haplogroup</small>"]

    %% Reporting
    clinical --> report["HTML Report<br/>+ MultiQC"]
    slivar --> report
    clinvar --> report
    pharmcat --> report
    cpsr --> report

    %% Styling
    classDef input fill:#0ea5e9,stroke:#0284c7,color:#fff
    classDef core fill:#8b5cf6,stroke:#7c3aed,color:#fff
    classDef analysis fill:#10b981,stroke:#059669,color:#fff
    classDef sv fill:#f59e0b,stroke:#d97706,color:#fff
    classDef annotation fill:#ec4899,stroke:#db2777,color:#fff
    classDef report fill:#ef4444,stroke:#dc2626,color:#fff

    class FASTQ,BAM,VCF input
    class fastp,align,DV core
    class clinvar,pharmcat,cpic,cpsr,eh,roh,prs,ancestry,pypgx,cyrius,telomere,coverage,mito,haplo analysis
    class manta,delly,cnvpytor,consensus,duphold,annotsv sv
    class vep,vcfanno,slivar,clinical annotation
    class report report
```

### All Steps

| # | Step | Tool | Docker Image | Runtime | Required? |
|---|---|---|---|---|---|
| 1 | [ORA to FASTQ](01-ora-to-fastq.md) | orad | `orad` binary | ~30 min | Only for Illumina ORA files |
| 1b | [QC & Trimming](01b-fastp-qc.md) | fastp | `FASTP_IMAGE` | ~15-30 min | Recommended |
| 2 | [Alignment](02-alignment.md) | minimap2 + samtools | `MINIMAP2_IMAGE` + `SAMTOOLS_IMAGE` | ~1-2 hr | Yes (if starting from FASTQ) |
| 3 | [Variant Calling](03-variant-calling.md) | DeepVariant | `DEEPVARIANT_IMAGE` | ~2-4 hr | Yes |
| 4 | [Structural Variants](04-structural-variants.md) | Manta | `MANTA_IMAGE` | ~20 min | Recommended |
| 5 | [SV Annotation](05-annotsv.md) | AnnotSV | `ANNOTSV_IMAGE` | ~10 min | If step 4 run |
| 6 | [ClinVar Screen](06-clinvar-screen.md) | bcftools isec | `BCFTOOLS_IMAGE` | ~5 min | Yes |
| 7 | [Pharmacogenomics](07-pharmacogenomics.md) | PharmCAT | `PHARMCAT_IMAGE` | ~10 min | Yes |
| 8 | [HLA Typing](08-hla-typing.md) | T1K | `T1K_IMAGE` | ~30 min | Optional |
| 9 | [STR Expansions](09-str-expansions.md) | ExpansionHunter | `EXPANSIONHUNTER_IMAGE` | ~15 min | Recommended |
| 9b | [STR Annotation](09b-stranger.md) | Stranger | `STRANGER_IMAGE` | ~1 min | If step 9 run |
| 10 | [Telomere Length](10-telomere-analysis.md) | TelomereHunter | `lgalarno/telomerehunter` (digest-pinned) | ~1 hr | Optional |
| 11 | [ROH Analysis](11-roh-analysis.md) | bcftools roh | `BCFTOOLS_IMAGE` | ~5 min | Recommended |
| 12 | [Mito Haplogroup](12-mito-haplogroup.md) | haplogrep3 | `jtb114/haplogrep3` (digest-pinned) | ~1 min | Optional |
| 13 | [VEP Annotation](13-vep-annotation.md) | VEP | `VEP_IMAGE` | ~2-4 hr | Recommended |
| 14 | [Imputation Prep](14-imputation-prep.md) | bcftools | `BCFTOOLS_IMAGE` | ~10 min | Optional |
| 15 | [SV Quality](15-duphold.md) | duphold | `DUPHOLD_IMAGE` | ~20 min | If step 4 run |
| 16 | [Coverage QC](16-indexcov.md) | indexcov | `GOLEFT_IMAGE` | ~5 sec | Recommended |
| 16b | [Coverage Stats](16b-mosdepth.md) | mosdepth | `MOSDEPTH_IMAGE` | ~10 min | Recommended |
| 17 | [Cancer Predisposition](17-cpsr.md) | CPSR | `PCGR_IMAGE` | ~30-60 min | Recommended |
| 18 | [CNV Calling](18-cnvpytor.md) | CNVpytor | `CNVPYTOR_IMAGE` | ~1-3 hr | Optional |
| 19 | [SV Calling (Delly)](19-delly.md) | Delly | `DELLY_IMAGE` | ~2-4 hr | Optional |
| 20 | [Mitochondrial](20-mtoolbox.md) | GATK Mutect2 | `GATK_IMAGE` | ~15-30 min | Optional |

#### Post-Processing Steps

These run after the core pipeline completes and combine outputs from earlier steps.

| # | Step | Tool | Docker Image | Runtime | Required? |
|---|---|---|---|---|---|
| 21 | [CYP2D6 Star Alleles](21-cyrius.md) | Cyrius | `PYTHON_IMAGE` | ~10 min | Experimental |
| 22 | [SV Consensus Merge](22-survivor-merge.md) | bcftools | `BCFTOOLS_IMAGE` | ~5 min | Experimental |
| 23 | [Clinical Filter](23-clinical-filter.md) | bcftools +split-vep | `BCFTOOLS_IMAGE` | ~5-10 min | If step 13 run |
| 24 | [HTML Report](24-html-report.md) | bash + bcftools | `BCFTOOLS_IMAGE` | ~1-3 min | Recommended |
| 25 | [Polygenic Risk Scores](25-prs.md) | plink2 | `PLINK2_IMAGE` | ~30 min | Exploratory |
| 26 | [Ancestry SNPs](26-ancestry.md) | plink2 | `PLINK2_IMAGE` | ~30-60 min | Experimental |
| 27 | [CPIC Recommendations](27-cpic-lookup.md) | Python + CPIC | `PYTHON_IMAGE` | ~5 min | If step 7 run |
| 28 | [MultiQC Report](28-multiqc.md) | MultiQC | `MULTIQC_IMAGE` | ~1 min | Recommended |
| 29 | [Somatic Variants](29-mutect2-somatic.md) | GATK Mutect2 | `GATK_IMAGE` | ~2-6 hr | Experimental |
| 30 | [Annotation Enrichment](30-vcfanno.md) | vcfanno | `VCFANNO_IMAGE` | ~5-15 min | If step 13 run |
| 31 | [Variant Prioritization](31-slivar.md) | slivar | `SLIVAR_IMAGE` | ~5-10 min | If step 13 run |
| 32 | [pypgx Pharmacogenomics](32-pypgx.md) | pypgx | `PYPGX_IMAGE` | ~20-40 min | Recommended |

**Minimum useful run:** Steps 2, 3, 6, 7 (alignment + variant calling + ClinVar + PharmCAT) = ~4-6 hours.
**Full analysis:** All 34 default steps = ~12-20 hours (step 29 somatic calling is opt-in via `SOMATIC=true`). Steps 4/18/19 and 10/12/20 can run in parallel.

#### Alternative Tools (Benchmarking)

No single variant caller is universally best. The pipeline includes alternative tools that output to separate directories so you can compare results without overwriting the defaults.

| Script | Tool | Alternative To | Output Directory |
|---|---|---|---|
| [02a](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/02a-alignment-bwamem2.sh) | BWA-MEM2 | minimap2 (step 2) | `aligned_bwamem2/` |
| [03a](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/03a-gatk-haplotypecaller.sh) | GATK HaplotypeCaller | DeepVariant (step 3) | `vcf_gatk/` |
| [03b](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/03b-freebayes.sh) | FreeBayes | DeepVariant (step 3) | `vcf_freebayes/` |
| [04a](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/04a-tiddit.sh) | TIDDIT | Manta (step 4) | `sv_tiddit/` |
| [03c](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/03c-strelka2-germline.sh) | Strelka2 | DeepVariant (step 3) | `vcf_strelka2/` |
| [03d](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/03d-octopus.sh) | Octopus | DeepVariant (step 3) | `vcf_octopus/` |
| [04b](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/04b-gridss.sh) | GRIDSS | Manta (step 4) | `sv_gridss/` |
| [benchmark](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/benchmark-variants.sh) | bcftools isec / hap.py | — | `benchmark/` |

See [docs/benchmarking.md](benchmarking.md) for how to run and interpret results, and [docs/tool-rationale.md](tool-rationale.md) for why each default was chosen.

