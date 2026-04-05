# Awesome-Bioinformatics Triage

Tools from [awesome-bioinformatics](https://github.com/danielecook/Awesome-Bioinformatics) evaluated against this pipeline's scope (germline WGS, Docker-based, single-sample default with multi-sample planned).

Last reviewed: 2026-04-03

---

## Strong Fits

| Tool | Category | Why It Fits | Integration Point | Blockers | Roadmap |
|---|---|---|---|---|---|
| **Somalier** | QC / Sample Identity | Ultra-fast relatedness and sample-swap detection from BAM/VCF. Replaces ad-hoc sex-check (indexcov) with proper identity QC. Multi-sample aware. | After alignment (step 2) | None — single binary, small Docker image | v0.3.0 |
| **GLNexus** | Joint Genotyping | Merges per-sample gVCFs into a joint-called cohort VCF. Required for proper multi-sample analysis (trios, families). DeepVariant already supports gVCF output. | New step after per-sample variant calling | Requires switching DV to `--output_gvcf`; significant architectural change | v0.3.0 |
| **vcfanno** | Annotation | Adds arbitrary annotation tracks (gnomAD, CADD, SpliceAI, custom BEDs) to VCFs via TOML config. Much faster than VEP for bulk annotation; complements VEP rather than replacing it. | Alongside or after VEP (step 13) | Need to curate annotation sources; large reference downloads | v0.4.0 |
| **mosdepth** | QC / Coverage | Fast per-base and per-region depth statistics from BAM. More detailed than indexcov (which only reads the index). Produces coverage distributions, thresholds, and BED outputs. | After alignment (step 2), alongside indexcov | None — drop-in Docker image | v0.6.0 |
| **FastQC + MultiQC** | QC / Reporting | FastQC: per-FASTQ quality metrics (adapter content, GC bias, duplication). MultiQC: aggregates all QC outputs (FastQC, samtools flagstat, mosdepth, etc.) into a single HTML report. | Before alignment (step 1.5) + after all steps | None — standard Docker images | v0.6.0 |
| **Octopus** | Variant Calling | Haplotype-aware Bayesian caller with built-in somatic mode. Could serve as a 5th caller for benchmarking, or replace FreeBayes (which is slow and memory-heavy). | Alternative to step 3 (like 03a/03b/03c) | Large Docker image; evaluate against DV/GATK first | v0.6.0 |
| **GRIDSS** | Structural Variants | Assembly-based SV caller. Better sensitivity for complex rearrangements than Manta/Delly. Would strengthen SV consensus (step 22). | Alternative to step 4, feeds into SV consensus | Heavy resource requirements; Java-based | v0.6.0 |

---

## Lower Priority

| Tool | Category | Notes |
|---|---|---|
| **BWA-FastAlign** | Alignment | Faster BWA-MEM2 fork. Marginal improvement — minimap2 is already fast enough for the default path, and BWA-MEM2 is already an alternative. |
| **Nextclade** | Viral | Not relevant for human germline WGS. |
| **iVar** | Viral | Amplicon-based variant calling for viral genomes. Out of scope. |
| **MACS3** | ChIP-seq | Peak calling for chromatin profiling. Not applicable. |
| **HTSeq / featureCounts** | RNA-seq | Read counting for gene expression. Not applicable. |
| **MultiQC alone** | Reporting | Only useful paired with FastQC or other QC tools that generate parseable logs. |

---

## Not a Fit

| Tool | Category | Reason |
|---|---|---|
| **STAR / HISAT2 / Salmon** | RNA-seq alignment/quantification | Pipeline is DNA-only |
| **DESeq2 / edgeR** | Differential expression | RNA-seq analysis |
| **Bismark / methylKit** | Methylation | Requires bisulfite sequencing data |
| **Kraken2 / MetaPhlAn** | Metagenomics | Not human WGS analysis |
| **SPAdes / MEGAHIT** | De novo assembly | Pipeline uses reference-based alignment |
| **Picard** | BAM utilities | Mostly redundant with samtools; CollectHsMetrics useful only for WES (v0.6.0 scope) |
| **Cromwell / WDL** | Workflow engine | Nextflow chosen for v0.5.0; adding a third engine creates maintenance burden |

---

## Notes

- This triage evaluates tools against the pipeline's current architecture (Docker containers, bash scripts, single-sample default). Tools that require fundamental architectural changes (e.g., GLNexus needing gVCF mode) are flagged.
- "Roadmap" column indicates the earliest version where the tool makes sense, not a commitment to include it.
- Tools already in the pipeline (minimap2, BWA-MEM2, DeepVariant, GATK, FreeBayes, Strelka2, Manta, Delly, VEP, PharmCAT, etc.) are not listed.
