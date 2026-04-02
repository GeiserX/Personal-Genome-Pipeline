# Roadmap

Where the pipeline is headed. Items are grouped by priority and roughly ordered within each tier.

---

## v0.2.0 — Variant caller benchmarking & alternative tools

The pipeline currently ships one tool per step with no way to compare. Bioinformaticians routinely benchmark callers against each other because no single tool is universally superior ([PLOS ONE comparison](https://doi.org/10.1371/journal.pone.0339891), [UMCCR BWA vs minimap2](https://umccr.org/blog/bwa-mem-vs-minimap2/)). Thanks to [@madmolecularman](https://github.com/madmolecularman) for pushing this direction.

- [ ] **Alternative aligner: BWA-MEM2** — add BWA-MEM2 alongside minimap2. BWA-MEM2 produces different SAM tags (XS suboptimal alignment score) that some callers depend on — Strelka2 in particular has reduced SNP precision with minimap2 alignments due to missing XS tags and doubled AS scores
- [ ] **Alternative SNP callers: GATK HaplotypeCaller + FreeBayes** — add as optional callers alongside DeepVariant. DeepVariant leads in precision and F1; GATK balances precision/recall; FreeBayes maximizes sensitivity at the cost of higher false positives. Different callers catch different variants — no single caller finds everything
- [ ] **Concordance benchmarking script** — given a GIAB truth set (HG002 or NA12878) and one or more caller VCFs, run Illumina's [hap.py](https://github.com/Illumina/hap.py) or RTG vcfeval to produce precision/recall/F1 stratified by variant type (SNP vs indel) and genomic region. Output a comparison table and optional Venn diagram of caller-unique and shared calls
- [ ] **Alternative SV callers: TIDDIT + Strelka2** — add alongside Manta/Delly for broader SV sensitivity. TIDDIT excels at large inversions and translocations; Strelka2 germline mode catches small indels that Manta misses
- [ ] **Documented tool rationale** — for each step, document why the default was chosen with references to benchmarking data, so users can make informed decisions about which callers to run

## v0.3.0 — Multi-sample & joint analysis

The biggest gap in v0.1.x is that every step runs on a single sample in isolation. v0.3.0 focuses on making the pipeline useful for families and cohorts.

- [ ] **Joint PCA with 1000 Genomes reference panel** — project sample PCs onto a reference PCA, replacing the current single-sample ancestry stub (step 26) with real population placement
- [ ] **Multi-sample SV merging** — merge Manta/Delly calls across 2+ samples (e.g., partners, parent-child) to identify shared and private structural variants
- [ ] **Carrier cross-check automation** — given two VCFs, automatically check shared autosomal recessive carrier status (currently manual; see `docs/multi-sample.md`)
- [ ] **PRS percentile estimation** — use a public reference cohort (e.g., UK Biobank summary stats) to convert raw PRS scores into approximate percentiles
- [ ] **Trio analysis support** — de novo variant calling and compound heterozygote phasing for parent-child trios, with GEMINI-style inheritance model queries (de novo, compound het, X-linked recessive, autosomal recessive)

## v0.4.0 — Expanded annotation & clinical interpretation

Current annotation is VEP-only with basic ClinVar screening. Clinical-grade interpretation benefits from deeper pathogenicity scoring and structured querying.

- [ ] **CADD scores** — integrate Combined Annotation Dependent Depletion scores for all variants (coding + non-coding), the most widely used deleteriousness metric
- [ ] **SpliceAI** — deep learning splice-site variant prediction. Catches pathogenic intronic variants that VEP's rule-based splice prediction misses
- [ ] **REVEL scores** — ensemble pathogenicity scoring for missense variants, combining 13 individual tools. Recommended by ClinGen for missense variant classification
- [ ] **AlphaMissense** — DeepMind's protein-structure-informed missense classifier. Complements REVEL with structural context
- [ ] **gnomAD v4 constraint metrics** — per-gene pLI, LOEUF, and missense Z-scores. Essential for interpreting novel variants in loss-of-function-intolerant genes
- [ ] **Variant database with inheritance queries** — load annotated VCFs into a queryable store (GEMINI or modern successor) supporting inheritance model filtering: de novo, compound het, X-linked, autosomal recessive/dominant
- [ ] **pypgx alongside PharmCAT** — broader PGx star allele calling including CYP2D6 structural variation detection from WGS reads, filling the gap left by Cyrius

## v0.5.0 — Workflow engine integration

The 27 bash scripts work but lack built-in parallelism, resume-on-failure, and HPC portability. The [nf-core](https://nf-co.re/) ecosystem (147 community pipelines including [sarek](https://nf-co.re/sarek) with 15 variant callers and [raredisease](https://github.com/nf-core/raredisease) for clinical genomics) demonstrates the community standard.

- [ ] **Nextflow DSL2 wrapper** — convert the pipeline into a Nextflow workflow with channels and processes, preserving the current Docker-based execution model
- [ ] **nf-core module compatibility** — use [nf-core/modules](https://github.com/nf-core/modules) where they exist (BWA-MEM2, DeepVariant, VEP, Manta, bcftools) for community-maintained containers and automated testing
- [ ] **Snakemake alternative** — optional Snakemake wrapper for HPC environments that prefer it over Nextflow
- [ ] This unlocks: automatic parallelism via DAG-based step ordering, resume on failure, Singularity/Apptainer for HPC clusters, and optional cloud portability

## v0.6.0 — Tool upgrades & expanded coverage

Upgrade pinned tools where the pipeline currently runs on older versions, and fill known coverage gaps.

- [ ] **ExpansionHunter v5.x** — upgrade from v2.5.5 to v5.x (new `--variant-catalog` flag, expanded locus catalog, better long-repeat estimation)
- [ ] **PharmCAT upgrade evaluation** — validate PharmCAT 3.x against the current v2.15.5 pin; update step 7 and step 27 together if JSON structure and preprocessor flags are stable
- [ ] **PCGR/CPSR data bundle refresh** — the current `grch38.20220203` bundle is from Feb 2022; upgrade to the latest available bundle
- [ ] **Long-read support** — add an optional ONT/PacBio alignment path (minimap2 `map-ont`, pbmm2 for HiFi), long-read-aware variant calling (Clair3 for SNPs, Sniffles2 and cuteSV for SVs), and [pb-StarPhase](https://github.com/PacificBiosciences/pb-StarPhase) for long-read pharmacogenomics
- [ ] **Whole exome sequencing (WES) entry path** — support for targeted/exome BAMs with coverage thresholds, on-target rate QC (Picard CollectHsMetrics), and adjusted variant caller parameters. WES is where most rare disease diagnosis and cancer research is moving. Thanks to [@madmolecularman](https://github.com/madmolecularman) for domain expertise here
- [ ] **CYP2D6 star allele calling** — evaluate StellarPGx or Aldy as replacements for Cyrius, which returns None/None on most short-read WGS samples
- [ ] **Somatic variant calling** — optional tumor-only mode with Mutect2 for users who have matched tumor/normal WGS (rare but requested)

## v0.7.0 — Reporting & user experience

Make results more accessible to non-bioinformaticians.

- [ ] **Interactive HTML dashboard** — single-page HTML report combining all step outputs with collapsible sections, variant tables, PRS charts, and pharmacogenomics summaries
- [ ] **PDF clinical summary** — one-page printable summary designed to hand to a healthcare provider (PharmCAT results, ClinVar pathogenic hits, key carrier findings)
- [ ] **Automated database update script** — `scripts/update-databases.sh` that downloads the latest ClinVar, VEP cache, and PGS Catalog scoring files, with version tracking
- [ ] **Progress dashboard** — real-time terminal UI showing step status, elapsed time, and resource usage during `run-all.sh` execution
- [ ] **Conda/Bioconda alternative** — offer a non-Docker installation path for HPC environments where Docker is not available

## v1.0.0 — Local-first distributed platform

The long-term vision: a self-hostable, open-source health data platform — the "Home Assistant of personal genomics."

- [ ] **Local worker agent** — lightweight daemon that runs the pipeline on the user's own hardware, reporting progress to a central coordinator. Users with powerful machines run their own analysis; users without can opt into a shared compute pool
- [ ] **Central web UI** — result visualization, step management, multi-sample comparison dashboard. No raw genomic data stored centrally — only aggregate results and metadata
- [ ] **Data sovereignty by design** — raw FASTQ/BAM/VCF never leaves the user's machine. Only aggregate, non-re-identifiable results (PRS scores, PGx diplotypes, carrier status) are optionally shared with the coordinator for cross-sample features
- [ ] **Health integration layer** — combine genomic findings with clinical history, lab results, and supplement plans in a unified personal health record. Structured data import from common lab formats (HL7 FHIR, PDF extraction)
- [ ] **Consumer DNA one-click import** — streamlined import from 23andMe, AncestryDNA, MyHeritage, and other consumer genotyping platforms, building on the existing `chip-to-vcf.sh` converter
- [ ] Open-source everything, self-hostable, no vendor lock-in

## Ongoing maintenance

These are not versioned milestones but continuous responsibilities.

- **ClinVar monthly refresh** — re-download on the first Thursday of each month; re-run step 6 on at least one sample to verify
- **VEP cache refresh** — upgrade with each Ensembl release (~every 6 months; release 116 expected Apr 2026)
- **PGS Catalog quarterly check** — verify scoring file versions haven't changed; treat any change as a result-changing event
- **CPIC guideline monitoring** — recheck when PharmCAT bumps versions or when a drug-gene pair used in step 27 gets updated upstream
- **CI health** — keep ShellCheck, markdown-links, and contract-validation passing; add new contract checks as steps are added

## Not planned

These are out of scope for this pipeline and unlikely to be added.

- **Cloud-only execution** — the pipeline is local-first by design. v1.0.0 envisions optional cloud portability via workflow engines, but the default path will always be local consumer hardware.
- **Clinical validation / CLIA compliance** — this pipeline is for personal exploration, not clinical diagnostics. It will never carry a clinical validation stamp.
- **Somatic-only (tumor without normal)** — requires a panel of normals and is methodologically fraught without matched germline.

---

Suggestions and contributions welcome via [issues](https://github.com/GeiserX/genomics-pipeline/issues) or [pull requests](https://github.com/GeiserX/genomics-pipeline/pulls).
