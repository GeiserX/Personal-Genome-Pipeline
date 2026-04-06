# Roadmap

Where the pipeline is headed. Items are grouped by priority and roughly ordered within each tier.

---

## v0.2.0 — Variant caller benchmarking & alternative tools ✅

Alternative callers, benchmarking infrastructure, and tool rationale documentation. Thanks to [@madmolecularman](https://github.com/madmolecularman) for pushing this direction.

- [x] **Alternative aligner: BWA-MEM2** (`scripts/02a-alignment-bwamem2.sh` → `aligned_bwamem2/`) — produces XS tags needed by Strelka2. All alternative caller scripts accept `ALIGN_DIR=aligned_bwamem2`
- [x] **Alternative SNP callers: GATK HaplotypeCaller + FreeBayes + Strelka2** — three optional callers alongside DeepVariant, each writing to isolated output directories (`vcf_gatk/`, `vcf_freebayes/`, `vcf_strelka2/`). Note: Strelka2 is a small variant caller (SNVs + indels ≤49bp), not an SV caller
- [x] **Concordance benchmarking script** (`scripts/benchmark-variants.sh`) — two modes: pairwise concordance (`bcftools isec` with PASS filter + normalization) and truth set benchmarking (`hap.py`). Auto-discovers all caller VCFs
- [x] **Alternative SV caller: TIDDIT** (`scripts/04a-tiddit.sh` → `sv_tiddit/`) — excels at large inversions and translocations; auto-detects BWA index for assembly mode
- [x] **Documented tool rationale** (`docs/tool-rationale.md`) — per-step rationale with references to benchmarking data and decision matrices

## v0.3.0 — Tool upgrades, QC, & expanded coverage

Upgrade pinned tools, add pre-alignment QC ([#14](https://github.com/GeiserX/Personal-Genome-Pipeline/issues/14)), fill coverage gaps, and add new sequencing platform support. Thanks to [@madmolecularman](https://github.com/madmolecularman) for driving the QC discussion.

- [x] **fastp QC + adapter trimming** (`scripts/01b-fastp-qc.sh` → `fastq_trimmed/`) — pre-alignment step: adapter removal, quality trimming, polyG tail removal, JSON+HTML reports. Default on, skippable with `SKIP_TRIM=true`. Addresses [#14](https://github.com/GeiserX/Personal-Genome-Pipeline/issues/14)
- [x] **mosdepth coverage statistics** (`scripts/16b-mosdepth.sh`) — fast per-base and per-region depth from BAM; coverage distributions, threshold reports, and WES on-target rate
- [x] **MultiQC aggregation** (`scripts/28-multiqc.sh`) — single HTML report combining fastp, samtools flagstat, mosdepth, and other QC outputs
- [x] **ExpansionHunter upgrade to v5.0.0** — new `--variant-catalog` flag replacing `--variant-catalog-format`, biocontainers image replacing deprecated weisburd image
- [x] **PharmCAT upgrade to 3.2.0** — preprocessor renamed (no `.py`), explicit reporter flags (`-reporterJson -reporterHtml`), `wildtypeAllele` → `referenceAllele` in JSON, NAT2 calling added. Step 7 and step 27 updated together
- [x] **PCGR/CPSR upgrade to 2.2.5** — upgraded from 1.4.1 to 2.2.5 with new ref data bundle (`20250314`), separate VEP cache mount, and completely rewritten CLI
- [x] **Octopus variant caller** (`scripts/03d-octopus.sh` → `vcf_octopus/`) — haplotype-aware Bayesian caller as a 5th benchmarking alternative. Auto-discovered by `benchmark-variants.sh`
- [x] **GRIDSS structural variant caller** (`scripts/04b-gridss.sh` → `sv_gridss/`) — assembly-based SV caller for complex rearrangements; strengthens SV consensus alongside Manta/Delly. Requires BWA index and 32GB RAM
- [x] **CYP2D6 star allele calling** — evaluated Aldy v4.8.3 (best CYP2D6 SV caller per Twesigomwe 2020), StellarPGx (broken Docker), and BCyrius (no public repo). Aldy documented as recommended optional replacement for Cyrius in `docs/21-cyrius.md`. Note: Aldy uses an academic-only license (IURTC) incompatible with GPL-3.0, so it cannot be a required dependency. pypgx (v0.4.0) remains the long-term GPL-compatible replacement
- [x] **Long-read support** (`scripts/02b-alignment-longread.sh`, `scripts/03e-clair3.sh`, `scripts/04c-sniffles2.sh`) — ONT and PacBio HiFi alignment (minimap2 `map-ont`/`map-hifi`), Clair3 v2.0.0 variant calling, Sniffles2 SV calling. Comprehensive guide in `docs/long-read-guide.md`
- [x] **Whole exome sequencing (WES) entry path** — comprehensive guide in `docs/wes-guide.md` covering per-step compatibility, capture BED files, `DATA_TYPE=WES` env var, coverage QC metrics, and limitations. Thanks to [@madmolecularman](https://github.com/madmolecularman) for domain expertise here
- [x] **Somatic variant calling** (`scripts/29-mutect2-somatic.sh`) — [EXPERIMENTAL] tumor-only Mutect2 mode with gnomAD germline resource and Panel of Normals filtering. Marked experimental due to high false positive rate without matched normal

## v0.4.0 — Expanded annotation & clinical interpretation

Current annotation is VEP-only with basic ClinVar screening. Clinical-grade interpretation benefits from deeper pathogenicity scoring and structured querying.

- [ ] **CADD scores** — integrate Combined Annotation Dependent Depletion scores for all variants (coding + non-coding), the most widely used deleteriousness metric
- [ ] **SpliceAI** — deep learning splice-site variant prediction. Catches pathogenic intronic variants that VEP's rule-based splice prediction misses
- [ ] **REVEL scores** — ensemble pathogenicity scoring for missense variants, combining 13 individual tools. Recommended by ClinGen for missense variant classification
- [ ] **AlphaMissense** — DeepMind's protein-structure-informed missense classifier. Complements REVEL with structural context
- [ ] **gnomAD v4 constraint metrics** — per-gene pLI, LOEUF, and missense Z-scores. Essential for interpreting novel variants in loss-of-function-intolerant genes
- [ ] **vcfanno annotation engine** — add arbitrary annotation tracks (gnomAD, CADD, SpliceAI, custom BEDs) to VCFs via TOML config; faster than VEP for bulk annotation, complements rather than replaces it
- [ ] **Variant database with inheritance queries** — load annotated VCFs into a queryable store (GEMINI or modern successor) supporting inheritance model filtering: de novo, compound het, X-linked, autosomal recessive/dominant
- [ ] **pypgx alongside PharmCAT** — broader PGx star allele calling including CYP2D6 structural variation detection from WGS reads, filling the gap left by Cyrius

## v0.5.0 — Workflow engine integration

The 27 bash scripts work but lack built-in parallelism, resume-on-failure, and HPC portability. The [nf-core](https://nf-co.re/) ecosystem (147 community pipelines including [sarek](https://nf-co.re/sarek) with 15 variant callers and [raredisease](https://github.com/nf-core/raredisease) for clinical genomics) demonstrates the community standard.

- [ ] **Nextflow DSL2 wrapper** — convert the pipeline into a Nextflow workflow with channels and processes, preserving the current Docker-based execution model
- [ ] **nf-core module compatibility** — use [nf-core/modules](https://github.com/nf-core/modules) where they exist (BWA-MEM2, DeepVariant, VEP, Manta, bcftools) for community-maintained containers and automated testing
- [ ] **Snakemake alternative** — optional Snakemake wrapper for HPC environments that prefer it over Nextflow
- [ ] This unlocks: automatic parallelism via DAG-based step ordering, resume on failure, Singularity/Apptainer for HPC clusters, and optional cloud portability

## v0.6.0 — Multi-sample & joint analysis

Every step currently runs on a single sample in isolation. v0.6.0 focuses on making the pipeline useful for families and cohorts.

- [ ] **Joint PCA with 1000 Genomes reference panel** — project sample PCs onto a reference PCA, replacing the current single-sample ancestry stub (step 26) with real population placement
- [ ] **Multi-sample SV merging** — merge Manta/Delly calls across 2+ samples (e.g., partners, parent-child) to identify shared and private structural variants
- [ ] **Carrier cross-check automation** — given two VCFs, automatically check shared autosomal recessive carrier status (currently manual; see `docs/multi-sample.md`)
- [ ] **PRS percentile estimation** — use a public reference cohort (e.g., UK Biobank summary stats) to convert raw PRS scores into approximate percentiles
- [ ] **Somalier sample identity QC** — ultra-fast relatedness and sample-swap detection from BAM/VCF; replaces ad-hoc sex-check with proper identity QC for multi-sample runs
- [ ] **GLNexus joint genotyping** — merge per-sample gVCFs into joint-called cohort VCFs; requires switching DeepVariant to `--output_gvcf` mode
- [ ] **Trio analysis support** — de novo variant calling and compound heterozygote phasing for parent-child trios, with GEMINI-style inheritance model queries (de novo, compound het, X-linked recessive, autosomal recessive)

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

Suggestions and contributions welcome via [issues](https://github.com/GeiserX/Personal-Genome-Pipeline/issues) or [pull requests](https://github.com/GeiserX/Personal-Genome-Pipeline/pulls).
