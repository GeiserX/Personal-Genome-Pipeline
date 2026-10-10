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

## v0.3.0 — Tool upgrades, QC, & expanded coverage ✅

Upgrade pinned tools, add pre-alignment QC ([#14](https://github.com/GeiserX/Personal-Genome-Pipeline/issues/14)), fill coverage gaps, and add new sequencing platform support. Thanks to [@madmolecularman](https://github.com/madmolecularman) for driving the QC discussion.

- [x] **fastp QC + adapter trimming** (`scripts/01b-fastp-qc.sh` → `fastq_trimmed/`) — pre-alignment step: adapter removal, quality trimming, polyG tail removal, JSON+HTML reports. Default on, skippable with `SKIP_TRIM=true`. Addresses [#14](https://github.com/GeiserX/Personal-Genome-Pipeline/issues/14)
- [x] **mosdepth coverage statistics** (`scripts/16b-mosdepth.sh`) — fast per-base and per-region depth from BAM; coverage distributions, threshold reports, and WES on-target rate
- [x] **MultiQC aggregation** (`scripts/28-multiqc.sh`) — single HTML report combining fastp, samtools flagstat, mosdepth, and other QC outputs
- [x] **ExpansionHunter upgrade to v5.0.0** — new `--variant-catalog` flag replacing `--variant-catalog-format`, biocontainers image replacing deprecated weisburd image
- [x] **PharmCAT upgrade to 3.2.0** — preprocessor renamed (no `.py`), explicit reporter flags (`-reporterJson -reporterHtml`), `wildtypeAllele` → `referenceAllele` in JSON, NAT2 calling added. Step 7 and step 27 updated together
- [x] **PCGR/CPSR upgrade to 2.2.5** — upgraded from 1.4.1 to 2.2.5 with new ref data bundle (`20250314`), separate VEP cache mount, and completely rewritten CLI
- [x] **Octopus variant caller** (`scripts/03d-octopus.sh` → `vcf_octopus/`) — haplotype-aware Bayesian caller as a 5th benchmarking alternative. Auto-discovered by `benchmark-variants.sh`
- [x] **GRIDSS structural variant caller** (`scripts/04b-gridss.sh` → `sv_gridss/`) — assembly-based SV caller for complex rearrangements; strengthens SV consensus alongside Manta/Delly. Requires BWA index and 32GB RAM
- [x] **CYP2D6 star allele calling** — evaluated Aldy v4.8.3 (best CYP2D6 SV caller per Twesigomwe 2020), StellarPGx (broken Docker), and BCyrius (no public repo). Aldy documented as recommended optional replacement for Cyrius in `docs/21-cyrius.md`. Note: Aldy uses an academic-only license (IURTC) incompatible with GPL-3.0, so it cannot be a required dependency. pypgx (step 32) is the GPL-compatible alternative included in the pipeline
- [x] **Long-read support** (`scripts/02b-alignment-longread.sh`, `scripts/03e-clair3.sh`, `scripts/04c-sniffles2.sh`) — ONT and PacBio HiFi alignment (minimap2 `map-ont`/`map-hifi`), Clair3 variant calling, Sniffles2 SV calling. Comprehensive guide in `docs/long-read-guide.md`
- [x] **Whole exome sequencing (WES) guide** — `docs/wes-guide.md` covers per-step compatibility, capture BED files, coverage QC metrics and limitations. WES data runs step by step; no script has a WES mode. Thanks to [@madmolecularman](https://github.com/madmolecularman) for domain expertise here
- [x] **Somatic variant calling** (`scripts/29-mutect2-somatic.sh`) — [EXPERIMENTAL] tumor-only Mutect2 mode with gnomAD germline resource and Panel of Normals filtering. Marked experimental due to high false positive rate without matched normal

## v0.4.0 — Expanded annotation & clinical interpretation ✅

Deep pathogenicity scoring, structured variant querying, and broader pharmacogenomics. All new annotation tracks are optional — scripts detect which databases are present and degrade gracefully.

- [x] **CADD scores** — Combined Annotation Dependent Depletion scores for all variants via vcfanno. Pre-scored whole-genome SNVs (~81.5 GB) + gnomAD indels (~1.2 GB). PHRED >= 20 flagged as clinically interesting
- [x] **SpliceAI** — deep learning splice-site variant predictions via vcfanno. Pre-scored masked files (~91 GB). Delta score >= 0.2 flagged for cryptic splice variants
- [x] **REVEL scores** — ensemble missense pathogenicity scoring via vcfanno (~0.6 GB). ClinGen-calibrated thresholds: >= 0.644 (PP3_Supporting), >= 0.773 (PP3_Moderate), >= 0.932 (PP3_Strong)
- [x] **AlphaMissense** — DeepMind's protein-structure-informed missense classifier via vcfanno (~613 MB). Class boundaries (AlphaMissense's own, not ACMG levels): < 0.34 likely benign, > 0.564 likely pathogenic
- [x] **gnomAD v4 constraint metrics** — per-gene pLI, LOEUF, and missense Z-scores (~91 MB). Integrated into clinical filter summary TSV and slivar output
- [x] **vcfanno annotation engine** (`scripts/30-vcfanno.sh`) — adds CADD, SpliceAI, REVEL, AlphaMissense to VEP VCFs via TOML config in a single pass. Handles CADD chr prefix mismatch with two-pass approach
- [x] **Variant prioritization with inheritance queries** (`scripts/31-slivar.sh`) — slivar (GEMINI successor) for streaming VCF filtering with JS expressions. Rare HIGH/MODERATE variants, ClinVar pathogenic, compound het detection, gene constraint enrichment
- [x] **pypgx alongside PharmCAT** (`scripts/32-pypgx.sh`) — 23-gene curated star allele calling including CYP2D6 structural variation from BAM read depth. Cross-validates with PharmCAT on shared genes

## v0.5.0 — Nextflow workflow engine ✅

The bash scripts work but lack built-in parallelism, resume-on-failure, and HPC portability. v0.5.0 adds a [Nextflow](https://www.nextflow.io/) DSL2 execution path alongside the existing bash scripts (which remain first-class).

### Why Nextflow over Snakemake?

Both are mature workflow engines. We chose Nextflow because:

- **nf-core ecosystem**: 147 community pipelines including [sarek](https://nf-co.re/sarek) (WGS variant calling) and [raredisease](https://github.com/nf-core/raredisease) (clinical genomics). Sarek's output is our primary input — channel compatibility matters.
- **nf-core/modules**: 1000+ reusable modules. We can both use existing modules and contribute novel ones (PharmCAT, pypgx, slivar) under MIT.
- **Container-first design**: Nextflow's `container` directive maps directly to our Docker-based architecture. Singularity support is automatic.
- **Resume**: Content-hash caching is more robust than file-existence checks.
- **Industry momentum**: Seqera/Nextflow has commercial backing; major sequencing centers standardize on Nextflow.

Snakemake's Python DSL and HPC scheduler integration are genuine strengths, but the nf-core ecosystem size and sarek compatibility are decisive.

### Scope

v0.5.0 started as post-processing: alignment and variant calling were left to nf-core/sarek, and the Nextflow pipeline took sarek's VCF and BAM. That is no longer the scope. The Nextflow pipeline now starts from FASTQ as well (fastp, minimap2 with read groups and duplicate marking, the indexcov sex check, DeepVariant with a gVCF), and a samplesheet row can also start from a BAM, a CRAM or a VCF, sarek's included. Nextflow is the pipeline; the numbered scripts are single steps that share its images and helpers ([docs/nextflow.md](docs/nextflow.md)).

### Delivery

v0.5.0 shipped 27 modules in 6 workflows, for interpretation after calling. CI ran their stub-testable subset; the tools that need a database were checked by hand.

Today the Nextflow path has 39 modules in 7 workflows and runs from FASTQ to the report. The e2e job runs the real tools of most steps on a slice of HG002 and compares the pipeline with the bash steps on the same reads ([docs/testing.md](docs/testing.md)). The Nextflow ExpansionHunter module is only stub-run. VEP, CPSR and AnnotSV need databases CI does not download, and GRIDSS needs a 31 GB Java heap, so CI does not run them. See [docs/nextflow.md](docs/nextflow.md) for known limitations.

- [x] **PR #17 — Full Nextflow pipeline** (v0.5.0): All 6 workflows (PGX, ANNOTATION, CLINICAL, BAM_ANALYSIS, SV, REPORTING) with 27 modules, `--tools` gating, stub CI, Docker + Singularity profiles

### Parallel track: nf-core module contributions

PharmCAT, pypgx, and slivar modules will be contributed to [nf-core/modules](https://github.com/nf-core/modules) under MIT license — independent of the pipeline's GPL-3.0 license. Once merged, these modules will be available to all nf-core pipelines.

### Bash scripts

`run-all.sh` now launches the Nextflow pipeline, and new analysis steps are Nextflow-first. The numbered scripts in `scripts/` stay as documented single-step commands that share the pipeline's images and helpers, and get bug fixes and tool version bumps. The steps that stay bash-only are listed with their reasons in [docs/nextflow.md](docs/nextflow.md#bash-vs-nextflow-parity).

## v0.6.0 to v0.8.2 — Hardening and tool upgrades ✅

Released on 2026-07-01. These releases went to hardening and upgrades instead of the multi-sample and reporting work first planned for v0.6.0 and v0.7.0, which moved to the planned sections below.

- [x] **v0.6.0** — Stranger STR annotation (step 9b), a CPIC step that fails loudly instead of silently, safe tool bumps, Nextflow 25.10 support, GitHub-hosted CI runners
- [x] **v0.7.0** — DeepVariant 1.10, VEP 116, Delly 2.1, and ACMG secondary findings in CPSR (step 17)
- [x] **v0.7.1** — the last two floating `:latest` images pinned by digest
- [x] **v0.8.0** — CNV calling moved from CNVnator to CNVpytor 1.3.2 (step 18)
- [x] **v0.8.1, v0.8.2** — CNVpytor restricted to the main chromosomes; the pypgx image pinned to 0.26.0 to match its bundle

## Planned — Multi-sample & joint analysis

Every step currently runs on a single sample in isolation. This work would make the pipeline useful for families and cohorts; no release is assigned to it yet.

- [x] **Projection onto the 1000 Genomes reference panel** — step 26 projects the sample onto pgsc_calc's panel and names the most similar population, replacing the single-sample PCA
- [ ] **Multi-sample SV merging** — merge Manta/Delly calls across 2+ samples (e.g., partners, parent-child) to identify shared and private structural variants
- [ ] **Carrier cross-check automation** — given two VCFs, automatically check shared autosomal recessive carrier status (currently manual; see `docs/multi-sample.md`)
- [x] **PRS percentile estimation** — step 25 runs pgsc_calc, which with the ancestry panel reports each score as a percentile among the most similar reference group
- [x] **Somalier sample identity QC** — step 33 (`sample_qc` in Nextflow): somalier's sex from the reads checked against the declared sex, the relatedness of every pair of samples in a run, and VerifyBamID2's contamination estimate
- [ ] **GLNexus joint genotyping** — merge per-sample gVCFs into joint-called cohort VCFs; requires switching DeepVariant to `--output_gvcf` mode
- [ ] **Trio analysis support** — de novo variant calling and compound heterozygote phasing for parent-child trios, with slivar inheritance model queries (de novo, compound het, X-linked recessive, autosomal recessive)

## Planned — Reporting & user experience

Make results more accessible to non-bioinformaticians.

- [ ] **Interactive HTML dashboard** — single-page HTML report combining all step outputs with collapsible sections, variant tables, PRS charts, and pharmacogenomics summaries
- [ ] **PDF clinical summary** — one-page printable summary designed to hand to a healthcare provider (PharmCAT results, ClinVar pathogenic hits, key carrier findings)
- [ ] **Database refresh** — `setup.sh --refresh clinvar` to fetch the latest ClinVar in place, and a monthly issue that lists which databases and images have a newer release
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
- **VEP cache refresh** — upgrade with each Ensembl release (~every 6 months); release 116 is the one pinned now
- **PGS Catalog quarterly check** — verify scoring file versions haven't changed; treat any change as a result-changing event
- **CPIC guideline monitoring** — recheck when PharmCAT bumps versions or when a drug-gene pair used in step 27 gets updated upstream
- **CI health** — keep ShellCheck, markdown-links, and contract-validation passing; add new contract checks as steps are added

## Not planned

These are out of scope for this pipeline and unlikely to be added.

- **Cloud-only execution** — the pipeline is local-first by design. v1.0.0 envisions optional cloud portability via workflow engines, but the default path will always be local consumer hardware.
- **Clinical validation / CLIA compliance** — this pipeline is for personal exploration, not clinical diagnostics. It will never carry a clinical validation stamp.

---

Suggestions and contributions welcome via [issues](https://github.com/GeiserX/Personal-Genome-Pipeline/issues) or [pull requests](https://github.com/GeiserX/Personal-Genome-Pipeline/pulls).
