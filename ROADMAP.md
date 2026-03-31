# Roadmap

Where the pipeline is headed. Items are grouped by priority and roughly ordered within each tier.

---

## v0.2.0 — Multi-sample & joint analysis

The biggest gap in v0.1.x is that every step runs on a single sample in isolation. v0.2.0 focuses on making the pipeline useful for families and cohorts.

- [ ] **Joint PCA with 1000 Genomes reference panel** — project sample PCs onto a reference PCA, replacing the current single-sample ancestry stub (step 26) with real population placement
- [ ] **Multi-sample SV merging** — merge Manta/Delly calls across 2+ samples (e.g., partners, parent-child) to identify shared and private structural variants
- [ ] **Carrier cross-check automation** — given two VCFs, automatically check shared autosomal recessive carrier status (currently manual; see `docs/multi-sample.md`)
- [ ] **PRS percentile estimation** — use a public reference cohort (e.g., UK Biobank summary stats) to convert raw PRS scores into approximate percentiles
- [ ] **Trio analysis support** — de novo variant calling and compound heterozygote phasing for parent-child trios

## v0.3.0 — Tool upgrades & expanded coverage

Upgrade pinned tools where the pipeline currently runs on older versions, and fill known coverage gaps.

- [ ] **ExpansionHunter v5.x** — upgrade from v2.5.5 to v5.x (new `--variant-catalog` flag, expanded locus catalog, better long-repeat estimation)
- [ ] **PharmCAT upgrade evaluation** — validate PharmCAT 3.x against the current v2.15.5 pin; update step 7 and step 27 together if JSON structure and preprocessor flags are stable
- [ ] **PCGR/CPSR data bundle refresh** — the current `grch38.20220203` bundle is from Feb 2022; upgrade to the latest available bundle
- [ ] **Long-read support** — add an optional ONT/PacBio alignment path (minimap2 already supports it) and long-read-aware SV calling (Sniffles2, cuteSV)
- [ ] **CYP2D6 star allele calling** — evaluate StellarPGx or Aldy as replacements for Cyrius, which returns None/None on most short-read WGS samples
- [ ] **Somatic variant calling** — optional tumor-only mode with Mutect2 for users who have matched tumor/normal WGS (rare but requested)

## v0.4.0 — Reporting & user experience

Make results more accessible to non-bioinformaticians.

- [ ] **Interactive HTML dashboard** — single-page HTML report combining all step outputs with collapsible sections, variant tables, PRS charts, and pharmacogenomics summaries
- [ ] **PDF clinical summary** — one-page printable summary designed to hand to a healthcare provider (PharmCAT results, ClinVar pathogenic hits, key carrier findings)
- [ ] **Automated database update script** — `scripts/update-databases.sh` that downloads the latest ClinVar, VEP cache, and PGS Catalog scoring files, with version tracking
- [ ] **Progress dashboard** — real-time terminal UI showing step status, elapsed time, and resource usage during `run-all.sh` execution
- [ ] **Conda/Bioconda alternative** — offer a non-Docker installation path for HPC environments where Docker is not available

## Ongoing maintenance

These are not versioned milestones but continuous responsibilities.

- **ClinVar monthly refresh** — re-download on the first Thursday of each month; re-run step 6 on at least one sample to verify
- **VEP cache refresh** — upgrade with each Ensembl release (~every 6 months; release 116 expected Apr 2026)
- **PGS Catalog quarterly check** — verify scoring file versions haven't changed; treat any change as a result-changing event
- **CPIC guideline monitoring** — recheck when PharmCAT bumps versions or when a drug-gene pair used in step 27 gets updated upstream
- **CI health** — keep ShellCheck, markdown-links, and contract-validation passing; add new contract checks as steps are added

## Not planned

These are out of scope for this pipeline and unlikely to be added.

- **Cloud execution** — the pipeline is designed for local consumer hardware. Cloud deployment adds complexity without clear benefit for the target audience.
- **GUI/web interface** — this is a CLI pipeline. If you want a GUI, consider wrapping it in Nextflow or Snakemake with their native UIs.
- **Clinical validation / CLIA compliance** — this pipeline is for personal exploration, not clinical diagnostics. It will never carry a clinical validation stamp.
- **Somatic-only (tumor without normal)** — requires a panel of normals and is methodologically fraught without matched germline.

---

Suggestions and contributions welcome via [issues](https://github.com/GeiserX/genomics-pipeline/issues) or [pull requests](https://github.com/GeiserX/genomics-pipeline/pulls).
