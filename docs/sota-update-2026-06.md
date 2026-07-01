# SOTA Update Roadmap — mid-2026

A point-in-time review of every tool/container/database against its latest upstream release, with the recommended action. Versions confirmed from each project's GitHub `releases/latest` or vendor page.

> Status: **planning roadmap.** The strict-parser config fix already landed on `main` (#30/#31). The bumps below — especially the variant-caller and VEP-cache changes — require a full re-run on a known sample to validate before merging (see [`lessons-learned.md`](lessons-learned.md) and the revalidation checklist in `CLAUDE.md`). None of them are applied here; this doc is the backlog.

## Priority actions
1. **Nextflow strict-syntax** — already addressed on `main` (#30/#31; `nextflow.config` no longer uses a top-level `def`). For NF 26.x also migrate `conf/base.config` `check_max()` → `process.resourceLimits`, and the `vcfanno` optional-input scope. The pipeline currently runs on **NF 25.10.4**.
2. **DeepVariant 1.6.0 → 1.10.0** — accuracy gains, pangenome-aware reassembly, native long-read phasing. Re-calls variants → full re-run + concordance check vs the previous VCF.
3. **VEP 112 → 116** + cache 112 → 116 (must match) + **dbNSFP 5.3.1** (one source for REVEL + AlphaMissense + CADD + MetaRNN).
4. **samtools/bcftools 1.20/1.21 → 1.23.1** — note the **CRAM 3.0 → 3.1 default flip at 1.22** (write `--output-fmt cram,version=3.0` if older readers must consume the CRAM).
5. **Add**: Cyrius (CYP2D6, already wired), AlphaMissense (VEP plugin), pgsc_calc (PRS), ACMG SF v3.3 secondary-findings (84 genes).
6. **Refresh DBs** (operational — user-supplied reference data, not pinned in-repo): ClinVar (monthly), gnomAD v4.1, CADD v1.7, AlphaMissense hg38, PGS Catalog, IPD-IMGT/HLA 3.60+.

## Tool versions

| Tool | Current | Latest | Action | Notes |
|---|---|---|---|---|
| minimap2 | 2.28 | 2.31 | bump | drop-in |
| samtools | 1.20 | 1.23.1 | bump | CRAM 3.0→3.1 default at 1.22; keep htslib/samtools/bcftools aligned |
| bcftools | 1.21 | 1.23.1 | bump | fixes silent output truncation |
| DeepVariant | 1.6.0 | 1.10.0 | bump | re-call; standard WGS path unchanged |
| fastp | 1.3.1 | 1.3.6 | bump | BGZF multithread hang fixes |
| Manta | 1.6.0 | 1.6.0 (EOL, archived Oct 2025) | keep | still de-facto short-read germline SV caller; no successor |
| Delly | 1.7.3 | 2.1.0 | bump (major) | short-read PE/SR path backward-compatible; re-test SV step |
| CNVnator | 0.4.1 | — | **replace → CNVpytor 1.3.1** | same lab, maintained Python rewrite, CRAM+BAF |
| GRIDSS | 2.13.2 | 2.13.2 (2022) | keep/optional | unmaintained; droppable for single-genome runs |
| AnnotSV | 3.4.4 | 3.5.10 | bump | refresh bundled annotations |
| VEP | release_112 | release_116 | bump | cache must match; dbNSFP 5.3.1 |
| PCGR/CPSR | 2.2.5 | 2.3.0 | bump | re-download data bundle; CPSR ACMG-SF mode |
| PharmCAT | 3.2.0 | 3.2.0 | **keep (latest)** | feed Cyrius CYP2D6 as outside-call |
| TelomereHunter | `latest` | 1.1.0 / TH2 | **pin 1.1.0** | stop using `latest`; eval TelomereHunter2 |
| haplogrep3 | `latest` | v3.3.2 | **pin v3.3.2** | |
| T1K | 1.0.9 | 1.0.9 | keep | refresh IPD-IMGT/HLA ref (3.60+) |
| ExpansionHunter | 5.0.0 | 5.0.0 | keep | add `stranger` for annotation |
| goleft | 0.2.4 | 0.2.6 | bump | drop-in |
| mosdepth | 0.3.13 | 0.3.14 | bump | drop-in |
| GATK | 4.6.1.0 | 4.6.2.0 | bump | |
| Picard | 3.4.0 | 3.4.0 | keep | |
| PLINK2 | 2.00a5.10 | 2.00a7.1 | bump | pin by build |
| MultiQC | 1.33 | 1.35 | bump | min Python 3.9 |
| vcfanno | 0.3.7 | 0.3.9 | bump | |
| slivar | 0.3.3 | 0.3.4 | bump | |
| pypgx | 0.26.0 | 0.27.0 | bump | PharmCAT stays primary |
| Clair3 | 2.0.0 | 2.0.2 | bump | long-read path only |
| Sniffles | 2.4 | 2.8.0 | bump | long-read path only |

> Biocontainer tags carry a build suffix (e.g. `…1.23.1--h96c455f_0`); the **version** is fixed as above — pick the newest `_N` at pin time. Keep `versions.env` and each module's `container` tag in sync (CI enforces this).

## Database refreshes
- **ClinVar** — latest weekly `vcf_GRCh38/clinvar.vcf.gz` (pin the dated file). Operational refresh (user-supplied reference data; not pinned in this repo).
- **VEP cache** → release 116 (`homo_sapiens_vep_116_GRCh38.tar.gz`). Use `wget -c` (resumable; the 26 GB download cannot be resumed by VEP's `INSTALL.pl`).
- **dbNSFP 5.3.1** — single VEP `--plugin dbNSFP` source for REVEL + AlphaMissense + CADD + MetaRNN (lets you retire standalone annotators).
- **gnomAD v4.1/v4.1.1** GRCh38 (constraint recalculated, AN bug fixed). No v5 yet.
- **CADD v1.7** GRCh38 (ESM-1v + regulatory CNN + Zoonomia).
- **AlphaMissense hg38** — `AlphaMissense_hg38.tsv.gz` → `tabix -s1 -b2 -e2 -S1`; VEP `--plugin AlphaMissense,file=…`.
- **PGS Catalog** — via `pgsc_calc` (don't hand-roll scoring files).

## Steps to add (ranked)
1. **Cyrius** CYP2D6 star-allele caller (CNV/hybrid alleles PharmCAT misses) — already wired into the default tool set; feed its diplotype into PharmCAT as an outside-call. *(pypgx resolves `*5/*5` deletions where Cyrius can return None and PharmCAT reports No Result — keep all three and reconcile; see lessons-learned.)*
2. **AlphaMissense** via VEP plugin — easy, high value.
3. **pgsc_calc** (Nextflow, NF-26 compatible) — SOTA polygenic scoring.
4. **ACMG SF v3.3** (2025, 84 genes) via CPSR secondary-findings mode (PCGR 2.3.0).
5. Consolidate missense annotation on **dbNSFP 5.3.1**.
6. **stranger** for repeat-expansion annotation (with ExpansionHunter).
7. **mtDNA heteroplasmy** (mutserve or GATK Mutect2 mito-mode) + VEP gnomADMT.

Long-read-only (skip for short-read Illumina WGS — no methylation signal in the data): modkit, vg-giraffe pangenome.

## Drop / replace
- CNVnator → **CNVpytor**.
- Manta → keep but **EOL/archived**; no drop-in better short-read germline SV caller.
- GRIDSS → drop candidate (marginal gain for a single genome).
- TelomereHunter `latest` → pin 1.1.0; haplogrep3 `latest` → pin v3.3.2.

## Validation before merging any of the above
Per `CLAUDE.md`: run the affected steps on a **known sample** and diff against the previous run — pathogenic hit counts (ClinVar/VEP), diplotypes + phenotypes (PharmCAT/CPIC), `variants_used/variants_total` + raw deltas (PGS). Treat a scoring-file or cache version change as a new baseline, not a directly comparable result.
