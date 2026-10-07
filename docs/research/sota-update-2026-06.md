# SOTA Update Roadmap — mid-2026

A point-in-time review of every tool/container/database against its latest upstream release, with the recommended action. Versions confirmed from each project's GitHub `releases/latest` or vendor page.

> Status: **a dated snapshot from June 2026, not kept up to date.** Most of the bumps below were applied in v0.6.0 to v0.8.2. The Status column says what had happened to each row by 2026-10-02: **applied** (the pin in `versions.env` today), **open** (not done) or **superseded** (a later decision replaced the row). The current pins are on [Image versions](../versions.md). Any bump still needs a re-run on a known sample before it merges (see [`lessons-learned.md`](../lessons-learned.md)).

## Priority actions
1. **Nextflow strict-syntax** — already addressed on `main` (#30/#31; `nextflow.config` no longer uses a top-level `def`). *Applied since:* `conf/base.config` uses `process.resourceLimits` instead of `check_max()`. Still open for NF 26.x: the `vcfanno` optional-input scope; no CI job runs 26.x yet. The pipeline ran on **NF 25.10.4** then; it is validated on 25.10.8 today (`NEXTFLOW_VERSION` in `versions.env`).
2. *(Applied in v0.7.0.)* **DeepVariant 1.6.0 → 1.10.0** — accuracy gains, pangenome-aware reassembly, native long-read phasing. Re-calls variants → full re-run + concordance check vs the previous VCF.
3. *(VEP 116 applied in v0.7.0; dbNSFP open.)* **VEP 112 → 116** + cache 112 → 116 (must match) + **dbNSFP 5.3.1** (one source for REVEL + AlphaMissense + CADD + MetaRNN).
4. *(Open.)* **samtools/bcftools 1.20/1.21 → 1.23.1** — note the **CRAM 3.0 → 3.1 default flip at 1.22** (write `--output-fmt cram,version=3.0` if older readers must consume the CRAM).
5. **Add**: Cyrius (CYP2D6, already wired), AlphaMissense (added through vcfanno, step 30, not the VEP plugin), pgsc_calc (PRS, open), ACMG SF v3.3 secondary-findings (84 genes; open: CPSR 2.2.5 carries v3.2 and the PCGR 2.3.x image is not on Docker Hub).
6. **Refresh DBs** (operational — user-supplied reference data, not pinned in-repo): ClinVar (monthly), gnomAD v4.1, CADD v1.7, AlphaMissense hg38, PGS Catalog, IPD-IMGT/HLA 3.60+.

## Tool versions

| Tool | Current | Latest | Action | Notes | Status (2026-10-02) |
|---|---|---|---|---|---|
| minimap2 | 2.28 | 2.31 | bump | drop-in | applied |
| samtools | 1.20 | 1.23.1 | bump | CRAM 3.0→3.1 default at 1.22; keep htslib/samtools/bcftools aligned | open |
| bcftools | 1.21 | 1.23.1 | bump | fixes silent output truncation | open |
| DeepVariant | 1.6.0 | 1.10.0 | bump | re-call; standard WGS path unchanged | applied |
| fastp | 1.3.1 | 1.3.6 | bump | BGZF multithread hang fixes | applied |
| Manta | 1.6.0 | 1.6.0 (EOL, archived Oct 2025) | keep | still de-facto short-read germline SV caller; no successor | applied (kept) |
| Delly | 1.7.3 | 2.1.0 | bump (major) | short-read PE/SR path backward-compatible; re-test SV step | applied |
| CNVnator | 0.4.1 | — | **replace → CNVpytor 1.3.1** | same lab, maintained Python rewrite, CRAM+BAF | applied (CNVpytor 1.3.2) |
| GRIDSS | 2.13.2 | 2.13.2 (2022) | keep/optional | unmaintained; droppable for single-genome runs | applied (kept, opt-in since) |
| AnnotSV | 3.4.4 | 3.5.10 | bump | refresh bundled annotations | applied |
| VEP | release_112 | release_116 | bump | cache must match; dbNSFP 5.3.1 | applied (dbNSFP open) |
| PCGR/CPSR | 2.2.5 | 2.3.0 | bump | re-download data bundle; CPSR ACMG-SF mode | open (no 2.3.x image on Docker Hub) |
| PharmCAT | 3.2.0 | 3.2.0 | **keep (latest)** | feed Cyrius CYP2D6 as outside-call | superseded: 3.4.0 is out (2026-07-14), so 3.2.0 is no longer the latest; bump open |
| TelomereHunter | `latest` | 1.1.0 / TH2 | **pin 1.1.0** | stop using `latest`; eval TelomereHunter2 | superseded: pinned by digest instead (no versioned tags) |
| haplogrep3 | `latest` | v3.3.2 | **pin v3.3.2** | | superseded: pinned by digest instead; that digest runs 3.2.1 |
| T1K | 1.0.9 | 1.0.9 | keep | refresh IPD-IMGT/HLA ref (3.60+) | superseded: v1.0.10 is out; bump open |
| ExpansionHunter | 5.0.0 | 5.0.0 | keep | add `stranger` for annotation | applied (stranger is step 9b) |
| goleft | 0.2.4 | 0.2.6 | bump | drop-in | applied |
| mosdepth | 0.3.13 | 0.3.14 | bump | drop-in | applied |
| GATK | 4.6.1.0 | 4.6.2.0 | bump | | applied |
| Picard | 3.4.0 | 3.4.0 | keep | | applied (kept) |
| PLINK2 | 2.00a5.10 | 2.00a7.1 | bump | pin by build | superseded: `pgscatalog/plink2` has no 2.00a7.1 image; held at 2.00a5.10, the version pgsc_calc pins |
| MultiQC | 1.33 | 1.35 | bump | min Python 3.9 | applied |
| vcfanno | 0.3.7 | 0.3.9 | bump | | applied |
| slivar | 0.3.3 | 0.3.4 | bump | | applied |
| pypgx | 0.26.0 | 0.27.0 | hold | 0.27.0 has no matching bundle and failed every gene (see lessons-learned); bump image and bundle together | applied (held at 0.26.0) |
| Clair3 | 2.0.0 | 2.0.2 | bump | long-read path only | applied |
| Sniffles | 2.4 | 2.8.0 | bump | long-read path only | applied |

> Biocontainer tags carry a build suffix (e.g. `…1.23.1--h96c455f_0`); the **version** is fixed as above — pick the newest `_N` at pin time. Keep `versions.env` and each module's `container` tag in sync (CI enforces this).

## Database refreshes
- **ClinVar** — latest weekly `vcf_GRCh38/clinvar.vcf.gz` (pin the dated file). Operational refresh (user-supplied reference data; not pinned in this repo).
- **VEP cache** → release 116 (`homo_sapiens_vep_116_GRCh38.tar.gz`). Use `wget -c` (resumable; the 26 GB download cannot be resumed by VEP's `INSTALL.pl`).
- **dbNSFP 5.3.1** — single VEP `--plugin dbNSFP` source for REVEL + AlphaMissense + CADD + MetaRNN (lets you retire standalone annotators).
- **gnomAD v4.1/v4.1.1** GRCh38 (constraint recalculated, AN bug fixed). No v5 yet.
- **CADD v1.7** GRCh38 (ESM-1v + regulatory CNN + Zoonomia).
- **AlphaMissense hg38** — `AlphaMissense_hg38.tsv.gz` → `tabix -s1 -b2 -e2 -S1`. The pipeline reads it through vcfanno at step 30, not the VEP plugin.
- **PGS Catalog** — via `pgsc_calc` (don't hand-roll scoring files).

## Steps to add (ranked)
1. **Cyrius** CYP2D6 star-allele caller (CNV/hybrid alleles PharmCAT misses) — already wired into the default tool set; feed its diplotype into PharmCAT as an outside-call. *(No single caller settles CYP2D6 copy number. On a reference with ALT contigs, depth at CYP2D6 drops and a depth-based caller (pypgx or Cyrius) can report a deletion that is not there: compare CYP2D6 depth with its flanks first, and report CYP2D6 only when two callers agree. See lessons-learned.)*
2. *(Applied through vcfanno at step 30, not the VEP plugin.)* **AlphaMissense** via VEP plugin — easy, high value.
3. **pgsc_calc** (Nextflow, NF-26 compatible) — SOTA polygenic scoring.
4. **ACMG SF v3.3** (2025, 84 genes) via CPSR secondary-findings mode (PCGR 2.3.0). *Done as well:* the reports list the ClinVar hits and rare HIGH-impact clinical-filter records in the 84 genes (step 23).
5. Consolidate missense annotation on **dbNSFP 5.3.1**.
6. **stranger** for repeat-expansion annotation (with ExpansionHunter).
7. **mtDNA heteroplasmy** (mutserve or GATK Mutect2 mito-mode) + VEP gnomADMT.

Long-read-only (skip for short-read Illumina WGS — no methylation signal in the data): modkit. *Corrected:* vg-giraffe and pangenome-aware DeepVariant are short-read methods, not long-read-only: DeepVariant 1.10's pangenome-aware case study runs on Illumina HG003 reads (SNP F1 0.9977 and INDEL F1 0.9972 on chr20), and DeepVariant publishes an image for it. No opt-in script is added until its memory need on a whole genome, with the HPRC graph loaded, has been measured on a real machine.

## Drop / replace
- CNVnator → **CNVpytor**.
- Manta → keep but **EOL/archived**; no drop-in better short-read germline SV caller.
- GRIDSS → drop candidate (marginal gain for a single genome).
- TelomereHunter `latest` → pin 1.1.0; haplogrep3 `latest` → pin v3.3.2. *Superseded:* both are pinned by digest instead.

## Validation before merging any of the above
Per `CLAUDE.md`: run the affected steps on a **known sample** and diff against the previous run — pathogenic hit counts (ClinVar/VEP), diplotypes + phenotypes (PharmCAT/CPIC), `Variants_Matched/Variants_Total` + raw deltas (PGS). Treat a scoring-file or cache version change as a new baseline, not a directly comparable result.
