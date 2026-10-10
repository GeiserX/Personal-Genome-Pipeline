# Step 17: Cancer Predisposition Screening with CPSR

## What This Does
Screens germline variants against curated cancer predisposition gene panels to identify clinically actionable cancer risk variants. CPSR uses its own panels sourced from Genomics England PanelApp and other curated databases — these are cancer-focused and distinct from the ACMG SF list of secondary findings (which also includes cardiac and metabolic genes not covered by CPSR).

## Why
ClinVar screening (step 6) finds known pathogenic variants, but CPSR applies ACMG/AMP classification criteria to novel or rare variants in cancer predisposition genes — catching variants ClinVar hasn't yet classified.

## Tool
- **CPSR** (Cancer Predisposition Sequencing Reporter), bundled inside the PCGR image

## Docker Image
- `PCGR_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

The image is PCGR 2.3.2 from `ghcr.io/sigven/pcgr` (Docker Hub's `sigven/pcgr` stops at 2.2.5). `cpsr` is on its `PATH`. It needs a separate ref data bundle (`PCGR_DATA_BUNDLE`, ~7 GB) and the VEP cache of `PCGR_VEP_CACHE_RELEASE`; the three move together.

## Prerequisites

### 1. VEP Cache
PCGR 2.3.2 bundles VEP 115, which requires the **release-115** cache (~24 GB). This is different from the release-116 cache used by step 13. Both coexist in the same `vep_cache/` directory under different subdirectories (`116_GRCh38/` and `115_GRCh38/`).
```bash
mkdir -p ${GENOME_DIR}/vep_cache
wget -c -P ${GENOME_DIR}/vep_cache https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf ${GENOME_DIR}/vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz -C ${GENOME_DIR}/vep_cache
```

### 2. PCGR Ref Data Bundle
PCGR 2.x uses a ref data bundle (~7 GB for `20260620`) separate from VEP:
```bash
mkdir -p ${GENOME_DIR}/pcgr_data
cd ${GENOME_DIR}/pcgr_data
wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20260620.grch38.tgz
tar xzf pcgr_ref_data.20260620.grch38.tgz
mkdir -p 20260620 && mv data/ 20260620/
```
This creates a `20260620/data/` directory with ClinVar, CancerMine, UniProt, and other databases. VEP cache is now mounted separately.

### Upgrading from PCGR 2.2.5
PCGR 2.2.5 read bundle `20250314` and the release-113 cache. 2.3.2 needs both new ones: about 31 GB to download (the 7 GB bundle and the 24 GB cache). Once step 17 runs on 2.3.2, `pcgr_data/20250314/` and `vep_cache/homo_sapiens/113_GRCh38/` are no longer read and can be deleted.

## Command
```bash
source versions.env   # from the repository root
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}/vep_cache:/mnt/.vep \
  -v ${GENOME_DIR}/pcgr_data/${PCGR_DATA_BUNDLE}:/mnt/bundle \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/mnt/inputs \
  -v ${GENOME_DIR}/${SAMPLE}/cpsr:/mnt/outputs \
  "${PCGR_IMAGE}" \
  cpsr \
    --input_vcf /mnt/inputs/${SAMPLE}.vcf.gz \
    --vep_dir /mnt/.vep \
    --refdata_dir /mnt/bundle \
    --output_dir /mnt/outputs \
    --genome_assembly grch38 \
    --sample_id ${SAMPLE} \
    --panel_id 0 \
    --secondary_findings \
    --force_overwrite
```

Run by hand like this, `--sample_id` must be 3 to 40 characters long. The script and the Nextflow module handle that for you (see Notes).

## Panel Options
`--panel_id` takes one or more of these ids, comma-separated. The list is the one PCGR 2.3.2 ships (`pcgr/pcgr_vars.py`); GEP means Genomics England PanelApp. 2.3 dropped 2.2.5's panel 38 (Rhabdoid tumour predisposition), so the ids from 38 up moved down by one.

| Panel ID | Description |
|---|---|
| 0 | CPSR exploratory cancer predisposition panel (PanelApp genes, TCGA's germline study, Cancer Gene Census, other sources), the widest panel and the one this step uses |
| 1 | Adult solid tumours cancer susceptibility (GEP) |
| 2 | Adult solid tumours for rare disease (GEP) |
| 3 | Bladder cancer pertinent cancer susceptibility (GEP) |
| 4 | Brain cancer pertinent cancer susceptibility (GEP) |
| 5 | Breast cancer pertinent cancer susceptibility (GEP) |
| 6 | Childhood solid tumours cancer susceptibility (GEP) |
| 7 | Colorectal cancer pertinent cancer susceptibility (GEP) |
| 8 | Endometrial cancer pertinent cancer susceptibility (GEP) |
| 9 | Familial Tumours Syndromes of the central & peripheral Nervous system (GEP) |
| 10 | Familial breast cancer (GEP) |
| 11 | Familial melanoma (GEP) |
| 12 | Familial prostate cancer (GEP) |
| 13 | Familial rhabdomyosarcoma (GEP) |
| 14 | GI tract tumours (GEP) |
| 15 | Genodermatoses with malignancies (GEP) |
| 16 | Haematological malignancies cancer susceptibility (GEP) |
| 17 | Haematological malignancies for rare disease (GEP) |
| 18 | Head and neck cancer pertinent cancer susceptibility (GEP) |
| 19 | Inherited MMR deficiency (Lynch syndrome) (GEP) |
| 20 | Inherited non-medullary thyroid cancer (GEP) |
| 21 | Inherited ovarian cancer (without breast cancer) (GEP) |
| 22 | Inherited pancreatic cancer (GEP) |
| 23 | Inherited polyposis and early onset colorectal cancer (GEP) |
| 24 | Inherited predisposition to acute myeloid leukaemia (AML) (GEP) |
| 25 | Inherited susceptibility to acute lymphoblastoid leukaemia (ALL) (GEP) |
| 26 | Inherited predisposition to GIST (GEP) |
| 27 | Inherited renal cancer (GEP) |
| 28 | Inherited phaeochromocytoma and paraganglioma (GEP) |
| 29 | Melanoma pertinent cancer susceptibility (GEP) |
| 30 | Multiple endocrine tumours (GEP) |
| 31 | Multiple monogenic benign skin tumours (GEP) |
| 32 | Neuroendocrine cancer pertinent cancer susceptibility (GEP) |
| 33 | Neurofibromatosis Type 1 (GEP) |
| 34 | Ovarian cancer pertinent cancer susceptibility (GEP) |
| 35 | Parathyroid Cancer (GEP) |
| 36 | Prostate cancer pertinent cancer susceptibility (GEP) |
| 37 | Renal cancer pertinent cancer susceptibility (GEP) |
| 38 | Sarcoma cancer susceptibility (GEP) |
| 39 | Sarcoma susceptibility (GEP) |
| 40 | Thyroid cancer pertinent cancer susceptibility (GEP) |
| 41 | Tumour predisposition - childhood onset (GEP) |
| 42 | Upper gastrointestinal cancer pertinent cancer susceptibility (GEP) |
| 43 | DNA repair genes pertinent cancer susceptibility (GEP) |

## Output
- `${SAMPLE}.cpsr.grch38.html` — Interactive HTML report with classified variants
- `${SAMPLE}.cpsr.grch38.classification.tsv.gz` — Tab-separated variant classifications (gzipped; read it with `zcat`)
- Every variant in the panel genes gets one of the five ACMG/AMP classes: Pathogenic, Likely pathogenic, VUS, Likely benign, Benign
- In the TSV, `CLASSIFICATION` is the final class and `ASSERTION_AUTHORITY` says whose it is: ClinVar's, unless ClinVar has no record or a conflicted one, then CPSR's own (`CPSR_CLASSIFICATION`, which every variant has). PCGR 2.2.5 named the final column `FINAL_CLASSIFICATION`; the reports (`bin/collect_summary.py`) read either. The file names are the same in both versions.

## Runtime
~30-60 minutes per genome (depends on variant count).

## Notes
- **Sample id length:** CPSR accepts sample ids of 3 to 40 characters (PCGR 2.3.2, `SAMPLE_ID_MIN_LENGTH` and `SAMPLE_ID_MAX_LENGTH` in `pcgr/pcgr_vars.py`); any other length stops it before it reads the VCF. When your id is outside that range, the pipeline gives CPSR a padded or shortened id (`S1` becomes `S1_cpsr`, a longer id is cut to its first 40 characters, both by `bin/cpsr_sample_id`) and renames its files back to `${SAMPLE}.cpsr.*`, so CPSR's own report title shows that id.
- The ref data bundle (~7 GB) and VEP cache only need to be downloaded once — shared across all samples.
- **PCGR 2.x breaking changes:** The CLI changed completely from 1.x. The old `--pcgr_dir` flag (which internally appended `/data`) is replaced by `--refdata_dir` and `--vep_dir` as separate mount points. The single monolithic data bundle is split into a smaller ref data bundle + the standard Ensembl VEP cache. Docker volume mounts changed from a single `:/genome` to four separate mounts for VEP, bundle, inputs, and outputs.
- **Data bundle freshness:** The `20260620` bundle dates from June 2026 (ClinVar 2026-06, GENCODE 49). Check the [PCGR releases page](https://github.com/sigven/pcgr/releases) periodically for updated bundles — newer bundles include more recent ClinVar classifications and gene-disease annotations.
- **VEP cache version:** PCGR 2.3.2 requires the VEP release-115 cache, while step 13 uses release-116. Both coexist in `vep_cache/homo_sapiens/` (subdirectories `116_GRCh38/` and `115_GRCh38/`). You need both if running both steps.
- Use `--panel_id 0` for the comprehensive cancer superpanel (500+ genes) — cancer-focused.
- CPSR 2.3 classifies every variant in the target genes by itself; 2.2's `--classify_all` is gone and 2.3 refuses it. `--clinvar_trust_level` (0 to 4, default 0) decides when CPSR's class replaces ClinVar's; the step keeps the default.
- `--secondary_findings` reports pathogenic / likely-pathogenic variants in the ACMG SF incidental-findings gene list (cardiac and metabolic genes beyond CPSR's cancer panels; the list itself comes with the PCGR data bundle) in a dedicated section of the report. `scripts/17-cpsr.sh` passes it unless you set `CPSR_SECONDARY_FINDINGS=false`, for when you do not want findings outside cancer predisposition: `CPSR_SECONDARY_FINDINGS=false ./scripts/17-cpsr.sh your_sample`. The Nextflow module always passes it.
- CPSR is complementary to ClinVar screening — ClinVar finds known variants, CPSR classifies novel ones.
- The same ref data bundle is used by PCGR for somatic analysis (not relevant for germline WGS).
