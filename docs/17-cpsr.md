# Step 17: Cancer Predisposition Screening with CPSR

## What This Does
Screens germline variants against curated cancer predisposition gene panels to identify clinically actionable cancer risk variants. CPSR uses its own panels sourced from Genomics England PanelApp and other curated databases — these are cancer-focused and distinct from the 81-gene ACMG SF v3.2 list (which also includes cardiac and metabolic genes not covered by CPSR).

## Why
ClinVar screening (step 6) finds known pathogenic variants, but CPSR applies ACMG/AMP classification criteria to novel or rare variants in cancer predisposition genes — catching variants ClinVar hasn't yet classified.

## Tool
- **CPSR** (Cancer Predisposition Sequencing Reporter), bundled inside the PCGR image

## Docker Image
- `PCGR_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

CPSR binary is at `/usr/local/bin/cpsr` inside this image. Requires a separate ref data bundle (~5 GB) and a VEP cache.

## Prerequisites

### 1. VEP Cache
PCGR 2.2.5 bundles VEP 113, which requires the **release-113** cache. This is different from the release-116 cache used by step 13. Both coexist in the same `vep_cache/` directory under different subdirectories (`116_GRCh38/` and `113_GRCh38/`).
```bash
mkdir -p ${GENOME_DIR}/vep_cache
wget -c -P ${GENOME_DIR}/vep_cache https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf ${GENOME_DIR}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz -C ${GENOME_DIR}/vep_cache
```

### 2. PCGR Ref Data Bundle
PCGR 2.x uses a new, smaller ref data bundle (~5 GB) separate from VEP:
```bash
mkdir -p ${GENOME_DIR}/pcgr_data
cd ${GENOME_DIR}/pcgr_data
wget -c https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20250314.grch38.tgz
tar xzf pcgr_ref_data.20250314.grch38.tgz
mkdir -p 20250314 && mv data/ 20250314/
```
This creates a `20250314/data/` directory with ClinVar, CancerMine, UniProt, and other databases. VEP cache is now mounted separately.

## Command
```bash
source versions.env   # from the repository root
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}/vep_cache:/mnt/.vep \
  -v ${GENOME_DIR}/pcgr_data/20250314:/mnt/bundle \
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
    --classify_all \
    --secondary_findings \
    --force_overwrite
```

## Panel Options
`--panel_id` takes one or more of these ids, comma-separated. The list is the one CPSR 2.2.5 ships (`pcgr/pcgr_vars.py`); GEP means Genomics England PanelApp.

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
| 38 | Rhabdoid tumour predisposition (GEP) |
| 39 | Sarcoma cancer susceptibility (GEP) |
| 40 | Sarcoma susceptibility (GEP) |
| 41 | Thyroid cancer pertinent cancer susceptibility (GEP) |
| 42 | Tumour predisposition - childhood onset (GEP) |
| 43 | Upper gastrointestinal cancer pertinent cancer susceptibility (GEP) |
| 44 | DNA repair genes pertinent cancer susceptibility (GEP) |

## Output
- `${SAMPLE}.cpsr.grch38.html` — Interactive HTML report with classified variants
- `${SAMPLE}.cpsr.grch38.classification.tsv.gz` — Tab-separated variant classifications (gzipped; read it with `zcat`)
- Every variant in the panel genes gets one of the five ACMG/AMP classes: Pathogenic, Likely pathogenic, VUS, Likely benign, Benign

## Runtime
~30-60 minutes per genome (depends on variant count).

## Notes
- The ref data bundle (~5 GB) and VEP cache only need to be downloaded once — shared across all samples.
- **PCGR 2.x breaking changes:** The CLI changed completely from 1.x. The old `--pcgr_dir` flag (which internally appended `/data`) is replaced by `--refdata_dir` and `--vep_dir` as separate mount points. The single monolithic data bundle is split into a smaller ref data bundle + the standard Ensembl VEP cache. Docker volume mounts changed from a single `:/genome` to four separate mounts for VEP, bundle, inputs, and outputs.
- **Data bundle freshness:** The `20250314` bundle dates from March 2025. Check the [PCGR releases page](https://github.com/sigven/pcgr/releases) periodically for updated bundles — newer bundles include more recent ClinVar classifications and gene-disease annotations.
- **VEP cache version:** PCGR 2.2.5 requires VEP release-113 cache, while step 13 uses release-116. Both coexist in `vep_cache/homo_sapiens/` (subdirectories `116_GRCh38/` and `113_GRCh38/`). You need both if running both steps.
- Use `--panel_id 0` for the comprehensive cancer superpanel (500+ genes) — cancer-focused.
- `--classify_all` ensures all variants in target genes get ACMG classification, not just known pathogenic.
- `--secondary_findings` reports pathogenic / likely-pathogenic variants in the **ACMG SF v3.2** incidental-findings gene list (81 genes — including cardiac and metabolic genes beyond CPSR's cancer panels) in a dedicated section of the report. `scripts/17-cpsr.sh` passes it unless you set `CPSR_SECONDARY_FINDINGS=false`, for when you do not want findings outside cancer predisposition: `CPSR_SECONDARY_FINDINGS=false ./scripts/17-cpsr.sh your_sample`. The Nextflow module always passes it. The fuller v3.3 (84-gene) list follows once PCGR's 2.3.0 container image is published (the 2.3.0 source release is out, but its Docker image is not yet available, so the pin stays at PCGR 2.2.5).
- CPSR is complementary to ClinVar screening — ClinVar finds known variants, CPSR classifies novel ones.
- The same ref data bundle is used by PCGR for somatic analysis (not relevant for germline WGS).
