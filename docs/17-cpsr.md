# Step 17: Cancer Predisposition Screening with CPSR

## What This Does
Screens germline variants against curated cancer predisposition gene panels to identify clinically actionable cancer risk variants. CPSR uses its own panels sourced from Genomics England PanelApp and other curated databases — these are cancer-focused and distinct from the 81-gene ACMG SF v3.2 list (which also includes cardiac and metabolic genes not covered by CPSR).

## Why
ClinVar screening (step 6) finds known pathogenic variants, but CPSR applies ACMG/AMP classification criteria to novel or rare variants in cancer predisposition genes — catching variants ClinVar hasn't yet classified.

## Tool
- **CPSR** (Cancer Predisposition Sequencing Reporter), bundled inside the PCGR image

## Docker Image
```
sigven/pcgr:2.2.5
```
CPSR binary is at `/usr/local/bin/cpsr` inside this image. Requires a separate ref data bundle (~5 GB) and a VEP cache.

## Prerequisites

### 1. VEP Cache
If you already have the VEP cache from step 13, reuse it. Otherwise:
```bash
mkdir -p ${GENOME_DIR}/vep_cache
wget -c -P ${GENOME_DIR}/vep_cache https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz
tar xzf ${GENOME_DIR}/vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz -C ${GENOME_DIR}/vep_cache
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
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}/vep_cache:/mnt/.vep \
  -v ${GENOME_DIR}/pcgr_data/20250314:/mnt/bundle \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/mnt/inputs \
  -v ${GENOME_DIR}/${SAMPLE}/cpsr:/mnt/outputs \
  sigven/pcgr:2.2.5 \
  cpsr \
    --input_vcf /mnt/inputs/${SAMPLE}.vcf.gz \
    --vep_dir /mnt/.vep \
    --refdata_dir /mnt/bundle \
    --output_dir /mnt/outputs \
    --genome_assembly grch38 \
    --sample_id ${SAMPLE} \
    --panel_id 0 \
    --classify_all \
    --force_overwrite
```

## Panel Options
| Panel ID | Description |
|---|---|
| 0 | Comprehensive cancer superpanel (500+ genes) — recommended |
| 1 | Adult-onset hereditary cancer |
| 2 | Childhood-onset hereditary cancer |
| 3 | Lynch syndrome |
| 4 | BRCA1/BRCA2 |

## Output
- `${SAMPLE}.cpsr.grch38.html` — Interactive HTML report with classified variants
- `${SAMPLE}.cpsr.grch38.snvs_indels.tiers.tsv` — Tab-separated variant classifications
- Variants classified into 5 tiers (Pathogenic → Benign) using ACMG/AMP criteria

## Runtime
~30-60 minutes per genome (depends on variant count).

## Notes
- The ref data bundle (~5 GB) and VEP cache only need to be downloaded once — shared across all samples.
- **PCGR 2.x breaking changes:** The CLI changed completely from 1.x. The old `--pcgr_dir` flag (which internally appended `/data`) is replaced by `--refdata_dir` and `--vep_dir` as separate mount points. The single monolithic data bundle is split into a smaller ref data bundle + the standard Ensembl VEP cache. Docker volume mounts changed from a single `:/genome` to four separate mounts for VEP, bundle, inputs, and outputs.
- **Data bundle freshness:** The `20250314` bundle dates from March 2025. Check the [PCGR releases page](https://github.com/sigven/pcgr/releases) periodically for updated bundles — newer bundles include more recent ClinVar classifications and gene-disease annotations.
- **VEP cache reuse:** If you already downloaded the VEP cache for step 13, the same directory works here. No duplicate download needed.
- Use `--panel_id 0` for the comprehensive cancer superpanel (500+ genes). Note: this is broader than ACMG SF but cancer-focused — it does not replace a full ACMG incidental-findings screen.
- `--classify_all` ensures all variants in target genes get ACMG classification, not just known pathogenic.
- CPSR is complementary to ClinVar screening — ClinVar finds known variants, CPSR classifies novel ones.
- The same ref data bundle is used by PCGR for somatic analysis (not relevant for germline WGS).
