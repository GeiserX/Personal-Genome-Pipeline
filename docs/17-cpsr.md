# Step 17: Cancer Predisposition Screening with CPSR

## What This Does
Screens germline variants against ACMG SF v3.2 (Secondary Findings) and curated cancer predisposition gene panels to identify clinically actionable cancer risk variants.

## Why
ClinVar screening (step 6) finds known pathogenic variants, but CPSR applies ACMG/AMP classification criteria to novel or rare variants in cancer predisposition genes — catching variants ClinVar hasn't yet classified.

## Tool
- **CPSR** (Cancer Predisposition Sequencing Reporter), bundled inside the PCGR image

## Docker Image
```
sigven/pcgr:1.4.1
```
CPSR binary is at `/usr/local/bin/cpsr` inside this image. Requires a separate data bundle (~21GB).

## Prerequisites
Download and extract the PCGR/CPSR data bundle:
```bash
mkdir -p ${GENOME_DIR}/pcgr_data
cd ${GENOME_DIR}/pcgr_data
wget -c http://insilico.hpc.uio.no/pcgr/pcgr.databundle.grch38.20220203.tgz
tar xzf pcgr.databundle.grch38.20220203.tgz
```
This creates a `data/` directory with VEP cache, ClinVar, CancerMine, UniProt, and other databases.

## Command
```bash
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  sigven/pcgr:1.4.1 \
  cpsr \
    --input_vcf /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    --pcgr_dir /genome/pcgr_data \
    --output_dir /genome/${SAMPLE}/cpsr \
    --genome_assembly grch38 \
    --sample_id ${SAMPLE} \
    --panel_id 0 \
    --classify_all \
    --force_overwrite
```

## Panel Options
| Panel ID | Description |
|---|---|
| 0 | Full ACMG SF v3.2 (81 genes) — recommended |
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
- The 21GB data bundle only needs to be downloaded once — shared across all samples.
- Use `--panel_id 0` for comprehensive screening (all ACMG SF genes).
- `--classify_all` ensures all variants in target genes get ACMG classification, not just known pathogenic.
- CPSR is complementary to ClinVar screening — ClinVar finds known variants, CPSR classifies novel ones.
- The same data bundle is used by PCGR for somatic analysis (not relevant for germline WGS).
