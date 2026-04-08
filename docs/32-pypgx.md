# Step 32: Comprehensive Pharmacogenomics (pypgx)

## What This Does

Comprehensive pharmacogenomic star allele calling for 23 clinically actionable genes, including structural variation (SV) detection from BAM read depth. Complements PharmCAT (step 7) by covering genes that VCF-only callers miss entirely, most notably CYP2D6 gene deletions and duplications.

## Why

PharmCAT (step 7) calls ~23 genes from VCF data alone. This works well for simple SNP-based star alleles but fails for genes with structural variation — CYP2D6, CYP2A6, GSTM1, and GSTT1 all have common whole-gene deletions and duplications that VCF callers cannot represent. CYP2D6 alone affects 25% of all prescribed drugs, and PharmCAT frequently returns "Not called" for it.

Cyrius (step 21) was designed specifically for CYP2D6 but can fail on some WGS samples due to CYP2D7 pseudogene homology. pypgx uses a different read-depth algorithm that handles this homology more robustly.

pypgx also calls genes absent from PharmCAT entirely: COMT, MTHFR, ABCB1, GSTM1, GSTT1, and IFNL3.

## Tool

- **pypgx** v0.26.0 (Sboner Lab, Weill Cornell Medicine)
- License: Apache-2.0 (GPL-3.0 compatible)
- Publication: [Lee et al., 2019](https://doi.org/10.1002/cpt.1552)

## Docker Image

```
quay.io/biocontainers/pypgx:0.26.0--pyh7e72e81_0
```

## Prerequisites

### pypgx-bundle (required, one-time download)

pypgx requires a companion data bundle containing 1000 Genomes phasing panels (for Beagle haplotype estimation) and CNV classifier models. The bundle is **not included** in the Docker image and must be downloaded separately (~370 MB):

```bash
cd ${GENOME_DIR}/reference
git clone --branch 0.26.0 --depth 1 https://github.com/sbslee/pypgx-bundle.git
```

The bundle version must match the pypgx version (`0.26.0`). The script validates the bundle exists at `${GENOME_DIR}/reference/pypgx-bundle/` and exits with download instructions if missing.

### Input files

- BAM from alignment (step 2): `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam` (+ `.bai` index)
- VCF from variant calling (step 3): `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz` (+ `.tbi` index)

## Command

```bash
./scripts/32-pypgx.sh your_name
```

## Gene List

The script calls 23 curated genes spanning CPIC Level A/B recommendations and key genes PharmCAT misses.

### BAM-based genes (structural variation detection)

These genes have common whole-gene deletions, duplications, or hybrid alleles that cannot be detected from VCF data alone. pypgx uses BAM read depth to identify copy number changes.

| Gene | Why BAM-based | Clinical Impact |
|---|---|---|
| CYP2D6 | Deletions (*5), duplications (*1x2, *2x2), CYP2D7 hybrids | 25% of all drugs; codeine, tamoxifen, SSRIs |
| CYP2A6 | Whole-gene deletion (*4), duplications | Nicotine metabolism, tegafur |
| GSTM1 | Homozygous deletion (null genotype) in ~50% of population | Detoxification, carcinogen metabolism |
| GSTT1 | Homozygous deletion (null genotype) in ~20% of population | Detoxification, drug conjugation |

### VCF-based genes (SNP/indel star alleles)

These are called using both BAM and VCF data. Star alleles are determined from variant calls.

| Gene | Key Drugs | CPIC Level | Notes |
|---|---|---|---|
| CYP1A2 | Caffeine, clozapine, theophylline | B | |
| CYP2B6 | Efavirenz, methadone | A | |
| CYP2C9 | Warfarin, phenytoin, NSAIDs | A | |
| CYP2C19 | Clopidogrel, SSRIs, PPIs | A | |
| CYP3A4 | Tacrolimus (with CYP3A5) | B | |
| CYP3A5 | Tacrolimus | A | |
| CYP4F2 | Warfarin (vitamin K cycle) | B | |
| DPYD | Fluorouracil, capecitabine | A | |
| TPMT | Azathioprine, mercaptopurine | A | |
| NUDT15 | Azathioprine, mercaptopurine | A | |
| UGT1A1 | Irinotecan, atazanavir, bilirubin | A | |
| SLCO1B1 | Simvastatin, statins | A | |
| VKORC1 | Warfarin | A | |
| NAT2 | Isoniazid, hydralazine | A | |
| COMT | Catecholamine metabolism | -- | Not in PharmCAT |
| MTHFR | Folate metabolism, methotrexate | -- | Not in PharmCAT |
| ABCB1 | Drug efflux (broad substrate range) | -- | Not in PharmCAT |
| G6PD | Rasburicase, primaquine, dapsone | A | |
| IFNL3 | Peginterferon (historical, DAAs replaced) | A | |

## CYP2D6 Structural Variation Detection

CYP2D6 is the most complex pharmacogene. It has a tandemly duplicated pseudogene (CYP2D7) with >90% sequence identity, and common structural variants in the general population:

- **Gene deletion (*5)**: Entire CYP2D6 removed. Homozygous = poor metabolizer.
- **Gene duplication (*1x2, *2x2, etc.)**: Extra functional copies. Can produce ultra-rapid metabolizer status.
- **CYP2D6/CYP2D7 hybrids (*36, *13, etc.)**: Recombination between gene and pseudogene.

pypgx detects these by analyzing read depth across the CYP2D6/CYP2D7 locus. A drop in coverage indicates deletion; elevated coverage indicates duplication. VCF data alone cannot represent these copy number changes, which is why PharmCAT returns "Not called" for CYP2D6 in most WGS samples.

The BAM-based calling uses both `--variants` and `--depth-of-coverage` per the upstream pypgx WGS workflow, combining SNV/haplotype information with read-depth SV detection for the most complete genotype call.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/pypgx/`.

| File | Contents |
|---|---|
| `<gene>/results.zip` | Per-gene pypgx archive with genotype data |
| `${SAMPLE}_pypgx_summary.tsv` | Consolidated: gene, diplotype, phenotype, SV flag, source |
| `${SAMPLE}_pharmcat_comparison.tsv` | Side-by-side comparison with PharmCAT (if step 7 was run) |

### Summary TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| Diplotype | Star allele call (e.g., *1/*4) |
| Phenotype | Metabolizer status (e.g., Intermediate Metabolizer) |
| SV_detected | Whether structural variation was detected |
| Source | BAM (SV genes) or VCF (variant-based genes) |

### PharmCAT comparison TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| PharmCAT_diplotype | Diplotype from step 7 |
| pypgx_diplotype | Diplotype from this step |
| Match | Yes, No, pypgx only, or PharmCAT only |

## Runtime

~20-40 minutes for 23 genes. All genes run sequentially inside a single Docker container to avoid repeated container startup overhead.

## Resource Requirements

```
--cpus 4 --memory 8g
```

The BAM-based SV detection (CYP2D6, CYP2A6, GSTM1, GSTT1) is the most memory-intensive phase due to read-depth calculation across the locus.

## Comparison with PharmCAT

| Aspect | PharmCAT (step 7) | pypgx (step 32) |
|---|---|---|
| Input | VCF only | BAM + VCF |
| CYP2D6 SVs | Cannot detect | Read-depth detection |
| GSTM1/GSTT1 | Not called | Deletion detection |
| COMT, MTHFR | Not covered | Covered |
| Drug recommendations | Yes (HTML/JSON report) | No (star alleles only) |
| CPIC integration | Built-in | Manual lookup via [CPIC guidelines](https://cpicpgx.org/guidelines/) |
| Validation | Widely used (research tool) | Less extensively validated |

The two tools are complementary. PharmCAT provides drug recommendations for the genes it covers; pypgx extends coverage to genes PharmCAT cannot call. When both tools call the same gene, concordance is expected for simple genotypes but discrepancies can occur for complex haplotypes. Neither tool is definitively "correct" in all cases — discrepancies should be investigated by examining the underlying variant calls. Neither tool constitutes a clinical test; results should be confirmed by a certified pharmacogenomics laboratory before making prescribing decisions.

## Note on Aldy

[Aldy](https://github.com/0xTCG/aldy) is widely considered the best CYP2D6 caller, with superior handling of complex structural rearrangements and hybrid alleles. However, Aldy is released under a custom academic-only license that is incompatible with GPL-3.0 redistribution. pypgx (Apache-2.0) is the GPL-compatible alternative used in this pipeline. If you are using this pipeline for personal/academic analysis and Aldy's license terms are acceptable, it can be run separately.

## Limitations

- pypgx gene coverage (88 total) is broader than the 23 curated here. The curated list focuses on CPIC Level A/B genes and key PharmCAT gaps. Edit the `BAM_GENES` and `VCF_GENES` variables in the script to add more.
- Star allele definitions evolve. pypgx 0.26.0 uses a specific PharmVar database snapshot that may not include the latest allele definitions.
- SV detection accuracy depends on sequencing depth. 30X WGS is adequate; lower depths produce less reliable copy number calls.
- pypgx does not produce drug recommendations directly. Consult [CPIC guidelines](https://cpicpgx.org/guidelines/) to translate diplotypes into clinical actions. Note: step 27 (CPIC lookup) currently parses PharmCAT output only and cannot read pypgx results.
- **BAM-based SV genes (CYP2D6, CYP2A6, GSTM1, GSTT1) may fail** with pypgx 0.26.0 due to a pandas 2.x compatibility bug in the genotyping module (`'Series' object has no attribute 'Haplotype1'`). In practice, GSTT1 often succeeds while CYP2D6, CYP2A6, and GSTM1 fail. All VCF-based genes are unaffected. For CYP2D6, Cyrius (step 21) or Aldy are alternatives. This is an upstream pypgx issue — monitor [pypgx releases](https://github.com/sbslee/pypgx/releases) for a fix.
- Individual gene failures do not stop the pipeline. However, if **all** genes fail, the script exits with status 1 before generating the summary TSV — this signals a systemic problem (e.g., wrong BAM path, corrupted index, missing pypgx-bundle). Rerun with verbose output to identify the root cause. For partial failures, check the summary TSV for "FAILED" entries.

## Maintenance

- pypgx is pinned to `0.26.0` in `versions.env`. Check [pypgx releases](https://github.com/sbslee/pypgx/releases) periodically for updates to star allele definitions or algorithm improvements.
- The pypgx-bundle must match the pypgx version. When updating pypgx, re-download the bundle: `cd ${GENOME_DIR}/reference && rm -rf pypgx-bundle && git clone --branch <new_version> --depth 1 https://github.com/sbslee/pypgx-bundle.git`
- If you update pypgx, rerun on a known sample and compare diplotype calls against the previous version before adopting the new results.
- The curated gene list should be reviewed against [CPIC guideline updates](https://cpicpgx.org/guidelines/) at least quarterly.

## Links

- [pypgx documentation](https://pypgx.readthedocs.io/)
- [pypgx GitHub](https://github.com/sbslee/pypgx)
- [PharmVar database](https://www.pharmvar.org/) (star allele definitions)
- [CPIC guidelines](https://cpicpgx.org/guidelines/)
- [Aldy](https://github.com/0xTCG/aldy) (academic-only alternative for CYP2D6)
