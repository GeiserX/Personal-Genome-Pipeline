# Step 32: Comprehensive Pharmacogenomics (pypgx)

## What This Does

Comprehensive pharmacogenomic star allele calling for 23 clinically actionable genes, including structural variation (SV) detection from BAM read depth. Complements PharmCAT (step 7) by covering genes that VCF-only callers miss entirely, most notably CYP2D6 gene deletions and duplications.

## Why

PharmCAT (step 7) calls ~23 genes from VCF data alone. This works well for simple SNP-based star alleles but fails for genes with structural variation — CYP2D6, CYP2A6, GSTM1, and GSTT1 all have common whole-gene deletions and duplications that VCF callers cannot represent. CYP2D6 alone affects 25% of all prescribed drugs, and PharmCAT frequently returns "Not called" for it.

Cyrius (step 21) was designed specifically for CYP2D6 but fails on most WGS samples due to CYP2D7 pseudogene homology. pypgx uses a different read-depth algorithm that handles this homology more robustly.

pypgx also calls genes absent from PharmCAT entirely: COMT, MTHFR, ABCB1, GSTM1, GSTT1, and IFNL3.

## Tool

- **pypgx** v0.26.0 (Sboner Lab, Weill Cornell Medicine)
- License: Apache-2.0 (GPL-3.0 compatible)
- Publication: [Lee et al., 2019](https://doi.org/10.1002/cpt.1552)

## Docker Image

```
quay.io/biocontainers/pypgx:0.26.0--pyh7e72e81_0
```

## Input

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

| Gene | Key Drugs | CPIC Level |
|---|---|---|
| CYP1A2 | Caffeine, clozapine, theophylline | B |
| CYP2B6 | Efavirenz, methadone | A |
| CYP2C9 | Warfarin, phenytoin, NSAIDs | A |
| CYP2C19 | Clopidogrel, SSRIs, PPIs | A |
| CYP3A4 | Tacrolimus (with CYP3A5) | B |
| CYP3A5 | Tacrolimus | A |
| CYP4F2 | Warfarin (vitamin K cycle) | B |
| DPYD | Fluorouracil, capecitabine | A |
| TPMT | Azathioprine, mercaptopurine | A |
| NUDT15 | Azathioprine, mercaptopurine | A |
| UGT1A1 | Irinotecan, atazanavir, bilirubin | A |
| SLCO1B1 | Simvastatin, statins | A |
| VKORC1 | Warfarin | A |
| NAT2 | Isoniazid, hydralazine | A |
| COMT | Catecholamine metabolism | Not in PharmCAT |
| MTHFR | Folate metabolism, methotrexate | Not in PharmCAT |
| ABCB1 | Drug efflux (broad substrate range) | Not in PharmCAT |
| G6PD | Rasburicase, primaquine, dapsone | A |
| IFNL3 | Peginterferon (historical, DAAs replaced) | A |

## CYP2D6 Structural Variation Detection

CYP2D6 is the most complex pharmacogene. It has a tandemly duplicated pseudogene (CYP2D7) with >90% sequence identity, and common structural variants in the general population:

- **Gene deletion (*5)**: Entire CYP2D6 removed. Homozygous = poor metabolizer.
- **Gene duplication (*1x2, *2x2, etc.)**: Extra functional copies. Can produce ultra-rapid metabolizer status.
- **CYP2D6/CYP2D7 hybrids (*36, *13, etc.)**: Recombination between gene and pseudogene.

pypgx detects these by analyzing read depth across the CYP2D6/CYP2D7 locus. A drop in coverage indicates deletion; elevated coverage indicates duplication. This is fundamentally impossible from VCF data alone, which is why PharmCAT returns "Not called" for CYP2D6 in most WGS samples.

The BAM-based calling for CYP2D6 omits the `--vcf` flag intentionally, so pypgx relies entirely on read-depth patterns for SV detection rather than mixing in potentially misleading VCF calls from the pseudogene-confounded region.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/pypgx/`.

| File | Contents |
|---|---|
| `<gene>/results.zip` | Per-gene pypgx archive with genotype data |
| `${SAMPLE}_pypgx_summary.tsv` | Consolidated: gene, diplotype, phenotype, activity score, SV flag, source |
| `${SAMPLE}_pharmcat_comparison.tsv` | Side-by-side comparison with PharmCAT (if step 7 was run) |

### Summary TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| Diplotype | Star allele call (e.g., *1/*4) |
| Phenotype | Metabolizer status (e.g., Intermediate Metabolizer) |
| ActivityScore | Numeric activity score where applicable |
| SV_detected | Whether structural variation was detected |
| Source | BAM (SV genes) or BAM+VCF |

### PharmCAT comparison TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| PharmCAT_diplotype | Diplotype from step 7 |
| pypgx_diplotype | Diplotype from this step |
| Match | Yes, No, pypgx only, or PharmCAT only |

## Runtime

~15-30 minutes for 23 genes. All genes run sequentially inside a single Docker container to avoid repeated container startup overhead.

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
| CPIC integration | Built-in | Requires step 27 |
| Clinical validation | Hospital-grade | Research-grade |

The two tools are complementary. PharmCAT provides clinical-grade drug recommendations for the genes it covers; pypgx extends coverage to genes PharmCAT cannot call. When both tools call the same gene, concordance is expected for simple genotypes but discrepancies can occur for complex haplotypes. Neither tool is definitively "correct" in all cases — discrepancies should be investigated by examining the underlying variant calls.

## Note on Aldy

[Aldy](https://github.com/inumanag/aldy) is widely considered the best CYP2D6 caller, with superior handling of complex structural rearrangements and hybrid alleles. However, Aldy is released under a custom academic-only license that is incompatible with GPL-3.0 redistribution. pypgx (Apache-2.0) is the GPL-compatible alternative used in this pipeline. If you are using this pipeline for personal/academic analysis and Aldy's license terms are acceptable, it can be run separately.

## Limitations

- pypgx gene coverage (88 total) is broader than the 23 curated here. The curated list focuses on CPIC Level A/B genes and key PharmCAT gaps. Edit the `BAM_GENES` and `VCF_GENES` variables in the script to add more.
- Star allele definitions evolve. pypgx 0.26.0 uses a specific PharmVar database snapshot that may not include the latest allele definitions.
- SV detection accuracy depends on sequencing depth. 30X WGS is adequate; lower depths produce less reliable copy number calls.
- pypgx does not produce drug recommendations directly. Use step 27 (CPIC lookup) or consult [CPIC guidelines](https://cpicpgx.org/guidelines/) to translate diplotypes into clinical actions.
- Individual gene failures do not stop the pipeline. Check the summary TSV for "FAILED" entries and investigate per-gene logs in the output directory.

## Maintenance

- pypgx is pinned to `0.26.0` in `versions.env`. Check [pypgx releases](https://github.com/sbslee/pypgx/releases) periodically for updates to star allele definitions or algorithm improvements.
- If you update pypgx, rerun on a known sample and compare diplotype calls against the previous version before adopting the new results.
- The curated gene list should be reviewed against [CPIC guideline updates](https://cpicpgx.org/guidelines/) at least quarterly.

## Links

- [pypgx documentation](https://pypgx.readthedocs.io/)
- [pypgx GitHub](https://github.com/sbslee/pypgx)
- [PharmVar database](https://www.pharmvar.org/) (star allele definitions)
- [CPIC guidelines](https://cpicpgx.org/guidelines/)
- [Aldy](https://github.com/inumanag/aldy) (academic-only alternative for CYP2D6)
