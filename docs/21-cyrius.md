# Step 21: CYP2D6 Star Allele Calling with Cyrius

> **EXPERIMENTAL:** Cyrius is installed via pip at runtime inside a generic Python container, which is fragile (network dependency, version drift). Results should be cross-referenced with PharmCAT's CYP2D6 call.

## What This Does

Calls CYP2D6 star alleles (diplotypes) from your WGS BAM using Illumina's Cyrius tool. CYP2D6 is the single most important pharmacogene — it metabolizes roughly 25% of clinically used drugs — but its highly homologous pseudogene (CYP2D7) and frequent structural rearrangements (deletions, duplications, gene-pseudogene hybrids) make it extremely difficult to genotype from short reads.

## Why

PharmCAT (step 7) handles most pharmacogenes well, but its internal CYP2D6 calling is limited for WGS data. Cyrius was purpose-built by Illumina to resolve CYP2D6 using read-depth patterns across the CYP2D6/CYP2D7 region. Running Cyrius separately gives you a CYP2D6 diplotype that you can cross-reference with PharmCAT's results.

## Tool

- **Cyrius** (Chen et al., Pharmacogenomics J 2021) -- Illumina's depth-based CYP2D6 caller

## Docker Image

```
python:3.11-slim
```

Cyrius is installed via `pip install cyrius` inside the container at runtime. No dedicated Cyrius Docker image is required.

## Input

- Sorted BAM with index from alignment (step 1):
  - `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam`
  - `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai`

## Command

```bash
./scripts/21-cyrius.sh your_name
```

## What the Script Does Internally

1. Validates that the sorted BAM and its index exist
2. Creates a manifest file listing the BAM path (Cyrius requires this)
3. Installs Cyrius in a Python 3.11 container and runs `star_caller` with `--genome 38` (GRCh38)
4. Parses the output TSV to display the called diplotype

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_cyp2d6.tsv` | Tab-delimited results with sample name, diplotype, and supporting evidence |

All output is written to `${GENOME_DIR}/${SAMPLE}/cyrius/`.

## Runtime

~5-15 minutes (includes pip install overhead on first run).

## Interpreting Results

The output TSV contains the CYP2D6 diplotype in star-allele notation, for example:

| Sample | Genotype |
|---|---|
| sergio | *1/*2 |

Common results and what they mean:

- `*1/*1` -- Normal metabolizer (two fully functional copies)
- `*1/*2` -- Normal metabolizer (*2 is also functional)
- `*1/*4` -- Intermediate metabolizer (*4 is non-functional)
- `*4/*4` -- Poor metabolizer (no functional copies)
- `*1/*1xN` -- Ultrarapid metabolizer (gene duplication, N extra copies)
- `*5/*5` -- Poor metabolizer (whole gene deletion)

Look up your specific diplotype at [PharmGKB CYP2D6](https://www.pharmgkb.org/gene/PA128) for the corresponding metabolizer phenotype and drug implications.

### Drugs affected by CYP2D6 status

Codeine, tramadol, oxycodone, tamoxifen, ondansetron, atomoxetine, most tricyclic antidepressants (amitriptyline, nortriptyline), and paroxetine -- among many others.

## Limitations

- Cyrius works best with 30X+ WGS data. Lower coverage may produce uncertain calls.
- Rare hybrid alleles (e.g., *36, *68) may not be resolved.
- The pip-install-at-runtime approach adds startup time and requires internet. If you run this frequently, consider building a custom Docker image with Cyrius pre-installed.
- Cyrius only calls CYP2D6. For other pharmacogenes, rely on PharmCAT (step 7).

## Notes

- The script creates a manifest file listing the BAM path, then runs Cyrius in a single container invocation.
- Cross-reference the Cyrius diplotype with PharmCAT's CYP2D6 call. If they disagree, Cyrius is generally more reliable for WGS data.
- Feed the Cyrius result into CPIC lookups (step 27) for actionable drug recommendations.

## Links

- [Cyrius GitHub](https://github.com/Illumina/Cyrius)
- [PharmGKB CYP2D6](https://www.pharmgkb.org/gene/PA128)
- [CPIC CYP2D6 guidelines](https://cpicpgx.org/genes-drugs/)
- [Chen et al. 2021 (Cyrius paper)](https://doi.org/10.1038/s41397-021-00244-y)
