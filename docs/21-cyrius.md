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

- Sorted BAM with index from alignment (step 2):
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
| sample1 | *1/*2 |

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
- **Cyrius has not been updated since May 2021** (v1.1.1). A 2025 study (BCyrius, PMID 39901590) found Cyrius fails to call or miscalls 50/360 simulated samples (13.9%) due to its outdated star allele database. Consider Aldy as an alternative (see below).

## Recommended Alternative: Aldy

[Aldy](https://github.com/0xTCG/aldy) v4.8.3 is a leading CYP2D6 caller for short-read WGS data. A systematic comparison (Twesigomwe et al. 2020, PMID 32789024) found Aldy was "the best performing algorithm in calling CYP2D6 structural variants." It identifies 92.2% of currently defined minor star alleles (vs 85.6% for Cyrius) and is actively maintained with the current PharmVar database.

Aldy also calls 37 additional pharmacogenes (CYP2C19, CYP2B6, UGT1A1, NAT2, DPYD, SLCO1B1, etc.), which can supplement PharmCAT results.

**To use Aldy instead of Cyrius:**

```bash
# Install in a Python container (one-time, or build a custom image)
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  python:3.11-slim \
  bash -c "
    pip install -q aldy==4.8.3 &&
    aldy genotype \
      -p illumina \
      -g CYP2D6 \
      -o /genome/${SAMPLE}/cyrius/${SAMPLE}_aldy_cyp2d6.aldy \
      /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
  "
```

> **License note:** Aldy uses an academic/non-commercial license (IURTC, Indiana University). It is free for personal and research use but is NOT compatible with GPL-3.0 redistribution. This is why it is documented here as an optional recommendation rather than replacing Cyrius in the pipeline script. A GPL-compatible alternative (pypgx) is available in step 32.

## Notes

- The script creates a manifest file listing the BAM path, then runs Cyrius in a single container invocation.
- Cross-reference the Cyrius diplotype with PharmCAT's CYP2D6 call and pypgx (step 32). Cyrius can fail on some WGS samples due to CYP2D7 pseudogene homology; pypgx handles this more robustly. If all three disagree, Aldy (see above) has the broadest star allele coverage among available callers.
- Step 27 (CPIC lookup) currently reads PharmCAT JSON only. To translate a Cyrius diplotype into drug recommendations, consult [CPIC guidelines](https://cpicpgx.org/guidelines/) manually.

## Links

- [Cyrius GitHub](https://github.com/Illumina/Cyrius)
- [Aldy GitHub](https://github.com/0xTCG/aldy) — recommended alternative
- [PharmGKB CYP2D6](https://www.pharmgkb.org/gene/PA128)
- [CPIC CYP2D6 guidelines](https://cpicpgx.org/genes-drugs/)
- [Chen et al. 2021 (Cyrius paper)](https://doi.org/10.1038/s41397-021-00244-y)
- [Twesigomwe et al. 2020 (CYP2D6 caller comparison)](https://doi.org/10.1038/s41525-020-0135-2)
