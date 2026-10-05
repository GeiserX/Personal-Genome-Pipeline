# Step 21: CYP2D6 Star Allele Calling with Cyrius

> **OPT-IN, NON-COMMERCIAL LICENCE:** a default run leaves this step out. Cyrius 1.1.1, as published on PyPI, is under the [PolyForm Strict License 1.0.0](https://polyformproject.org/licenses/strict/1.0.0): use for a non-commercial purpose only, and no distribution of it or of changed copies. (Its source files still carry GPL-3.0 headers from before Illumina changed the licence in December 2022; the licence of the package is PolyForm Strict.) It has had no release since May 2021. You install it yourself, once, with `./scripts/setup.sh --cyrius ${GENOME_DIR}`; the pipeline does not ship it.

## What This Does

Calls CYP2D6 star alleles (diplotypes) from your WGS BAM using Illumina's Cyrius tool. CYP2D6 is the single most important pharmacogene — it metabolizes roughly 25% of clinically used drugs — but its highly homologous pseudogene (CYP2D7) and frequent structural rearrangements (deletions, duplications, gene-pseudogene hybrids) make it extremely difficult to genotype from short reads.

## Why

PharmCAT (step 7) calls no CYP2D6 from a VCF. Cyrius was purpose-built by Illumina to resolve CYP2D6 using read-depth patterns across the CYP2D6/CYP2D7 region. It is the second CYP2D6 caller beside pypgx (step 32): [step 36](36-pgx-consensus.md) passes a CYP2D6 call to PharmCAT only when the two agree. Without Cyrius, no CYP2D6 call reaches PharmCAT, because one depth-based caller alone is not enough to steer drug guidance.

## Tool

- **Cyrius** (Chen et al., Pharmacogenomics J 2021) -- Illumina's depth-based CYP2D6 caller

## Docker Image

- `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

Cyrius is not in any image. `./scripts/setup.sh --cyrius ${GENOME_DIR}` installs it once, in this image, into `${GENOME_DIR}/tools/cyrius-1.1.1/`:

```bash
pip install --require-hashes --no-deps --only-binary :all: --target <dir> -r scripts/cyrius-constraints.txt
```

`scripts/cyrius-constraints.txt` pins Cyrius and every package it needs (pysam, numpy, scipy, statsmodels and theirs) with the sha256 of each wheel, so pip installs exactly those files: nothing resolved, nothing built from source. This install is the only part that uses the network. The step then runs `python3 -m cyrius` from that directory with no network, like every other step. The install records the image and the lock file it came from; when either changes, the step asks you to run `setup.sh --cyrius` again. No Cyrius image is published: one opt-in step under a non-commercial licence does not justify a registry image to maintain.

## Input

- Sorted BAM with index from alignment (step 2):
  - `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam`
  - `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam.bai`
- The Cyrius install of `./scripts/setup.sh --cyrius`

## Command

```bash
./scripts/setup.sh --cyrius "$GENOME_DIR"   # once
./scripts/21-cyrius.sh your_name
```

With `run-all.sh`, name it: `TOOLS=...,cyrius`. With Nextflow: `--tools ...,cyrius --cyrius_install ${GENOME_DIR}/tools/cyrius-1.1.1`.

## What the Script Does Internally

1. Validates that the sorted BAM, its index and the Cyrius install exist
2. Checks the depth at CYP2D6 against its flanks ([depth check](32-pypgx.md#cyp2d6-depth-check)), with `${MOSDEPTH_IMAGE}`
3. Creates a manifest file listing the BAM path (Cyrius requires this) and runs `python3 -m cyrius --genome 38` (GRCh38) from the install, with no network
4. When the depth check found multi-mapped reads, sets the call's Filter to `CYP2D6_depth_unreliable`: [step 36](36-pgx-consensus.md) then does not pass it on
5. Parses the output TSV to display the called diplotype

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_cyp2d6.tsv` | Tab-delimited results with sample name, diplotype, and Filter (`PASS`, a reason Cyrius gives, or `CYP2D6_depth_unreliable`) |
| `${SAMPLE}_cyp2d6_depth_check.tsv` | The CYP2D6 depth check: status (`ok` or `unreliable`), its message and the four depths |

All output is written to `${GENOME_DIR}/${SAMPLE}/cyrius/`.

## Runtime

~5-15 minutes. The one-time install (`setup.sh --cyrius`) takes about a minute.

## Interpreting Results

The output TSV contains the CYP2D6 diplotype in star-allele notation, for example:

| Sample | Genotype |
|---|---|
| `<sample>` | `*x/*y` |

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
- Cyrius only calls CYP2D6. For other pharmacogenes, rely on PharmCAT (step 7).
- **Cyrius has not been updated since May 2021** (v1.1.1). A 2025 study (BCyrius, PMID 39901590) found Cyrius fails to call or miscalls 50/360 simulated samples (13.9%) due to its outdated star allele database. Consider Aldy as an alternative (see below).

## Recommended Alternative: Aldy

[Aldy](https://github.com/0xTCG/aldy) v4.8.3 is a leading CYP2D6 caller for short-read WGS data. A systematic comparison (Twesigomwe et al. 2020, PMID 32789024) found Aldy was "the best performing algorithm in calling CYP2D6 structural variants." It identifies 92.2% of currently defined minor star alleles (vs 85.6% for Cyrius) and is actively maintained with the current PharmVar database.

Aldy also calls 37 additional pharmacogenes (CYP2C19, CYP2B6, UGT1A1, NAT2, DPYD, SLCO1B1, etc.), which can supplement PharmCAT results.

**To use Aldy instead of Cyrius:**

```bash
source versions.env   # from the repository root
# Install in a Python container (one-time, or build a custom image)
docker run --rm --user root \
  -v ${GENOME_DIR}:/genome \
  "${PYTHON_IMAGE}" \
  bash -c "
    pip install -q aldy==4.8.3 &&
    aldy genotype \
      -p illumina \
      -g CYP2D6 \
      -o /genome/${SAMPLE}/cyrius/${SAMPLE}_aldy_cyp2d6.aldy \
      /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
  "
```

> **License note:** Aldy uses an academic/non-commercial license (IURTC, Indiana University). It is free for personal and research use but is NOT compatible with GPL-3.0 redistribution, the same kind of restriction as Cyrius's PolyForm Strict licence. That is why both stay outside a default run: Aldy is documented here only, and Cyrius is opt-in. A GPL-compatible caller (pypgx) is step 32.

## Notes

- The script creates a manifest file listing the BAM path, then runs Cyrius in a single container invocation. It exits with an error when Cyrius writes no result file.
- [Step 36](36-pgx-consensus.md) compares the Cyrius diplotype with pypgx's (step 32). Only when both give the same diplotype and the depth check passed does PharmCAT get it, and with it the CPIC drug recommendations of step 27. Otherwise CYP2D6 is `indeterminate`, and the CPIC report says what each caller said. Cyrius can fail on some WGS samples due to CYP2D7 pseudogene homology. If the callers disagree, Aldy (see above) has the broadest star allele coverage among available callers.

## Links

- [Cyrius GitHub](https://github.com/Illumina/Cyrius)
- [Aldy GitHub](https://github.com/0xTCG/aldy) — recommended alternative
- [PharmGKB CYP2D6](https://www.pharmgkb.org/gene/PA128)
- [CPIC CYP2D6 guidelines](https://cpicpgx.org/genes-drugs/)
- [Chen et al. 2021 (Cyrius paper)](https://doi.org/10.1038/s41397-021-00244-y)
- [Twesigomwe et al. 2020 (CYP2D6 caller comparison)](https://doi.org/10.1038/s41525-020-0135-2)
