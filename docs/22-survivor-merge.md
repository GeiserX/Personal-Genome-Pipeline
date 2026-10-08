# Step 22: Structural Variant Consensus Merge

## What This Does

Keeps the structural variants (SVs) that two or more independent callers agree on. SURVIVOR merge pairs two calls when both of their breakpoints lie within 1,000 bp of each other, their SV type and strands agree, and the event is at least 50 bp long. Each consensus record says how many callers support it and which.

## Why

Individual SV callers each have distinct biases and false-positive profiles:

- **Manta**: fast, sensitive for smaller SVs (paired-end and split reads)
- **Delly**: strongest for inversions and balanced translocations (paired-end, split reads and depth)
- **CNVpytor**: best for large CNVs (read depth only)

An SV seen by two callers that use different signals is more likely to be real.

The step used to bin calls by chromosome, `int(POS/1000)` and SV type. That split one deletion called at positions 999 and 1001 over two bins, and counted two deletions of very different size that start in one bin as agreement. SURVIVOR compares both breakpoints, so neither happens: the e2e case `tests/e2e/sv-mito-telomere-steps-1-sv-merge.sh` runs both situations on the synthetic three-caller set in `tests/fixtures/sv/`.

## Tool

- **SURVIVOR** 1.0.7 (Jeffares et al., Nat Commun 2017), `SURVIVOR merge`
- **bcftools** for the PASS filter before the merge and the sorted, indexed output

Not chosen: `truvari collapse`, which also compares sequence similarity. It would need a comparison run against SURVIVOR on real calls first.

## Docker Image

- `SURVIVOR_IMAGE` (the merge), `BCFTOOLS_IMAGE` (filter, sort, index)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

At least two of the following (the script uses every one it finds):

| Caller | Expected path |
|---|---|
| Manta (step 4) | `${GENOME_DIR}/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz` |
| Delly (step 19) | `${GENOME_DIR}/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz` |
| CNVpytor (step 18) | `${GENOME_DIR}/${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.vcf.gz`, or `_cnvs.txt`, which the script turns into that VCF first |
| TIDDIT (script 4a) | `${GENOME_DIR}/${SAMPLE}/sv_tiddit/${SAMPLE}_sv.vcf.gz` |
| Sniffles2 (script 4c, long reads) | `${GENOME_DIR}/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz` |

GRIDSS (step 4b) is left out, and the script says so when its VCF exists. GRIDSS reports every event as a pair of breakends (`SVTYPE=BND`), which never match the DEL, DUP and INV records of the other callers.

The Nextflow pipeline (`survivor_merge` in `--tools`, with at least two of `manta`, `delly`, `cnvpytor`) merges the callers it ran.

## Command

```bash
./scripts/22-survivor-merge.sh your_name
```

## What the Script Does Internally

1. Turns CNVpytor's table into a VCF when step 18 left no VCF
2. Writes each caller's PASS (or unfiltered) records as plain VCF under `sv_merged/inputs/`, with one sample column named after the caller. SURVIVOR names its output columns after the input samples, so three inputs that all name the sample would give one name three times
3. Runs `SURVIVOR merge inputs/sv_files.txt 1000 2 1 1 0 50`: maximum breakpoint distance 1,000 bp, support from 2 callers, same type, same strands, no size-scaled distance, minimum size 50 bp
4. Checks that SURVIVOR wrote a VCF (it exits 0 when it cannot open an input), then sorts, compresses and indexes it
5. Prints the count, and the count per type and caller combination

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_sv_consensus.vcf.gz` | Consensus SVs supported by 2+ callers, with `SUPP` and `SUPP_VEC` in INFO |
| `${SAMPLE}_sv_consensus.vcf.gz.tbi` | Tabix index |
| `inputs/sv_files.txt`, `inputs/<caller>.vcf` | The PASS records SURVIVOR read, in `SUPP_VEC` order |

All output is written to `${GENOME_DIR}/${SAMPLE}/sv_merged/`.

`SUPP_VEC` has one digit per caller, in the order of `sv_files.txt` (manta, delly, cnvpytor, tiddit, sniffles2, for the ones found): `110` is Manta and Delly. Each sample column holds that caller's genotype and its original record ID.

## Runtime

A few minutes (mostly reading the input VCFs).

## Interpreting Results

A typical 30X WGS genome produces:

- **Manta**: 7,000-9,000 SVs (see [step 4](04-structural-variants.md))
- **Delly**: 5,000-15,000 SVs
- **CNVpytor**: 3,000-4,000 CNVs, 1,500-2,000 of them with e-value < 0.01 (see [interpreting results](interpreting-results.md#cnvpytor-results-step-18))

Multi-caller SVs have lower false-positive rates than single-caller calls. The count after the merge has not been measured on a real genome with this version of the step.

SV types in the output: **DEL** (deletion), **DUP** (duplication), **INV** (inversion), **INS** (insertion), **TRA** (a translocation, SURVIVOR's name for a breakend pair between two chromosomes).

### Quick inspection

```bash
source versions.env   # from the repository root
# Consensus SVs by type and by which callers support them
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools query -f '%INFO/SVTYPE\t%INFO/SUPP_VEC\n' \
    /genome/${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz | sort | uniq -c | sort -rn
```

## Limitations

- Single-caller SVs are discarded even if they are real. If you suspect a specific SV, check the individual caller outputs directly.
- CNVpytor's calls are depth-only and their breakpoints are coarse (its bin size), so a CNVpytor call often lies more than 1 kb from the Manta or Delly breakpoint of the same event and does not pair with it.
- The CNVpytor VCF made from its table carries no genotype and no quality.
- Breakend (BND) records pair only with breakends of another caller.

## Notes

- Run this step only after completing at least two of: step 4 (Manta), step 19 (Delly), step 18 (CNVpytor).
- The consensus VCF can be annotated or loaded into IGV for visual inspection.

## Links

- [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
- [Manta](https://github.com/Illumina/manta)
- [Delly](https://github.com/dellytools/delly)
- [CNVpytor](https://github.com/abyzovlab/CNVpytor)
