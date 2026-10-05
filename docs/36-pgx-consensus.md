# Step 36: PGx Consensus (PharmCAT Outside Calls)

## What This Does

Decides what the BAM-based callers tell PharmCAT. PharmCAT (step 7) reads a VCF, and from a VCF it types neither the HLA genes nor CYP2D6, so its drug guidance for them (abacavir, carbamazepine, allopurinol, codeine, tramadol, tamoxifen and many more) could never fire. This step writes PharmCAT's [outside-call file](https://pharmcat.org/using/Outside-Call-Format/):

- **HLA-A and HLA-B** from T1K (step 8), cut to two fields (`*57:01`), when every allele has a T1K quality above 0. One allele reads as homozygous, the way T1K reports a homozygous gene.
- **CYP2D6** only when pypgx (step 32) and Cyrius (step 21, opt-in) give the same diplotype **and** the [CYP2D6 depth check](32-pypgx.md#cyp2d6-depth-check) passed. One caller alone, a disagreement, or multi-mapped reads leave CYP2D6 `indeterminate`, and nothing reaches PharmCAT.

It also writes a consensus table that says, for each of the three genes, what every caller said and what was passed on. Step 27 copies it into the CPIC recommendations.

## Why

A single depth-based CYP2D6 call can be wrong in ways that change a prescription: on a BAM aligned to a reference with ALT contigs, the depth at CYP2D6 drops and pypgx can report a whole-gene deletion that is not there ([lessons learned](lessons-learned.md)). Two independent callers agreeing, on depth that passed its check, is the bar for letting a call steer drug guidance. PharmCAT 3.4 makes no CYP2D6 call of its own from a VCF (its report gives CYP2D6 `callSource` `NONE`), so it cannot be the second caller.

## Tool

- `bin/pgx_outside_calls.py` (standard library only), the same code as the `PGX_CONSENSUS` Nextflow module

## Docker Image

- `PYTHON_IMAGE`

## Input

Each one is read when it exists, and the table says when it does not:

| File | From |
|---|---|
| `${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv` | step 8 |
| `${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv` | step 32 |
| `${SAMPLE}/cyrius/${SAMPLE}_cyp2d6.tsv` | step 21 (opt-in) |
| `${SAMPLE}/pypgx/${SAMPLE}_cyp2d6_depth_check.tsv` (else step 21's) | the CYP2D6 depth check of step 32 or 21 |

## Command

```bash
./scripts/08-hla-typing.sh your_name
./scripts/32-pypgx.sh your_name
./scripts/21-cyrius.sh your_name          # opt-in: setup.sh --cyrius first
./scripts/36-pgx-consensus.sh your_name
./scripts/07-pharmacogenomics.sh your_name   # reads the outside calls
./scripts/27-cpic-lookup.sh your_name        # lists them
```

The Nextflow pipeline (and so `run-all.sh`) runs `PGX_CONSENSUS` for every sample with a BAM when `pharmcat` and any of `hla_typing`, `pypgx` or `cyrius` are in `--tools`, and PharmCAT waits for it.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/pgx_consensus/`.

| File | Contents |
|---|---|
| `${SAMPLE}_outside_calls.tsv` | PharmCAT's outside-call file: `gene<TAB>diplotype`, one line per gene passed on; empty when none was. Step 7 passes it with `-po` when it is not empty. |
| `${SAMPLE}_pgx_consensus.tsv` | `Gene`, `Result` (the diplotype, `indeterminate` or `not typed`), `Outside_call` (`yes` or `no`), `Reason`, `Evidence` (what each caller said) |

For CYP2D6 the `Reason` is one of: the depth check's message (multi-mapped reads), `no caller made a call`, `one caller only: a second caller must agree (Cyrius, --tools cyrius)`, `pypgx and Cyrius disagree`, `the CYP2D6 depth check did not run`, or `pypgx and Cyrius agree; the depth check passed`.

## Runtime

Seconds. PharmCAT (step 7) runs again afterwards, a few minutes.

## Notes

- Without Cyrius, CYP2D6 is always `indeterminate`: that is the point, not a fault. Install Cyrius (`setup.sh --cyrius`, non-commercial licence) and add it to the tools to get an agreed call.
- A Cyrius call counts only with Filter `PASS`; `None`, an ambiguous call (two diplotypes) or `CYP2D6_depth_unreliable` is no call.
- The two diplotypes must be the same up to the order of the alleles. The one passed to PharmCAT is pypgx's.
- PharmCAT turns the HLA alleles into its own phenotypes (`*57:01 positive`, `*58:01 negative` and so on) and gives the drug guidance for them; step 27 marks those genes `[outside call]`.
- Step 7 run before this step gives a report without the outside calls. Run step 7 again after it, as the command order above does.
- Steps 8, 21 and 32 delete this step's two files when they start: a call made from an earlier result never reaches PharmCAT. After running any of them again, run this step and step 7 again.
