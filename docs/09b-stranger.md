# Step 9b: STR Clinical Annotation (Stranger)

## What This Does

Annotates the ExpansionHunter VCF (step 9) with clinical pathogenicity status for each repeat locus. Raw repeat counts become labelled calls: **normal**, **pre_mutation**, or **full_mutation**, with disease name, OMIM number, inheritance mode, and the specific repeat-size thresholds used for classification.

This turns step 9 output from a table of numbers into actionable clinical calls.

## Why

ExpansionHunter reports the number of repeats at each locus but applies no pathogenicity judgement. Stranger adds that judgement using the ClinGen/OMIM short tandem repeat database, so you can immediately see whether a repeat count falls in the normal, pre-mutation, or full-mutation range without manually cross-referencing disease thresholds.

## Tool

- **Stranger** v0.10.2 (Clinical Genomics Stockholm) — annotates STR VCFs with pathogenicity labels from a curated repeat catalog covering HTT, FMR1, C9orf72, DMPK, ATXN*, RFC1, and ~40 other loci

## Docker Image

```
quay.io/biocontainers/stranger:0.10.2--pyhdfd78af_0
```

- Binary: `stranger` (on PATH)
- Bundled repeat catalog: clinical ClinGen/OMIM database (installed inside the container)

## Command

```bash
./scripts/09b-stranger.sh your_name
```

Step 09 must run first:

```bash
./scripts/09-expansion-hunter.sh your_name male   # or female
./scripts/09b-stranger.sh your_name
```

A custom repeat catalog (TSV) can be supplied via the `STRANGER_REPEATS` environment variable; the bundled catalog is used when it is not set.

## Output

| File | Description |
|---|---|
| `expansion_hunter/<sample>_eh_stranger.vcf` | Annotated VCF with `STR_STATUS` and disease metadata in INFO fields |

Key INFO fields added by Stranger:

| Field | Values | Meaning |
|---|---|---|
| `STR_STATUS` | `normal`, `pre_mutation`, `full_mutation` | Pathogenicity call for this allele |
| `Disease` | e.g. `Huntingtons disease` | Disease associated with this locus |
| `OMIM` | e.g. `143100` | OMIM disease identifier |
| `Inheritance` | `AD`, `AR`, `XD`, `XR` | Inheritance mode |
| `NormalMax` | integer | Upper bound of normal repeat range |
| `PathologicMin` | integer | Lower bound of clearly pathogenic range |

## Runtime

Under 1 minute. Stranger is a pure-Python VCF annotator — it reads the VCF once and writes to stdout.

## Notes

- Requires step 09 (ExpansionHunter) to have run first. The script exits cleanly with an informational message if the EH VCF is absent.
- The bundled catalog covers ~40 STR loci with established clinical thresholds. Custom catalogs can be used via `STRANGER_REPEATS=/path/to/catalog.tsv`.
- `STR_STATUS` is per-allele: a heterozygous locus may have one normal and one pre-mutation allele.
- Pre-mutation alleles at FMR1 (55-200 CGG) and ATXN1 carry carrier risk even without current disease. See `docs/interpreting-results.md` for guidance.
- Short-read WGS has limited ability to size very large expansions (>150 repeats) accurately; `full_mutation` calls at loci like FMR1 and C9orf72 should be confirmed with orthogonal methods.
- **RFC1 (CANVAS) is a special case — treat any flag as uninterpretable, not a diagnosis.** CANVAS requires the **AAGGG** motif specifically, **biallelic**, at **~400–2000+** repeats. Short-read ExpansionHunter reports only the degenerate `AARRG` motif and **cannot distinguish pathogenic AAGGG from the common benign AAAAG**, and the catalog's `STR_PATHOLOGIC_MIN` for RFC1 is far below the clinical threshold — so a modest expansion (e.g. 51/73) is over-called as `full_mutation`. Confirm only with motif-aware / flanking-PCR / repeat-primed-PCR testing, and only if cerebellar-ataxia/neuropathy/vestibular symptoms are present.
