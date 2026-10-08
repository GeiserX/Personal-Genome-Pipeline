# Step 9: Short Tandem Repeat (STR) Expansion Screening

## What This Does
Screens for pathogenic repeat expansions — a class of mutations invisible to both DeepVariant and Manta. In these diseases, a short DNA sequence (3-6 bases) gets repeated too many times.

## Why
STR expansions cause ~40 known neurological/neuromuscular diseases including Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, myotonic dystrophy, and multiple spinocerebellar ataxias.

## Tool
- **ExpansionHunter** v5.0.0 (Illumina) — upgraded from v2.5.5. Adds multithreading, improved long-repeat estimation, and a bundled GRCh38 variant catalog (31 pathogenic loci)

## Docker Image
- `EXPANSIONHUNTER_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

- Binary: `ExpansionHunter` (on PATH)
- GRCh38 catalog: `/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json` (31 pathogenic loci)

## Key Disease Thresholds
| Disease | Gene | Repeat Unit | Normal | Pathogenic |
|---|---|---|---|---|
| Huntington's | HTT | CAG | <27 | >35 |
| Fragile X | FMR1 | CGG | <45 | ≥55 (premutation) / >200 (full) |
| Friedreich's Ataxia | FXN | GAA | <33 | >66 |
| ALS/FTD | C9ORF72 | GGCCCC | <24 | >30 |
| Myotonic Dystrophy 1 | DMPK | CTG | <35 | >50 |
| SCA1 | ATXN1 | CAG | <33 | >39 |
| SCA2 | ATXN2 | CAG | <22 | >33 |

## FMR1 Clinical Zones

FMR1 (Fragile X) has four distinct clinical zones — the intermediate zone (45-54 repeats) is often omitted but clinically relevant:

| Zone | Repeats | Clinical Significance |
|---|---|---|
| Normal | <45 | No risk |
| Intermediate (gray zone) | 45-54 | Not affected, but repeats may expand in offspring. Genetic counseling recommended for carriers. |
| Premutation | 55-200 | Risk of FXTAS (tremor/ataxia, males >50), FXPOI (premature ovarian insufficiency). Offspring at risk of full expansion. |
| Full mutation | >200 | Fragile X syndrome (intellectual disability, behavioral features). Penetrance varies by sex and methylation. |

## Command

```bash
./scripts/09-expansion-hunter.sh your_name male
# or: ./scripts/09-expansion-hunter.sh your_name female
```

The second argument (`male`/`female`) is **required** — it affects X-linked loci (FMR1, AR): males have one allele, females have two.

## Notes
- Uses ExpansionHunter **v5.0.0** (`EXPANSIONHUNTER_IMAGE`)
- v5 CLI: `--reads`, `--reference`, `--variant-catalog` (JSON file), `--output-prefix` (auto-generates .vcf, .json)
- The 31-locus GRCh38 variant catalog is bundled inside the container at `/usr/local/share/ExpansionHunter/variant_catalog/grch38/variant_catalog.json`
- **Why the bundled catalog, not a larger one.** Illumina has published no newer catalog for ExpansionHunter 5.0.0 (its last release, 2021). The larger candidate is Stranger's own GRCh38 catalog, 51 loci with their normal and pathologic thresholds, which would keep step 9b's thresholds in step by construction. It is not used yet: it lacks two of the bundled loci (NIPA1 and NOTCH2NL), and ExpansionHunter has not been run on it here. Switching means passing it with `EH_CATALOG` (the script) and `--expansion_catalog` (the pipeline) after a real run; until then Stranger has thresholds for 29 of the 31 bundled loci, none for NIPA1 and NOTCH2NL.
- Short-read WGS can reliably detect expansions up to ~150 repeats; very large expansions (>1000) are less accurate
