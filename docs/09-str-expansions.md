# Step 9: Short Tandem Repeat (STR) Expansion Screening

## What This Does
Screens for pathogenic repeat expansions — a class of mutations invisible to both DeepVariant and Manta. In these diseases, a short DNA sequence (3-6 bases) gets repeated too many times.

## Why
STR expansions cause ~40 known neurological/neuromuscular diseases including Huntington's, Fragile X, Friedreich's ataxia, ALS/FTD, myotonic dystrophy, and multiple spinocerebellar ataxias.

## Tool
- **ExpansionHunter** v2.5.5 (Illumina) — note: v5.x is available upstream with an expanded catalog and `--variant-catalog` flag, but the Docker image used here ships v2.5.5

## Docker Image
```
weisburd/expansionhunter:latest
```
- Binary: `/ExpansionHunter/bin/ExpansionHunter`
- GRCh38 catalogs: `/pathogenic_repeats/GRCh38/` (38 disease loci)

## Key Disease Thresholds
| Disease | Gene | Repeat Unit | Normal | Pathogenic |
|---|---|---|---|---|
| Huntington's | HTT | CAG | <27 | >35 |
| Fragile X | FMR1 | CGG | <45 | >54 (intermediate) / >200 (full) |
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

## Notes
- This is v2.5.5, NOT v5.x. The flag is `--repeat-specs` (directory), not `--variant-catalog`
- `--log` is REQUIRED (missing it causes silent exit)
- `--sex` affects X-linked loci (FMR1, AR): males have one allele, females have two
- Short-read WGS can reliably detect expansions up to ~150 repeats; very large expansions (>1000) are less accurate
