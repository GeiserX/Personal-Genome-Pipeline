# Step 27: CPIC Drug-Gene Recommendation Lookup

## What This Does

Parses your PharmCAT results (step 7) to extract gene diplotypes and metabolizer phenotypes, then maps them against CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines to produce a plain-text report of medications that may require dosing adjustments based on your pharmacogenomic profile.

## Why

PharmCAT produces detailed JSON and HTML reports, but digging through them to find actionable drug recommendations takes time. This step distills the clinically relevant parts into a simple, readable report: for each gene where you are NOT a normal metabolizer, it lists the affected medications and points you to the corresponding CPIC guideline.

## Tool

- **Python 3.11** for JSON parsing of PharmCAT output
- **CPIC gene-drug pairs** hard-coded in the script (based on published CPIC guidelines)

## Docker Image

```
python:3.11-slim
```

## Input

- PharmCAT JSON report from step 7. The script searches for it in:
  - `${GENOME_DIR}/${SAMPLE}/pharmcat/${SAMPLE}.report.json`
  - `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.report.json`
  - `${GENOME_DIR}/${SAMPLE}/pharmcat/report.json`
  - Any `.json` file in the `pharmcat/` or `vcf/` directories

## Command

```bash
./scripts/27-cpic-lookup.sh your_name
```

## What the Script Does Internally

1. Locates the PharmCAT JSON report by scanning expected paths
2. Runs a Python script (inside Docker) to parse the JSON and extract each gene's diplotype and phenotype
3. Writes a formatted gene results table to the output report
4. For each gene where the phenotype is NOT normal/typical/extensive:
   - Looks up the gene in a built-in CPIC drug table
   - Lists all medications affected by that gene's altered function
   - Adds a link to the CPIC guidelines page
5. Appends a disclaimer

## Genes and Drugs Covered

The script includes CPIC-level gene-drug pairs for 16 pharmacogenes:

| Gene | Example affected drugs |
|---|---|
| CYP2C19 | Clopidogrel, escitalopram, omeprazole, sertraline |
| CYP2C9 | Warfarin, phenytoin, celecoxib, ibuprofen |
| CYP2D6 | Codeine, tramadol, tamoxifen, paroxetine, ondansetron |
| CYP3A5 | Tacrolimus |
| CYP2B6 | Efavirenz |
| DPYD | Fluorouracil, capecitabine |
| TPMT / NUDT15 | Azathioprine, mercaptopurine, thioguanine |
| UGT1A1 | Atazanavir, irinotecan |
| SLCO1B1 | Simvastatin, atorvastatin, rosuvastatin |
| VKORC1 | Warfarin |
| HLA-A / HLA-B | Carbamazepine, abacavir, allopurinol, phenytoin |
| IFNL3 | Peginterferon alfa |
| RYR1 / CACNA1S | Volatile anesthetics, succinylcholine |
| G6PD | Rasburicase |
| MT-RNR1 | Aminoglycosides |

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_cpic_recommendations.txt` | Human-readable report with gene results and drug recommendations |
| `${SAMPLE}_phenotypes.tsv` | Intermediate gene/diplotype/phenotype table parsed from PharmCAT |

All output is written to `${GENOME_DIR}/${SAMPLE}/cpic/`.

## Runtime

~1-2 minutes (mostly Docker startup overhead).

## Interpreting Results

The report has two sections:

### Gene Results Table

Lists every pharmacogene with its called diplotype and phenotype. For example:

```
Gene         Diplotype                 Phenotype
CYP2C19      *1/*17                    Rapid Metabolizer
CYP2D6       *1/*2                     Normal Metabolizer
UGT1A1       *28/*28                   Poor Metabolizer
```

### Affected Medications

Only genes where you are NOT a normal metabolizer appear here. For each, the report lists:
- Your phenotype and diplotype
- All drugs with CPIC recommendations for that gene
- A link to the CPIC guidelines

**Normal/typical/extensive metabolizers are intentionally skipped** -- they require no dosing changes.

### What to do with the results

1. Check if you currently take (or might be prescribed) any of the listed medications
2. For any matches, read the full CPIC guideline at [cpicpgx.org/guidelines](https://cpicpgx.org/guidelines/)
3. Share the report with your prescribing physician or pharmacist

## Limitations

- The CPIC drug list is hard-coded in the script. New CPIC guidelines published after the script was written will not be included until the script is updated.
- The script does not query the live CPIC API -- it uses a static lookup table. This is intentional for reproducibility and offline use.
- PharmCAT JSON format varies between versions. If parsing fails, the script outputs a warning and produces a partial report.
- CYP2D6 results from PharmCAT may be less reliable than Cyrius (step 21). Cross-reference both before acting on CYP2D6 recommendations.
- This is NOT medical advice. Always consult a healthcare professional before making medication changes.

## Notes

- Run this step after PharmCAT (step 7). For CYP2D6, also run Cyrius (step 21) and manually compare.
- The output report is printed to stdout as well as written to file.
- You can add or modify gene-drug pairs by editing the `CPIC_DRUGS` associative array in the script.
- For maintenance, review the hard-coded CPIC table at least quarterly or whenever you bump PharmCAT, so the lookup stays aligned with current guideline pairs.
- For the most up-to-date CPIC recommendations, always check [cpicpgx.org](https://cpicpgx.org/) directly.

## Links

- [CPIC Guidelines](https://cpicpgx.org/guidelines/)
- [PharmCAT](https://pharmcat.org/)
- [PharmGKB](https://www.pharmgkb.org/)
- [CPIC Gene-Drug Pairs Table](https://cpicpgx.org/genes-drugs/)
