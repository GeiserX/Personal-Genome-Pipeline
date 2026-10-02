# Step 24: HTML Report

## What This Does

Generates a self-contained HTML dashboard summarizing all pipeline results. Open it in any browser — no internet connection needed.

## Why

The pipeline produces output across many directories in different formats (VCF, TSV, HTML, TXT). This step consolidates everything into a single visual report with color-coded status indicators, variant counts, and key findings.

## Tool

bash + bcftools (for extracting counts from VCF files)

## Docker Image

`staphb/bcftools:1.21` (already used by other steps)

## Input

All output directories from previous pipeline steps. The script automatically detects which steps have been run.

## Command

```bash
./scripts/24-html-report.sh your_name
```

## Output

A single file: `${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html`

Typically 10-30 KB. Contains:
- **Variant Calling** — total variants, PASS count, SNPs, indels
- **ClinVar Screening** — hit count with top 20 detailed in a table
- **Pharmacogenomics** — PharmCAT report status
- **Structural Variants** — Manta, Delly, CNVpytor counts
- **Cancer Predisposition** — CPSR report status
- **Repeat Expansions** — key loci repeat counts (HTT, FMR1, C9orf72, etc.)
- **Ancestry & Identity** — haplogroup, ROH, telomere content
- **Mitochondrial** — chrM variant and heteroplasmy counts
- **Clinical Filter** — interesting variant counts from step 23

## Runtime

1-3 minutes (mostly Docker startup time for bcftools queries)

## How to Open

```bash
# macOS
open ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html

# Linux
xdg-open ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html

# Windows (WSL)
start ${GENOME_DIR}/${SAMPLE}/${SAMPLE}_report.html
```

## Notes

- The report is completely self-contained — all CSS is inline, no external dependencies
- Works offline in any modern browser
- Responsive layout (works on mobile/tablet)
- Steps that were not run show "N/A" or "Not run" — this is expected
- The report contains health findings: the ClinVar table lists your pathogenic and likely pathogenic variants, and other sections summarise pharmacogenomics, CPSR and the other steps. It holds no raw reads or full VCF, but share it only as you would share a medical record
- Re-run this script anytime to update the report after running additional steps
