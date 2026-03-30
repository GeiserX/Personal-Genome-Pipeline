# Step 23: Clinical Variant Filter

## What This Does

Extracts the small subset of clinically interesting variants from your VEP-annotated VCF. Instead of manually searching through 4-5 million variants, this step produces a focused list of ~200-500 variants that are rare AND functionally impactful.

## Why

The biggest challenge after running VEP annotation is: "I have millions of variants, what do I look at?" This step solves that by applying conservative filters to surface the variants most likely to be medically relevant.

## Tool

bcftools (already installed — no new Docker image needed)

## Docker Image

`staphb/bcftools:1.21`

## Input

- VEP-annotated VCF from step 13: `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf` (or `.vcf.gz`)

## Command

```bash
./scripts/23-clinical-filter.sh your_name
```

## What Gets Filtered

The script produces two variant sets that are merged:

### HIGH Impact Variants
- Stop-gained (premature stop codon — breaks the protein)
- Frameshift insertions/deletions (shifts reading frame — breaks the protein)
- Splice donor/acceptor (disrupts splicing — breaks the protein)
- Start-lost (no translation initiation)

Expected count: 100-200 per genome.

### Rare MODERATE Impact Variants
- Missense variants (amino acid change) with gnomAD allele frequency < 1%
- In-frame insertions/deletions with gnomAD AF < 1%

Expected count: 200-400 per genome after frequency filtering.

## Output

| File | Contents | Size |
|---|---|---|
| `${SAMPLE}_clinical.vcf.gz` | Combined clinically interesting VCF | < 5 MB |
| `${SAMPLE}_clinical_summary.tsv` | Human-readable tab-delimited table | < 1 MB |
| `${SAMPLE}_high_impact.vcf.gz` | HIGH impact variants only | < 2 MB |
| `${SAMPLE}_rare_moderate.vcf.gz` | Rare MODERATE variants only | < 3 MB |

## Runtime

~5-10 minutes (I/O-bound, reading the large VEP VCF)

## How to Use the Output

### Quick look at the summary

```bash
# View the most important variants
column -t ${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv | head -20
```

### Cross-reference with ClinVar

```bash
# Find which clinical variants are also in ClinVar
docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools isec -n=2 -w1 \
    /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz \
    /genome/clinvar/clinvar.vcf.gz \
    -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_clinvar.vcf.gz
```

### Load in a genome browser

The `_clinical.vcf.gz` file is small enough to load in [IGV Web](https://igv.org/app/) or [gene.iobio](https://gene.iobio.io/) for visual inspection.

## Limitations

- This is a **computational filter**, not a clinical interpretation
- Some pathogenic variants are LOW impact (e.g., synonymous variants affecting splicing, regulatory variants) and will be missed by this filter
- gnomAD frequency filtering depends on VEP having annotated the gnomAD fields correctly
- Always cross-reference findings with ClinVar and consult a professional for clinical decisions

## Notes

- No additional Docker images required — uses the same bcftools image as other steps
- The VEP VCF is compressed and indexed automatically if needed
- PASS filter is applied to exclude low-quality variant calls
