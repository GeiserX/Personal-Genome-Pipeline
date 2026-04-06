# Step 6: ClinVar Pathogenic Variant Screening

## What This Does
Intersects your sample VCF against the ClinVar database of known pathogenic variants, identifying any positions where your genome carries a clinically reported disease variant.

## Why
ClinVar is the gold standard for known disease-causing variants. This screen catches pathogenic SNPs and indels that have been reported in clinical settings — carrier status, dominant disease risk, and pharmacogenomic flags.

## Tool
- **bcftools isec** — VCF intersection to find overlapping variants

## Docker Image
```
staphb/bcftools:1.21
```

## Prerequisites
- Sample VCF from DeepVariant (step 3)
- `clinvar_pathogenic_chr.vcf.gz` from reference setup (step 00) — chr-prefixed, filtered to Pathogenic/Likely_pathogenic only

## Command
```bash
export GENOME_DIR=/path/to/data
./scripts/06-clinvar-screen.sh <sample_name>

# For long-read Clair3 output:
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh <sample_name>
```

### What the Script Does

1. Filters the sample VCF to PASS variants only (`bcftools view -f PASS`)
2. Intersects the PASS VCF against the ClinVar pathogenic subset using `bcftools isec -p`
3. Reports the count of shared variants (positions in both the sample and ClinVar pathogenic)

## Output

| File | Description |
|---|---|
| `clinvar/${SAMPLE}_pass.vcf.gz` | PASS-only subset of the sample VCF (intermediate) |
| `clinvar/isec/0000.vcf` | Variants unique to the sample |
| `clinvar/isec/0001.vcf` | Variants unique to ClinVar pathogenic |
| `clinvar/isec/0002.vcf` | **Shared variants — your pathogenic hits** |
| `clinvar/isec/0003.vcf` | Shared variants (ClinVar's perspective) |

## Interpreting Results

This step screens against **Pathogenic and Likely_pathogenic variants only** — benign and VUS entries are excluded at the database level (see step 00 reference setup). Every hit in the output is at a position ClinVar classifies as disease-associated.

| Scenario | Meaning | Action |
|---|---|---|
| Heterozygous + autosomal recessive | Healthy carrier | Note for family planning only |
| Homozygous + autosomal recessive | Affected | Investigate — confirm with phenotype |
| Any genotype + autosomal dominant | Potentially affected | Investigate — check penetrance and phenotype |
| Compound het (two variants, same gene) | Potentially affected (recessive) | Check if variants are on different alleles (phasing) |

## Important Notes
- Most hits will be **heterozygous carriers of recessive conditions** — this is normal and expected
- A typical 30X WGS shows 0-10 pathogenic/likely pathogenic overlaps. The majority are benign carrier states for recessive conditions
- Focus review on: homozygous pathogenic, any autosomal dominant pathogenic, and compound heterozygous variants in the same gene
- The ClinVar pathogenic database must be chr-prefixed to match the BAM/VCF coordinate system (done in step 00)
- ClinVar is updated monthly — re-download periodically to catch newly classified variants
