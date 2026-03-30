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
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 1: Intersect sample VCF with ClinVar pathogenic database
docker run --rm \
  --cpus 2 --memory 4g \
  -v ${GENOME_DIR}:/genome \
  staphb/bcftools:1.21 \
  bcftools isec \
    -n =2 -w 1 \
    /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    /genome/clinvar/clinvar_pathogenic_chr.vcf.gz \
    -Oz -o /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf.gz

# Step 2: Index the result
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/bcftools:1.21 \
  bcftools index -t /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf.gz

# Step 3: Extract human-readable summary
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/bcftools:1.21 \
  bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
    /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf.gz \
    > /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_summary.tsv
```

## Output
- `${SAMPLE}_clinvar_hits.vcf.gz` — VCF of sample variants overlapping ClinVar pathogenic entries
- `${SAMPLE}_clinvar_summary.tsv` — tab-separated summary with columns: CHROM, POS, REF, ALT, CLNSIG, CLNDN

## Interpreting Results
| Scenario | Meaning | Action |
|---|---|---|
| Heterozygous + autosomal recessive | Healthy carrier | Note for family planning only |
| Homozygous + autosomal recessive | Affected | Investigate — confirm with phenotype |
| Any genotype + autosomal dominant | Potentially affected | Investigate — check penetrance and phenotype |
| Compound het (two variants, same gene) | Potentially affected (recessive) | Check if variants are on different alleles (phasing) |

## Important Notes
- Most hits will be **heterozygous carriers of recessive conditions** — this is normal and expected
- A typical 30X WGS will show 20-50 ClinVar pathogenic overlaps, the vast majority being benign carrier states
- Focus review on: homozygous pathogenic, any autosomal dominant pathogenic, and compound heterozygous variants in the same gene
- The ClinVar pathogenic database must be chr-prefixed to match the BAM/VCF coordinate system (done in step 00)
- ClinVar is updated monthly — re-download periodically to catch newly classified variants
