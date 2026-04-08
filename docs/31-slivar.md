# Step 31: Variant Prioritization (slivar)

Prioritizes clinically interesting variants using tiered filters and detects compound heterozygote candidates. Optionally annotates results with gnomAD gene constraint metrics.

## Why

A typical VEP-annotated VCF contains thousands of MODERATE/HIGH impact variants. Most are common population variants or benign polymorphisms. This step filters down to the variants most likely to be clinically relevant using a three-tier approach, then flags potential compound heterozygotes that could cause autosomal recessive disease.

slivar (by Brent Pedersen, author of vcfanno, mosdepth, duphold) is a streaming VCF filter that replaces the unmaintained GEMINI database approach. It uses JS expressions for flexible filtering without loading variants into a database.

## Prerequisites

- VEP-annotated VCF from step 13, or vcfanno-enriched VCF from step 30 (recommended)
- Optional: gnomAD v4.1 gene constraint TSV at `${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv`

## Docker Images

```text
quay.io/biocontainers/slivar:0.3.3--h5f107b1_0    # compound het detection
staphb/bcftools:1.21                                # variant filtering via split-vep
```

## Usage

```bash
export GENOME_DIR=/path/to/your/data
./scripts/31-slivar.sh your_sample
```

## Filter Tiers

### Tier 1: rare_high
- PASS variants with HIGH VEP impact (stop-gain, frameshift, splice donor/acceptor)
- gnomAD allele frequency < 1% (or missing)

### Tier 2: rare_moderate_deleterious
- PASS variants with MODERATE VEP impact (missense, in-frame indel)
- gnomAD allele frequency < 1% (or missing)
- At least one deleterious predictor hit (if vcfanno annotations available):
  - CADD PHRED >= 20 (SNV and/or indel tags, whichever are present)
  - REVEL >= 0.5
  - AlphaMissense "likely_pathogenic"
  - SpliceAI annotation present (presence check only — bcftools cannot parse the pipe-delimited delta scores; threshold filtering at >= 0.2 is done in step 23)
- Without vcfanno: all rare MODERATE variants included (same as step 23)
- Only annotation tags that exist in the VCF header are referenced — partial installs (e.g., CADD SNVs only) work correctly

### Tier 3: clinvar_pathogenic
- PASS variants with ClinVar pathogenic or likely_pathogenic (from VEP CLIN_SIG field)
- Excludes `conflicting_interpretations_of_pathogenicity` (the substring match `~"pathogenic"` would otherwise include these)
- No frequency filter (pathogenic variants can be common carriers)

All tiers are merged and deduplicated into a single prioritized VCF.

## Compound Heterozygote Detection

slivar's `compound-hets` command groups heterozygous variants by gene from the prioritized VCF and reports pairs that could form compound heterozygotes (two different damaging variants in the same gene, one from each parent). It requires a PED file (`--ped`) describing sample relationships. For singleton samples (no trio), the `--allow-non-trios` flag is required.

The command outputs VCF to stdout with `INFO/slivar_comphet` annotations linking partner variants. Each VCF record represents a unique variant; the `slivar_comphet` field lists all its compound-het partners (format: `sample/GENE/PAIR_ID/chrom/pos/ref/alt`, comma-separated). A gene with N variants produces up to C(N,2) pairs but only N VCF records. The script counts unique pair IDs from this field and exports a human-readable TSV with columns: GENE, CHROM, POS, REF, ALT, IMPACT, Consequence, GT -- sorted by gene so that compound-het partners appear in consecutive rows.

**Important:** With single-sample unphased data, these are *candidates only*. The two variants might be on the same haplotype (cis) rather than different haplotypes (trans). Trio data or read-backed phasing is needed to confirm true compound hets.

## Gene Constraint Enrichment

If `gnomad_v4.1_constraint.tsv` is available, the summary TSV is enriched with per-gene constraint metrics:

| Column | Description | Threshold |
|---|---|---|
| LOEUF | Loss-of-function observed/expected upper bound | < 0.35 = constrained |
| pLI | Probability of LoF intolerance | > 0.9 = constrained |
| mis_z | Missense Z-score | > 3.09 = constrained |
| CONSTRAINED | YES if LOEUF < 0.35 or pLI > 0.9 | Flag column |

Variants in constrained genes are more likely to be pathogenic -- these genes are under strong purifying selection against damaging variants.

## Output

| File | Description |
|---|---|
| `slivar/${SAMPLE}_prioritized.vcf.gz` | All prioritized variants (merged, deduplicated) |
| `slivar/${SAMPLE}_slivar_summary.tsv` | Human-readable table with gene constraint |
| `slivar/${SAMPLE}_compound_hets.vcf.gz` | Compound het candidate variants (VCF) |
| `slivar/${SAMPLE}_compound_hets.tsv` | Compound het candidates (human-readable TSV) |
| `slivar/${SAMPLE}_rare_high.vcf.gz` | Tier 1: HIGH impact |
| `slivar/${SAMPLE}_rare_moderate_del.vcf.gz` | Tier 2: MODERATE + deleterious |
| `slivar/${SAMPLE}_clinvar_path.vcf.gz` | Tier 3: ClinVar P/LP |

## Runtime

~5-10 minutes. Most time is spent on bcftools split-vep filtering.

## Interpretation

A typical 30X WGS sample produces:
- **rare_high:** 50-200 variants (loss-of-function in rare alleles)
- **rare_moderate_deleterious:** 200-1,000 variants (depends on predictor availability)
- **clinvar_pathogenic:** 0-10 variants (most are heterozygous carriers)
- **compound het candidates:** 1,000-2,000 pairs across 100-200 genes (combinatorial: N variants in a gene = N*(N-1)/2 pairs; most are false positives in unphased data)

Focus review on:
1. Variants in **constrained genes** (CONSTRAINED=YES)
2. **Homozygous** rare HIGH/MODERATE variants (potential recessive disease)
3. **Compound het pairs** in known disease genes
4. Any **ClinVar pathogenic** variants, especially in dominant disease genes

See [docs/interpreting-results.md](interpreting-results.md) for pathogenicity score thresholds.
