# Step 30: Annotation Enrichment (vcfanno)

Adds pathogenicity scores (CADD, SpliceAI, REVEL, AlphaMissense) to VEP-annotated VCFs using vcfanno's TOML-driven annotation engine.

## Why

VEP (step 13) provides functional impact predictions, but thorough variant interpretation benefits from additional pathogenicity scores. CADD captures deleteriousness across coding and non-coding variants; SpliceAI finds cryptic splice-disrupting variants missed by rule-based methods; REVEL and AlphaMissense score missense variants using ensemble and structural approaches respectively.

vcfanno annotates in a single streaming pass per database, making it much faster than re-running VEP with plugins.

## Prerequisites

- VEP-annotated VCF from step 13 (`${SAMPLE}_vep.vcf.gz`)
- One or more annotation databases in `${GENOME_DIR}/annotations/` (all optional)

See [docs/00-reference-setup.md](00-reference-setup.md) for download instructions.

## Annotation Databases

| Database | File | Size | License |
|---|---|---|---|
| CADD v1.7 SNVs | `whole_genome_SNVs.tsv.gz` + `.tbi` | ~81.5 GB | Non-commercial |
| CADD v1.7 indels | `gnomad.genomes.r4.0.indel.tsv.gz` + `.tbi` | ~1.2 GB | Non-commercial |
| SpliceAI SNVs | `spliceai_scores.raw.snv.hg38.vcf.gz` + `.tbi` | ~16 GB | Apache 2.0 |
| SpliceAI indels | `spliceai_scores.raw.indel.hg38.vcf.gz` + `.tbi` | ~4 GB | Apache 2.0 |
| REVEL v1.3 | `revel_grch38.tsv.gz` + `.tbi` | ~526 MB | Free for research |
| AlphaMissense | `AlphaMissense_hg38.tsv.gz` + `.tbi` | ~613 MB | CC BY-NC-SA 4.0 |

All databases are optional. The script detects which files are present and annotates accordingly. Missing databases are silently skipped.

## Chromosome Naming Mismatch

CADD files use bare chromosome names (`1`, `2`, `3`) while the pipeline VCFs and other databases use chr-prefixed names (`chr1`, `chr2`, `chr3`). The script handles this with a two-pass approach:

1. **Pass 1 (CADD):** Strip `chr` prefix from VCF, annotate with CADD, re-add `chr` prefix
2. **Pass 2 (others):** Annotate with SpliceAI, REVEL, AlphaMissense (all chr-prefixed)

If only chr-prefixed databases are present (no CADD), a single pass is used.

## Docker Image

```
quay.io/biocontainers/vcfanno:0.3.7--he881be0_0
```

Also uses the bcftools image for bgzip/tabix/chr renaming operations.

## Usage

```bash
export GENOME_DIR=/path/to/your/data
./scripts/30-vcfanno.sh your_sample
```

## Output

| File | Description |
|---|---|
| `vep/${SAMPLE}_annotated.vcf.gz` | VCF with CADD_PHRED, SpliceAI, REVEL, AM_pathogenicity, AM_class in INFO |
| `vep/${SAMPLE}_annotated.vcf.gz.tbi` | Tabix index |

## Runtime

~5-15 minutes depending on which databases are present. The two-pass CADD approach adds ~2-3 minutes for the chromosome renaming steps.

## INFO Fields Added

| Field | Source | Description |
|---|---|---|
| `CADD_PHRED` | CADD v1.7 | PHRED-scaled deleteriousness score (higher = more deleterious) |
| `CADD_PHRED_indel` | CADD v1.7 | CADD score for indels specifically |
| `SpliceAI` | SpliceAI | Splice impact prediction (delta scores for AG/AL/DG/DL) |
| `REVEL` | REVEL v1.3 | Ensemble missense pathogenicity score (0-1) |
| `AM_pathogenicity` | AlphaMissense | Structure-informed missense pathogenicity (0-1) |
| `AM_class` | AlphaMissense | Classification: benign, ambiguous, or likely_pathogenic |

## Interpretation

See [docs/interpreting-results.md](interpreting-results.md) for score thresholds and clinical interpretation guidance.

## Query Examples

```bash
# Variants with CADD PHRED >= 20 (top 1% most deleterious)
bcftools view -i 'INFO/CADD_PHRED>=20' ${SAMPLE}_annotated.vcf.gz | head

# AlphaMissense likely pathogenic variants
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AM_class\n' \
  -i 'INFO/AM_class="likely_pathogenic"' ${SAMPLE}_annotated.vcf.gz

# High REVEL score missense variants (ClinGen moderate evidence)
bcftools view -i 'INFO/REVEL>=0.644' ${SAMPLE}_annotated.vcf.gz
```

## Downstream Steps

- **Step 23 (Clinical Filter):** Automatically uses the enriched VCF if available, adding CADD/SpliceAI/REVEL/AlphaMissense filter tiers
- **Step 31 (slivar):** Uses vcfanno annotations for deleterious predictor filtering in the rare_moderate tier
