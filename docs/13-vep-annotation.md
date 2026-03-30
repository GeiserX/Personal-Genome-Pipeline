# Step 13: Variant Effect Predictor (VEP) Annotation

## What This Does
Annotates every variant in the VCF with gene name, consequence type, predicted impact, pathogenicity scores (SIFT, PolyPhen), and population allele frequencies. This is the most comprehensive single annotation step in the pipeline.

## Why
Raw VCF variants are just genomic coordinates and genotypes. VEP transforms them into biologically interpretable annotations — which gene is affected, what the functional consequence is, how rare the variant is in the population, and whether it is predicted damaging.

## Tool
- **Ensembl VEP** release 112 (European Bioinformatics Institute)

## Docker Image
```
ensemblorg/ensembl-vep:release_112.0
```

## Prerequisites
- Offline VEP cache must be downloaded first (see step 00-reference-setup)
- Cache size: ~26 GB for GRCh38 homo_sapiens

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v ${GENOME_DIR}/vep_cache:/vep_cache \
  ensemblorg/ensembl-vep:release_112.0 \
  vep \
    --input_file /data/${SAMPLE}.vcf.gz \
    --output_file /data/${SAMPLE}_vep.vcf \
    --vcf \
    --offline \
    --cache \
    --dir_cache /vep_cache \
    --assembly GRCh38 \
    --fork 4 \
    --sift b \
    --polyphen b \
    --af \
    --af_gnomade \
    --canonical \
    --symbol \
    --force_overwrite

# Output: VCF with CSQ INFO field containing all annotations
```

## Output Format
- Default: VCF with `CSQ` INFO field (pipe-delimited sub-fields)
- Alternative: add `--tab` instead of `--vcf` for tab-delimited output (easier to parse manually)
- Key CSQ sub-fields: `SYMBOL`, `Consequence`, `IMPACT`, `SIFT`, `PolyPhen`, `gnomADe_AF`, `CANONICAL`

## Filtering for Clinical Relevance
After annotation, filter to actionable variants:
```bash
# Extract HIGH and MODERATE impact variants with AF <1%
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  staphb/bcftools:1.21 \
  bcftools view -i 'INFO/CSQ[*] ~ "HIGH" || INFO/CSQ[*] ~ "MODERATE"' \
    /data/${SAMPLE}_vep.vcf \
    > /data/${SAMPLE}_vep_filtered.vcf
```

## Important Notes
- Full WGS annotation takes **2-4 hours** depending on CPU and variant count (~5M variants)
- `--fork 4` enables parallelism — increase if more cores are available
- `--canonical` restricts to canonical transcripts, reducing redundant annotations per variant
- `--sift b` and `--polyphen b` output both prediction and score (e.g., `deleterious(0.01)`)
- `--af_gnomade` adds gnomAD exomes frequency — use `--af_gnomad` (without 'e') for gnomAD genomes
- Variants with IMPACT=HIGH and gnomAD AF <0.01 (or absent) are the highest priority for clinical review
- VEP does NOT assess variant pathogenicity in ClinVar context — combine with ClinVar annotation for full picture
