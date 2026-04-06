# Step 13: Variant Effect Predictor (VEP) Annotation

## What This Does
Annotates every variant in the VCF with gene name, consequence type, predicted impact, pathogenicity scores (SIFT, PolyPhen), population allele frequencies (gnomAD), ClinVar significance, and more. This is the most comprehensive single annotation step in the pipeline.

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
- Cache size: ~17 GB for GRCh38 homo_sapiens

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  --cpus 4 --memory 8g \
  --user root \
  -v ${GENOME_DIR}:/genome \
  -v ${GENOME_DIR}/vep_cache:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_112.0 \
  vep \
    --input_file /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    --output_file /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf \
    --vcf \
    --cache \
    --dir_cache /opt/vep/.vep \
    --offline \
    --assembly GRCh38 \
    --everything \
    --force_overwrite \
    --fork 4

# Output: VCF with CSQ INFO field containing all annotations
```

## Output Format
- Default: VCF with `CSQ` INFO field (pipe-delimited sub-fields)
- Alternative: add `--tab` instead of `--vcf` for tab-delimited output (easier to parse manually)
- The `--everything` flag enables all available annotations including:
  `SYMBOL`, `Consequence`, `IMPACT`, `SIFT`, `PolyPhen`, `gnomADe_AF`, `gnomADg_AF`, `MAX_AF`, `CLIN_SIG`, `CANONICAL`, `MANE_SELECT`, `BIOTYPE`, `Regulatory`, and many more

## Filtering for Clinical Relevance
After annotation, use step 23 (clinical filter) which automatically detects available CSQ fields and filters accordingly:
- HIGH impact variants (stop-gain, frameshift, splice)
- Rare MODERATE variants (gnomAD AF < 1%)
- ClinVar pathogenic/likely pathogenic hits

## Important Notes
- Full WGS annotation takes **2-4 hours** depending on CPU and variant count (~5M variants)
- `--fork 4` enables parallelism — increase if more cores are available
- `--everything` replaces individual flags (`--sift b`, `--polyphen b`, `--canonical`, `--af_gnomade`, etc.) with a single comprehensive flag
- `--dir_cache /opt/vep/.vep` is required when running as `--user root` (VEP looks in `/root/.vep` by default)
- Running `--offline` without a FASTA file disables HGVS notation (`INFO: Disabling --hgvs`). Add `--fasta /genome/reference/Homo_sapiens_assembly38.fasta` if HGVS is needed
- VEP does NOT assess variant pathogenicity in ClinVar context — combine with step 6 (ClinVar screen) for full picture
