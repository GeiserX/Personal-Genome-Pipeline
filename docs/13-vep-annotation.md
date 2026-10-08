# Step 13: Variant Effect Predictor (VEP) Annotation

## What This Does
Annotates every variant in the VCF with gene name, consequence type, predicted impact, pathogenicity scores (SIFT, PolyPhen), population allele frequencies (gnomAD), ClinVar significance, and more. This is the most comprehensive single annotation step in the pipeline.

## Why
Raw VCF variants are just genomic coordinates and genotypes. VEP transforms them into biologically interpretable annotations — which gene is affected, what the functional consequence is, how rare the variant is in the population, and whether it is predicted damaging.

## Tool
- **Ensembl VEP** release 116 (European Bioinformatics Institute)

## Docker Image
- `VEP_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- The VCF from step 3 **and its `.tbi`**: the step refuses a VCF without its index, which may be half written
- Offline VEP cache: step 13 downloads, checks and unpacks it the first time it runs (see step 00-reference-setup)
- Cache size: see [Hardware and storage requirements](hardware-requirements.md#shared-reference-data-one-time) (about 26 GB to download, 30 GB unpacked)

## Command
```bash
export GENOME_DIR=/path/to/your/data
THREADS=8 ./scripts/13-vep-annotation.sh your_sample
```

`THREADS` (default 8) sets the container's CPUs and VEP's `--fork`; memory is 2 GB per fork, 8 GB at least. What the script runs, with the same annotation fields as the Nextflow VEP module:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see "The reference path on every page" in 00-reference-setup.md

docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOME_DIR}:/genome \
  -v ${GENOME_DIR}/vep_cache:/opt/vep/.vep \
  "${VEP_IMAGE}" \
  vep \
    --input_file /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    -o /genome/${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz \
    --vcf \
    --compress_output bgzip \
    --cache \
    --cache_version 116 \
    --dir_cache /opt/vep/.vep \
    --offline \
    --assembly GRCh38 \
    --fasta /genome/${REF_FASTA} \
    --everything \
    --force_overwrite \
    --fork 8 \
    --custom file=/genome/clinvar/clinvar_pathogenic_chr.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
```

The last line is added when `clinvar/clinvar_pathogenic_chr.vcf.gz` is installed (`setup.sh`): it is the ClinVar file step 6 screens against, and VEP copies its `CLNSIG`, `CLNREVSTAT` and `CLNDN` for each exact match into the CSQ fields `ClinVar_CLNSIG`, `ClinVar_CLNREVSTAT` and `ClinVar_CLNDN`. Step 23's ClinVar tier reads `ClinVar_CLNSIG` before VEP's own `CLIN_SIG`, which comes from the cache release, so a ClinVar refresh reaches the tier after this step runs again. The Nextflow VEP module adds the same `--custom` file when `--clinvar` is set.

The script writes VEP's output under a temporary name, renames it to `${SAMPLE}_vep.vcf.gz` only when VEP succeeded, and indexes it.

## Output
- `${SAMPLE}/vep/${SAMPLE}_vep.vcf.gz` (+ `.tbi`) — the VCF with the `CSQ` INFO field, bgzip-compressed
- `${SAMPLE}/vep/${SAMPLE}_vep_summary.html` and `${SAMPLE}_vep_warnings.txt` — VEP's run statistics and warnings

A finished run removes the files built from an older annotation: vcfanno's `${SAMPLE}_annotated.vcf.gz` (+ `.tbi`, step 30) and the uncompressed `${SAMPLE}_vep.vcf` earlier versions wrote. Steps 30, 23 and 31 prefer those files when they exist, so after a VEP or cache update they would otherwise keep reading the old annotation; rerun step 30 after step 13.

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
- `--fork` follows `THREADS`; each fork loads its own copy of the cache index
- `--everything` replaces individual flags (`--sift b`, `--polyphen b`, `--canonical`, `--af_gnomade`, etc.) with a single comprehensive flag
- `--dir_cache /opt/vep/.vep`: the cache is mounted there, not in the home directory VEP looks in by default
- `--fasta`: running `--offline` without a FASTA file disables HGVS notation (`INFO: Disabling --hgvs`), so the script always passes the reference ([`REF_FASTA`](00-reference-setup.md#the-reference-path-on-every-page) is the reference path)
- VEP does NOT assess variant pathogenicity in ClinVar context — combine with step 6 (ClinVar screen) for full picture
- **Upgrading from an older release:** VEP reads the cache directory named after the release passed in `--cache_version` (`homo_sapiens/116_GRCh38/` for the pinned release 116). When that directory is missing, VEP stops with an error instead of annotating from an older cache such as `112_GRCh38/`; step 13 downloads the matching cache on its first run.
