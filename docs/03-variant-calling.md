# Step 3: Variant Calling (BAM to VCF)

## What This Does
Identifies all positions where the sample's DNA differs from the reference genome: SNPs (single nucleotide changes) and small indels (insertions/deletions <50bp).

## Why
The VCF file is the foundation for ALL downstream analyses: ClinVar screening, pharmacogenomics, PRS, ROH, etc.

## Tool
- **DeepVariant** v1.10.0 — Google's deep learning variant caller (state-of-the-art accuracy)

## Docker Image
- `DEEPVARIANT_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sorted, indexed BAM file
- GRCh38 reference genome (FASTA + FAI)

## Command
```bash
source versions.env   # from the repository root
REF_FASTA=reference/Homo_sapiens_assembly38.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  --cpus 8 --memory 32g \
  -v ${GENOME_DIR}:/genome \
  "${DEEPVARIANT_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/genome/${REF_FASTA}" \
    --reads=/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --output_vcf=/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    --sample_name="${SAMPLE}" \
    --num_shards=8

# For WES data, use MODEL_TYPE=WES:
# MODEL_TYPE=WES ./scripts/03-deepvariant.sh your_sample

# To call only some regions, set INTERVALS (space-separated, passed to --regions):
# INTERVALS="chr20:10000001-10500000" ./scripts/03-deepvariant.sh your_sample

# Output: ~93MB VCF with ~5.5M total variants (~4.6M PASS)
```

## Resource Requirements
- CPU: the script uses 8 (`--cpus 8`, `--num_shards=8`); more shards scale well if you run the command by hand on more cores
- RAM: 32GB recommended
- GPU: optional, and only `call_variants` uses it (see [troubleshooting](troubleshooting.md#step-3-deepvariant-gpu-acceleration-not-worth-it))
- Time: see [Hardware and storage requirements](hardware-requirements.md#runtime-per-step)

## Output Interpretation
- **PASS** variants: high-confidence calls (~4.6M per 30X WGS sample)
- **RefCall**: site looks like reference (not a variant)
- **LowQual**: low confidence
- GQ (Genotype Quality): higher = more confident
- DP (Read Depth): typical 25-35x for a 30X WGS sample

## Notes
- DeepVariant is the gold standard for SNPs/indels but does NOT detect:
  - Structural variants >50bp (use Manta, step 4)
  - Repeat expansions (use ExpansionHunter, step 9)
  - Copy number variants (use CNVpytor or GATK gCNV)
- bcftools can also call variants but is significantly less accurate than DeepVariant
