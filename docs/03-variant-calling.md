# Step 3: Variant Calling (BAM to VCF)

## What This Does
Identifies all positions where the sample's DNA differs from the reference genome: SNPs (single nucleotide changes) and small indels (insertions/deletions <50bp).

## Why
The VCF file is the foundation for ALL downstream analyses: ClinVar screening, pharmacogenomics, PRS, ROH, etc.

## Tool
- **DeepVariant** v1.6.0 — Google's deep learning variant caller (state-of-the-art accuracy)

## Docker Image
```
google/deepvariant:1.6.0
```

## Prerequisites
- Sorted, indexed BAM file
- GRCh38 reference genome (FASTA + FAI)

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  --cpus 16 --memory 32g \
  -v ${GENOME_DIR}:/genome \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=/genome/reference/Homo_sapiens_assembly38.fasta \
    --reads=/genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --output_vcf=/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    --num_shards=16

# For WES data, use MODEL_TYPE=WES:
# MODEL_TYPE=WES ./scripts/03-deepvariant.sh your_sample

# Output: ~93MB VCF with ~5.5M total variants (~4.6M PASS)
```

## Resource Requirements
- CPU: 16+ cores (scales well with --num_shards)
- RAM: 32GB recommended
- GPU: Optional (significantly faster with NVIDIA GPU)
- Time: 4-8 hours on CPU, 1-2 hours with GPU

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
  - Copy number variants (use CNVnator or GATK gCNV)
- bcftools can also call variants but is significantly less accurate than DeepVariant
