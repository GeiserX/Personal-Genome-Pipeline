# Step 12: Mitochondrial Haplogroup

## What This Does
Determines maternal lineage ancestry from mitochondrial DNA variants. Also screens for pathogenic mtDNA mutations.

## Why
Mitochondrial haplogroup reveals deep maternal ancestry and can identify mtDNA disease variants. Some haplogroups have known health associations (e.g., longevity, metabolic traits).

## Tool
- **haplogrep3** (Medical University of Innsbruck)

## Docker Image
```
jtb114/haplogrep3
```

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 1: Extract chrM variants from VCF
docker run --rm -v ${GENOME_DIR}/${SAMPLE}/vcf:/genome/${SAMPLE}/vcf staphb/bcftools:1.21 \
  bcftools view -r chrM /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz -Oz -o /genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz

# Step 2: Run haplogrep3
docker run --rm -v ${GENOME_DIR}/${SAMPLE}:/genome/${SAMPLE} jtb114/haplogrep3 \
  classify \
    --tree phylotree-fu-rcrs@1.2 \
    --input /genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz \
    --output /genome/${SAMPLE}/mito/${SAMPLE}_haplogroup.txt \
    --extend-report

# Output: haplogroup classification with quality score
```

## Interpretation
- Common European haplogroups: H, U, J, T, K, V, W, X
- Output includes quality score (0-1): >0.9 = high confidence
- Discordant variants may indicate heteroplasmy (mixture of mtDNA types)
