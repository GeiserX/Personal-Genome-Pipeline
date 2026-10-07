# Step 12: Mitochondrial Haplogroup

## What This Does
Determines maternal lineage ancestry from mitochondrial DNA variants. Also screens for pathogenic mtDNA mutations.

## Why
Mitochondrial haplogroup reveals deep maternal ancestry and can identify mtDNA disease variants. Some haplogroups have known health associations (e.g., longevity, metabolic traits).

## Tool
- **haplogrep3** (Medical University of Innsbruck)

## Docker Image
- `HAPLOGREP3_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

> The image is the Bioconda build of haplogrep3 3.2.2. It classifies without a network: the image test runs it with `--network none`. It replaced a digest-pinned build of 3.2.1 from a personal Docker Hub account; on the image test's chrM calls of the HG002 fixture both give H5a7 with quality 1.0000 and the same found and remaining polymorphisms.

## Command
```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 1: Extract chrM variants from VCF
docker run --rm -v ${GENOME_DIR}/${SAMPLE}/vcf:/genome/${SAMPLE}/vcf "${BCFTOOLS_IMAGE}" \
  bcftools view -r chrM /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz -Oz -o /genome/${SAMPLE}/vcf/${SAMPLE}_chrM.vcf.gz

# Step 2: Run haplogrep3
docker run --rm --network none -v ${GENOME_DIR}/${SAMPLE}:/genome/${SAMPLE} "${HAPLOGREP3_IMAGE}" \
  haplogrep3 classify \
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
