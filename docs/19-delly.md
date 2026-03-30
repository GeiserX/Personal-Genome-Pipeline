# Step 19: Structural Variant Calling with Delly

## What This Does
Third structural variant caller — combines paired-end, split-read, and read-depth signals for comprehensive SV detection including deletions, duplications, inversions, translocations, and insertions.

## Why
Using multiple SV callers and intersecting their results dramatically reduces false positives:
- **Manta** (step 4): Fast, sensitive for smaller SVs and indels
- **CNVnator** (step 18): Best for large CNVs via read-depth only
- **Delly**: Most balanced — uses all three signal types, especially strong for inversions and translocations

SVs called by 2+ callers are high-confidence. This consensus approach is standard in clinical WGS pipelines.

## Tool
- **Delly** (Rausch et al., Bioinformatics 2012)

## Docker Image
```
quay.io/biocontainers/delly:1.7.3--hd6466ae_0
```

## Command
```bash
# SV calling (all SV types)
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/delly:1.7.3--hd6466ae_0 \
  delly call \
    -g /genome/reference/Homo_sapiens_assembly38.fasta \
    -o /genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf \
    /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Convert BCF to VCF for downstream tools
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/bcftools:1.21 \
  bcftools view \
    /genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf \
    -Oz -o /genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz

# Index
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/bcftools:1.21 \
  bcftools index -t \
    /genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz
```

## Optional: Dedicated CNV Calling
Delly also has a dedicated CNV mode using read-depth only (similar to CNVnator):
```bash
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/delly:1.7.3--hd6466ae_0 \
  delly cnv \
    -g /genome/reference/Homo_sapiens_assembly38.fasta \
    -o /genome/${SAMPLE}/delly/${SAMPLE}_cnv.bcf \
    /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

## Output
- `${SAMPLE}_sv.bcf` / `${SAMPLE}_sv.vcf.gz` — SV calls in VCF format
- Each SV has type (DEL, DUP, INV, BND, INS), quality, genotype, and supporting read counts

## Filtering
```bash
# Keep only PASS variants
bcftools view -f PASS ${SAMPLE}_sv.vcf.gz

# Filter by SV type
bcftools view -i 'INFO/SVTYPE="DEL"' ${SAMPLE}_sv.vcf.gz
bcftools view -i 'INFO/SVTYPE="INV"' ${SAMPLE}_sv.vcf.gz
```

## Runtime
~2-4 hours per 30X WGS genome.

## Notes
- Delly outputs BCF by default (not VCF). Convert with `bcftools view` for compatibility.
- For consensus SV calling, use SURVIVOR or bcftools to merge calls from Manta + Delly + CNVnator.
- Delly is the most accurate caller for inversions and balanced translocations.
- The `delly cnv` mode is optional if you already run CNVnator — it provides similar depth-based CNV calls.
- Can be run in parallel with Manta and CNVnator (all independent after alignment).
