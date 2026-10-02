# Step 4: Structural Variant Calling (Manta)

## What This Does
Detects large DNA changes (>50bp) that DeepVariant misses: deletions, duplications, inversions, translocations, and insertions.

## Why
Structural variants cause ~25% of all genetic disease but are invisible to standard SNP/indel callers.

## Tool
- **Manta** (Illumina) — structural variant and indel caller

## Docker Image
- `MANTA_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 1: Configure Manta
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOME_DIR}:/genome \
  "${MANTA_IMAGE}" \
  configManta.py \
    --bam /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --referenceFasta "/genome/${REF_FASTA}" \
    --runDir /genome/${SAMPLE}/manta

# Step 2: Run Manta
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOME_DIR}:/genome \
  "${MANTA_IMAGE}" \
  /genome/${SAMPLE}/manta/runWorkflow.py -j 8

# Output: diploidSV.vcf.gz (~7-9K structural variants)
```

## Output
- `results/variants/diploidSV.vcf.gz` — main output (all SV calls)
- `results/variants/candidateSV.vcf.gz` — unfiltered candidates
- `results/variants/candidateSmallIndels.vcf.gz` — small indels

## Important Notes
- Raw Manta output is UNFILTERED — most calls are benign
- **Must run AnnotSV (step 5)** to classify pathogenicity
- SVs >5MB in short-read WGS are usually artifacts from segmental duplications
- Typical results: 7,000-9,000 SVs per 30X WGS sample
