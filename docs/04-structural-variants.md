# Step 4: Structural Variant Calling (Manta)

## What This Does
Detects large DNA changes (>50bp) that DeepVariant misses: deletions, duplications, inversions, translocations, and insertions.

## Why
Structural variants cause ~25% of all genetic disease but are invisible to standard SNP/indel callers.

## Tool
- **Manta** (Illumina) — structural variant and indel caller

## Docker Image
```
quay.io/biocontainers/manta:1.6.0--h9ee0642_2
```

## Command
```bash
SAMPLE=your_sample
GENOMA_DIR=/path/to/genome/data

# Step 1: Configure Manta
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOMA_DIR}:/genoma \
  quay.io/biocontainers/manta:1.6.0--h9ee0642_2 \
  configManta.py \
    --bam /genoma/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --referenceFasta /genoma/reference/Homo_sapiens_assembly38.fasta \
    --runDir /genoma/${SAMPLE}/manta

# Step 2: Run Manta
docker run --rm \
  --cpus 8 --memory 16g \
  -v ${GENOMA_DIR}:/genoma \
  quay.io/biocontainers/manta:1.6.0--h9ee0642_2 \
  /genoma/${SAMPLE}/manta/runWorkflow.py -j 8

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
