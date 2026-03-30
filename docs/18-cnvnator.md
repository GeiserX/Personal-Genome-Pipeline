# Step 18: Depth-Based CNV Calling with CNVnator

## What This Does
Detects copy number variants (CNVs) using read-depth analysis — complementary to Manta's paired-end/split-read approach. Especially effective for large CNVs (>1 kb) that Manta may miss.

## Why
Manta (step 4) detects SVs from discordant read pairs and split reads, which works well for balanced SVs (inversions, translocations) and smaller deletions/duplications. CNVnator uses read-depth signal only, making it better for:
- Large duplications/deletions (>10 kb)
- Tandem duplications where breakpoints are ambiguous
- Validating Manta calls with orthogonal evidence

## Tool
- **CNVnator** (Abyzov et al., Genome Research 2011)

## Docker Image
```
quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11
```

## Command
```bash
# Step 1: Extract read mapping from BAM
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 \
  cnvnator \
    -root /genome/${SAMPLE}/cnvnator/${SAMPLE}.root \
    -tree /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Step 2: Generate read-depth histogram
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 \
  cnvnator \
    -root /genome/${SAMPLE}/cnvnator/${SAMPLE}.root \
    -his 1000 \
    -fasta /genome/reference/Homo_sapiens_assembly38.fasta

# Step 3: Statistics
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 \
  cnvnator \
    -root /genome/${SAMPLE}/cnvnator/${SAMPLE}.root \
    -stat 1000

# Step 4: Partition
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 \
  cnvnator \
    -root /genome/${SAMPLE}/cnvnator/${SAMPLE}.root \
    -partition 1000

# Step 5: Call CNVs
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/cnvnator:0.4.1--py312h99c8fb2_11 \
  cnvnator \
    -root /genome/${SAMPLE}/cnvnator/${SAMPLE}.root \
    -call 1000 \
    > ${GENOME_DIR}/${SAMPLE}/cnvnator/${SAMPLE}_cnvs.txt
```

## Bin Size
The `1000` parameter is the bin size in base pairs. Use:
- `1000` for 30X WGS (recommended)
- `500` for higher coverage (>50X)
- `100` for targeted/exome data

## Output
- `${SAMPLE}.root` — ROOT file with read-depth data (intermediate, can be deleted)
- `${SAMPLE}_cnvs.txt` — Tab-separated CNV calls with columns:
  - Type (deletion/duplication)
  - Coordinates (chr:start-end)
  - Size
  - Normalized read depth
  - e-value (statistical significance)
  - q0 (fraction of reads with mapping quality 0)

## Filtering
```bash
# Keep only significant CNVs (e-value < 0.01, size > 1kb)
awk '$5 < 0.01 && $3 > 1000' ${SAMPLE}_cnvs.txt
```

## Runtime
~2-4 hours per 30X WGS genome.

## Notes
- Run AFTER alignment (step 2). Independent of Manta — can run in parallel.
- CNVnator and Manta overlap on ~60-70% of true CNVs. Variants called by both are high-confidence.
- The ROOT file format is from CERN's ROOT framework — the Docker image includes all dependencies.
- For maximum sensitivity, intersect CNVnator calls with Manta and Delly (step 19) for a consensus call set.
