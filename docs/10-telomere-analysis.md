# Step 10: Telomere Length Estimation

## What This Does
Estimates telomere content from WGS BAM files by quantifying telomeric repeat reads (TTAGGG/CCCTAA). Provides a biological age proxy based on telomere attrition.

## Why
Telomere length correlates with cellular aging and is a biomarker for age-related disease risk. Shorter telomeres are associated with cardiovascular disease, cancer predisposition, and reduced lifespan. Comparing telomere content between individuals of similar age provides a relative aging benchmark.

## Tool
- **TelomereHunter** (German Cancer Research Center)

## Docker Image
```
lgalarno/telomerehunter:latest
```

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v ${GENOME_DIR}:/genome \
  lgalarno/telomerehunter:latest \
  telomerehunter \
    -ibt /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    -o /genome/${SAMPLE}/telomere/${SAMPLE} \
    -p ${SAMPLE}

# Output: telomere content report in ${GENOME_DIR}/${SAMPLE}/telomere/${SAMPLE}/
```

TelomereHunter uses `-ibt` (input BAM tumor) for a single-sample analysis. No `--tumor_only` flag is needed — when only `-ibt` is provided (without `-ibc` for a matched control BAM), TelomereHunter runs in single-sample mode automatically.

## Key Metric
- **`tel_content`** — GC-corrected telomeric reads per million mapped reads
- Higher values indicate longer/more abundant telomeres
- Compare between samples of known age for relative ranking

## Important Notes
- `--user root` is REQUIRED — Docker container cannot write output files without root permissions
- Short-read WGS (150bp reads) systematically underestimates true telomere length because reads cannot span long repetitive regions
- Results are useful as a **relative comparison** between samples, NOT as an absolute telomere length measurement
- Long-read sequencing (PacBio/ONT) provides more accurate telomere length if absolute values are needed
