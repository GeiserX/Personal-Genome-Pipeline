# Step 2: Alignment (FASTQ to BAM)

## What This Does
Aligns raw sequencing reads against the GRCh38 human reference genome. Produces a sorted, indexed BAM file.

## Why
Alignment maps each 150bp sequencing read to its position in the human genome. Required for all downstream variant calling.

## Tools
- **minimap2** — fast aligner (preferred for WGS)
- **samtools** — sort + index the alignment

## Docker Images
- `staphb/samtools:1.20` (includes samtools)
- minimap2 is typically run natively or via `quay.io/biocontainers/minimap2`

## Prerequisites
- GRCh38 reference genome (`Homo_sapiens_assembly38.fasta`)
- minimap2 index (`.mmi` file, ~7GB, generated once)
- Paired-end FASTQ files

## Commands
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
REF=${GENOME_DIR}/reference/Homo_sapiens_assembly38.fasta

# Step 1: Create minimap2 index (one-time, ~30 min)
minimap2 -d ${GENOME_DIR}/reference/GRCh38.mmi $REF

# Step 2: Align + sort (4-8 hours for 30X WGS)
minimap2 -a -x map-hifi -t 16 \
  ${GENOME_DIR}/reference/GRCh38.mmi \
  ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
  ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz \
| samtools sort -@ 8 -o ${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Step 3: Index BAM
samtools index ${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Output: ~30-40GB BAM + ~9MB BAI index
```

## Resource Requirements
- CPU: 16+ cores recommended
- RAM: 16GB+ (minimap2 loads full index into memory)
- Disk: ~30-40GB per sample (BAM file)
- Time: 4-8 hours for 30X WGS

## Notes
- Use `-x map-hifi` for Illumina short reads (despite the name, works well)
- Alternative: `bwa-mem2` is equally valid but minimap2 is faster
- The BAM index (`.bai`) must always accompany the BAM file
