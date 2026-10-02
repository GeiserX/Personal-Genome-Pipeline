# Step 2: Alignment (FASTQ to BAM)

## What This Does
Aligns raw sequencing reads against the GRCh38 human reference genome. Produces a sorted, indexed BAM file.

## Why
Alignment maps each 150bp sequencing read to its position in the human genome. Required for all downstream variant calling.

## Tools
- **minimap2** — fast aligner (preferred for WGS)
- **samtools** — sort + index the alignment

## Docker Images
- `MINIMAP2_IMAGE` (minimap2 aligner)
- `SAMTOOLS_IMAGE` (samtools sort + index)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- GRCh38 reference genome (`${REF_FASTA}`, see [reference setup](00-reference-setup.md#the-reference-path-on-every-page))
- minimap2 index (`.mmi` file, ~7GB, generated once)
- Paired-end FASTQ files

## Commands
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
REF="${GENOME_DIR}/${REF_FASTA}"

# Step 1: Create minimap2 index (one-time, ~30 min)
minimap2 -d ${GENOME_DIR}/reference/GRCh38.mmi $REF

# Step 2: Align + sort (1-2 hours for 30X WGS)
minimap2 -a -x sr -t 16 \
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
- Time: 1-2 hours for 30X WGS

## Notes
- Use `-x sr` for Illumina short reads (short-read preset)
- Alternative: `bwa-mem2` is equally valid but minimap2 is faster (see `scripts/02a-alignment-bwamem2.sh`)
- The BAM index (`.bai`) must always accompany the BAM file
