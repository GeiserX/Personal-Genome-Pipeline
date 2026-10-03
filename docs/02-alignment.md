# Step 2: Alignment (FASTQ to BAM)

## What This Does
Aligns raw sequencing reads against the GRCh38 human reference genome and marks duplicate reads. Produces a sorted, duplicate-marked, indexed BAM file.

## Why
Alignment maps each 150bp sequencing read to its position in the human genome. Required for all downstream variant calling.

Duplicates are copies of one DNA fragment (PCR or optical copies). They are not independent evidence. DeepVariant copes with them, but GATK, FreeBayes, Octopus and the structural-variant and depth steps would count each copy as another read. `samtools markdup` sets the duplicate flag (0x400) on all copies but one, and those tools then skip them. No read is removed.

## Tools
- **minimap2**: fast aligner (preferred for WGS)
- **samtools**: fixmate, sort, markdup and index

## Docker Images
- `MINIMAP2_IMAGE` (minimap2 aligner)
- `SAMTOOLS_IMAGE` (samtools fixmate, sort, markdup + index)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- GRCh38 reference genome (`${REF_FASTA}`, see [reference setup](00-reference-setup.md#the-reference-path-on-every-page))
- minimap2 index (`.sr.mmi` file next to the reference, ~7GB, generated once by the script)
- Paired-end FASTQ files

## Commands
```bash
./scripts/02-alignment.sh your_sample        # THREADS=16 to use 16 CPUs (default 8)
```

What the script runs, written as plain commands:

```bash
REF_FASTA=reference/Homo_sapiens_assembly38.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
REF="${GENOME_DIR}/${REF_FASTA}"
MMI="${REF%.fasta}.sr.mmi"
OUT="${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam"

# Step 1: Create the minimap2 index with the short-read preset (one-time, ~30 min).
# The preset sets the k-mer and window size of the index (k21, w11 for -x sr);
# a plain `minimap2 -d` builds k15, w10, which does not match `-x sr` mapping.
minimap2 -x sr -d "$MMI" "$REF"

# Step 2: Align, mark duplicates and sort (1-2 hours for 30X WGS).
# fixmate -m adds the mate tags markdup needs; it reads the pairs minimap2
# writes next to each other, then sort orders by position for markdup.
minimap2 -a -x sr -t 16 \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
  "$MMI" \
  ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
  ${GENOME_DIR}/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz \
| samtools fixmate -u -m - - \
| samtools sort -u -@ 16 -m 1G - \
| samtools markdup -@ 16 - "${OUT%.bam}.tmp.bam"

# Step 3: Index, check, then rename into place
samtools index "${OUT%.bam}.tmp.bam"
samtools quickcheck "${OUT%.bam}.tmp.bam"
mv "${OUT%.bam}.tmp.bam" "$OUT" && mv "${OUT%.bam}.tmp.bam.bai" "${OUT}.bai"

# Output: ~80-120 GB BAM + ~9 MB BAI index
```

The script writes the BAM and the index under temporary names and renames them only after `samtools quickcheck` passes, so a run that is killed leaves no `${SAMPLE}_sorted.bam` behind. `run-all.sh` skips alignment only when the BAM, its `.bai` and quickcheck are all good. An index built by an older version (`reference/GRCh38.mmi`, default preset) is no longer used and can be deleted.

## Resource Requirements
- CPU: 16+ cores recommended (`THREADS`)
- RAM: 16GB+ (minimap2 loads full index into memory); samtools sort takes 1 GB per thread
- Disk and time: see [Hardware and storage requirements](hardware-requirements.md#runtime-per-step) (the BAM is about 80-120 GB). Sort spills go to `aligned/${SAMPLE}.sort_tmp/` in the sample directory and are removed at the end.

## Notes
- Use `-x sr` for Illumina short reads (short-read preset), for the index and for mapping
- Alternative: `bwa-mem2` is equally valid but minimap2 is faster (see `scripts/02a-alignment-bwamem2.sh`). It writes `aligned_bwamem2/` with the same fixmate, sort and markdup steps, piped, with no SAM file on disk. Its one-time index build needs about 90 GB of RAM for GRCh38 (28 GB per Gbp of reference, see [reference setup](00-reference-setup.md)); the script stops with that message when the build is killed for lack of memory.
- The BAM index (`.bai`) must always accompany the BAM file
