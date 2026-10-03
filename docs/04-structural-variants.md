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
REF_FASTA=reference/Homo_sapiens_assembly38.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

THREADS=8

# Step 1: Configure Manta
docker run --rm \
  --cpus ${THREADS} --memory 16g \
  -v ${GENOME_DIR}:/genome \
  "${MANTA_IMAGE}" \
  configManta.py \
    --bam /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    --referenceFasta "/genome/${REF_FASTA}" \
    --runDir /genome/${SAMPLE}/manta

# Step 2: Run Manta
docker run --rm \
  --cpus ${THREADS} --memory 16g \
  -v ${GENOME_DIR}:/genome \
  "${MANTA_IMAGE}" \
  /genome/${SAMPLE}/manta/runWorkflow.py -j ${THREADS}

# Step 3: Turn inversion breakend pairs into SVTYPE=INV records with the
# convertInversion.py, samtools, bgzip and tabix inside the Manta image.
# Manta's own file is kept as diploidSV.raw.vcf.gz.
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  "${MANTA_IMAGE}" \
  bash -c 'set -euo pipefail
    L="$(dirname "$(readlink -f "$(command -v configManta.py)")")/../libexec"
    cd "$2"
    mv diploidSV.vcf.gz diploidSV.raw.vcf.gz
    mv diploidSV.vcf.gz.tbi diploidSV.raw.vcf.gz.tbi
    "$L/convertInversion.py" "$L/samtools" "$1" diploidSV.raw.vcf.gz | "$L/bgzip" -c > diploidSV.vcf.gz
    "$L/tabix" -f -p vcf diploidSV.vcf.gz' \
  _ "/genome/${REF_FASTA}" "/genome/${SAMPLE}/manta/results/variants"

# Output: diploidSV.vcf.gz (~7-9K structural variants)
```

`scripts/04-manta.sh` runs the same three steps. It reads the CPU count from `THREADS` (default 8) and the BAM from `${SAMPLE}/${ALIGN_DIR}/` (default `aligned`). With `MANTA_CALL_REGIONS` set to a bgzipped BED under `GENOME_DIR` (its `.tbi` beside it), Manta calls only those regions (`configManta.py --callRegions`), for example chr1-22, X and Y without the ALT and decoy contigs. When Manta reported no inversion, the script copies `diploidSV.raw.vcf.gz` to `diploidSV.vcf.gz` instead of starting a container, and says so. Run on a folder that holds Manta's results from before the inversion step, it converts them without calling again.

## Output
- `results/variants/diploidSV.vcf.gz` — main output (all SV calls), with each inversion as one `SVTYPE=INV` record; steps 05, 15 and 22 read this file
- `results/variants/diploidSV.raw.vcf.gz` — the same calls as Manta wrote them, where an inversion is a pair of breakend (`SVTYPE=BND`) records
- `results/variants/candidateSV.vcf.gz` — unfiltered candidates
- `results/variants/candidateSmallIndels.vcf.gz` — small indels

## Important Notes
- Raw Manta output is UNFILTERED — most calls are benign
- **Must run AnnotSV (step 5)** to classify pathogenicity
- SVs >5MB in short-read WGS are usually artifacts from segmental duplications
- Typical results: 7,000-9,000 SVs per 30X WGS sample
