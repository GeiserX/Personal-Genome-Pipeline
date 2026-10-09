# Step 19: Structural Variant Calling with Delly

## What This Does
Third structural variant caller — combines paired-end, split-read, and read-depth signals for comprehensive SV detection including deletions, duplications, inversions, translocations, and insertions.

## Why
Using multiple SV callers and intersecting their results dramatically reduces false positives:
- **Manta** (step 4): Fast, sensitive for smaller SVs and indels
- **CNVpytor** (step 18): Best for large CNVs via read-depth only
- **Delly**: Most balanced — uses all three signal types, especially strong for inversions and translocations

SVs called by 2+ callers have lower false-positive rates than single-caller calls. Multi-caller intersection is a common strategy in WGS pipelines, though dedicated tools like SURVIVOR or Jasmine provide more precise breakpoint-aware merging than simple position overlap.

## Tool
- **Delly** (Rausch et al., Bioinformatics 2012)

## Docker Image
- `DELLY_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/19-delly.sh your_sample
```

The script passes Delly's GRCh38 exclude map (`-x`): telomeres, centromeres and every contig beyond chr1-22, X, Y and M (on the default no-ALT reference, the unplaced scaffolds and `chrEBV`; the map also names the ALT and decoy contigs of a full reference). `setup.sh` installs it from a pinned commit of the Delly repository as `reference/delly_human.hg38.excl.tsv` (see [reference setup](00-reference-setup.md#small-pinned-data-files)). Without it Delly spends hours in those regions and calls artefacts there; the script then runs without `-x` and says so.

What the script runs. Delly 2.3.0 renamed the short-read caller from `delly call` to `delly sr`; the pinned 2.7.0 answers `Unrecognized command` to `delly call`.

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
# SV calling (all SV types)
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  "${DELLY_IMAGE}" \
  delly sr \
    -g "/genome/${REF_FASTA}" \
    -x /genome/reference/delly_human.hg38.excl.tsv \
    -o /genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf \
    /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam

# Convert BCF to VCF for downstream tools
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  "${BCFTOOLS_IMAGE}" \
  bcftools view \
    /genome/${SAMPLE}/delly/${SAMPLE}_sv.bcf \
    -Oz -o /genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz

# Index
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  "${BCFTOOLS_IMAGE}" \
  bcftools index -t \
    /genome/${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz
```

## Optional: Dedicated CNV Calling
Delly also has a dedicated CNV mode using read-depth only (similar to CNVpytor):
```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  "${DELLY_IMAGE}" \
  delly cnv \
    -g "/genome/${REF_FASTA}" \
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
- For consensus SV calling, use SURVIVOR or bcftools to merge calls from Manta + Delly + CNVpytor.
- Delly is the most accurate caller for inversions and balanced translocations.
- The `delly cnv` mode is optional if you already run CNVpytor — it provides similar depth-based CNV calls.
- Can be run in parallel with Manta and CNVpytor (all independent after alignment).
