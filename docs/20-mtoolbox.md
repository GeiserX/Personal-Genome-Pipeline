# Step 20: Mitochondrial Variant Calling and Heteroplasmy Detection

## What This Does
Calls mitochondrial DNA variants with heteroplasmy fractions — detecting variants present in only a fraction of your mitochondrial copies. Uses GATK Mutect2 in mitochondrial mode.

## Why
Step 12 (haplogrep3) assigns your mitochondrial haplogroup from chrM variants already in the main VCF. This step goes deeper:
- **Heteroplasmy detection**: Identifies variants present in only a fraction of mtDNA copies (clinically important for mitochondrial diseases)
- **Dedicated mitochondrial calling**: Mutect2's mitochondrial mode handles the unique properties of mtDNA (high copy number, circular genome, no recombination)
- **Somatic-grade sensitivity**: Detects variants at allele fractions as low as 1-3%

## Tool
- **GATK Mutect2** (Broad Institute) in `--mitochondria-mode`

> **Note:** This step was originally planned for MToolBox, but no working Docker image exists for MToolBox (see [lessons-learned.md](lessons-learned.md#mtoolbox-no-working-docker-image-exists)). GATK Mutect2 is the standard clinical alternative. The script is still called `scripts/20-mtoolbox.sh`; the name is historical.

## Docker Image
- `GATK_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
# Extract chrM reads
samtools view -b sorted.bam chrM > chrM.bam
samtools index chrM.bam

# Call variants in mitochondrial mode
gatk Mutect2 \
  -R reference.fasta \
  -I chrM.bam \
  -L chrM \
  --mitochondria-mode \
  --max-mnp-distance 0 \
  -O chrM_mutect2.vcf.gz

# Filter
gatk FilterMutectCalls \
  -R reference.fasta \
  -V chrM_mutect2.vcf.gz \
  --mitochondria-mode \
  -O chrM_mutect2_filtered.vcf.gz

# Mark possible NuMTs, given the median autosomal depth
gatk NuMTFilterTool \
  -R reference.fasta \
  -V chrM_mutect2_filtered.vcf.gz \
  --autosomal-coverage 30 \
  -O chrM_filtered.vcf.gz
```

The script runs these steps: `./scripts/20-mtoolbox.sh your_sample`.

### NuMTs

NuMTs are copies of mitochondrial DNA in the nuclear genome. Their reads can map to chrM and look like low-level heteroplasmy. GATK's NuMTFilterTool marks an allele `possible_numt` when its depth is no more than such copies could give at the sample's autosomal depth. The script reads the median autosomal depth from step 16b's mosdepth output (`mosdepth/<sample>.mosdepth.summary.txt` and `.global.dist.txt`), so run step 16b first; `AUTOSOMAL_COVERAGE=30` sets it by hand. Without either, the filter runs at depth 0 and marks nothing, and the step says so.

## Output
- `${SAMPLE}_chrM_mutect2.vcf.gz` — Raw mitochondrial variant calls
- `${SAMPLE}_chrM_mutect2_filtered.vcf.gz` — after FilterMutectCalls
- `${SAMPLE}_chrM_filtered.vcf.gz` — after FilterMutectCalls and NuMTFilterTool, with PASS or the reasons a call failed (`possible_numt` among them); the file the reports read
- Each variant includes an `AF` (allele fraction) field indicating heteroplasmy level

## Interpreting Heteroplasmy
| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — effectively fixed variant |
| 0.10-0.95 | Heteroplasmic — mixed population, clinically significant threshold varies |
| 0.03-0.10 | Low-level heteroplasmy — may be age-related somatic |
| <0.03 | Near detection limit |

Count only PASS calls: a call marked `possible_numt` or with another filter is not evidence of heteroplasmy. Calls near the two ends of the control region (about chrM:1-500 and 16,000-16,569) are less reliable: the reference is circular, so reads that span its end align poorly at both edges.

## Runtime
~15-30 minutes per sample.

## Notes
- `--mitochondria-mode` disables several filters inappropriate for mtDNA: no germline filtering, adjusted LOD thresholds, handles high copy number.
- `--max-mnp-distance 0` prevents merging nearby variants into multi-nucleotide polymorphisms.
- The GATK Docker image is large (~2.2 GB) but well-maintained and versioned.
- For disease annotation of mitochondrial variants, cross-reference with [MitoMap](https://www.mitomap.org/) or the Ensembl VEP output from step 13.
- Some mitochondrial diseases require heteroplasmy above a tissue-specific threshold (e.g., m.3243A>G MELAS requires >60% in blood).
