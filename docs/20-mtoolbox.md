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

> **Note:** This step was originally planned for MToolBox, but no working Docker image exists for MToolBox (see [lessons-learned.md](lessons-learned.md#mtoolbox-no-working-docker-image-exists)). GATK Mutect2 is the standard clinical alternative.

## Docker Image
```
broadinstitute/gatk:4.6.1.0
```

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
  -O chrM_filtered.vcf.gz
```

## Output
- `${SAMPLE}_chrM_mutect2.vcf.gz` — Raw mitochondrial variant calls
- `${SAMPLE}_chrM_filtered.vcf.gz` — Filtered calls with PASS/FAIL status
- Each variant includes an `AF` (allele fraction) field indicating heteroplasmy level

## Interpreting Heteroplasmy
| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — effectively fixed variant |
| 0.10-0.95 | Heteroplasmic — mixed population, clinically significant threshold varies |
| 0.03-0.10 | Low-level heteroplasmy — may be age-related somatic |
| <0.03 | Near detection limit |

## Runtime
~15-30 minutes per sample.

## Notes
- `--mitochondria-mode` disables several filters inappropriate for mtDNA: no germline filtering, adjusted LOD thresholds, handles high copy number.
- `--max-mnp-distance 0` prevents merging nearby variants into multi-nucleotide polymorphisms.
- The GATK Docker image is large (~2.2 GB) but well-maintained and versioned.
- For disease annotation of mitochondrial variants, cross-reference with [MitoMap](https://www.mitomap.org/) or the Ensembl VEP output from step 13.
- Some mitochondrial diseases require heteroplasmy above a tissue-specific threshold (e.g., m.3243A>G MELAS requires >60% in blood).
