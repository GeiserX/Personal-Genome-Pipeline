# Step 20: Mitochondrial Analysis with MToolBox

## What This Does
Comprehensive mitochondrial DNA analysis: variant calling, heteroplasmy detection, haplogroup assignment, and disease-associated variant annotation — all from the WGS BAM.

## Why
Step 12 (haplogrep3) assigns mitochondrial haplogroups from chrM variants already in the VCF. MToolBox goes deeper:
- **Heteroplasmy detection**: Identifies variants present in only a fraction of mtDNA copies (important for mitochondrial diseases)
- **Disease annotation**: Cross-references MitoMap and HMTDB databases
- **Prioritization**: Ranks variants by pathogenicity score
- **Independent calling**: Uses its own variant caller optimized for the circular mitochondrial genome

## Tool
- **MToolBox** (Calabrese et al., Bioinformatics 2014)

## Docker Image
```
robertopreste/mtoolbox:latest
```

## Input Preparation
MToolBox works from FASTQ reads mapping to the mitochondrial genome. Extract them from the BAM:
```bash
# Extract chrM reads to FASTQ
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.21 \
  bash -c "
    samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chrM | \
    samtools sort -n - | \
    samtools fastq -1 /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_R1.fastq.gz \
                   -2 /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_R2.fastq.gz \
                   -s /genome/${SAMPLE}/mtoolbox/${SAMPLE}_chrM_singleton.fastq.gz -
  "
```

## Command
```bash
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v ${GENOME_DIR}/${SAMPLE}/mtoolbox:/input \
  -v ${GENOME_DIR}/${SAMPLE}/mtoolbox/output:/output \
  robertopreste/mtoolbox:latest \
  MToolBox.sh \
    -i /input/${SAMPLE}_chrM_R1.fastq.gz \
    -I /input/${SAMPLE}_chrM_R2.fastq.gz \
    -o /output \
    -m "-t 4"
```

## Output
- `prioritized_variants.txt` — Variants ranked by pathogenicity with disease annotations
- `annotation.csv` — Full annotation with heteroplasmy levels, MitoMap references, functional impact
- `mt_classification_best_results.csv` — Haplogroup assignment (independent of haplogrep3)
- `OUT.vcf` — Mitochondrial VCF with heteroplasmy fractions in FORMAT field

## Interpreting Heteroplasmy
| Level | Meaning |
|---|---|
| >90% | Homoplasmic — effectively fixed variant |
| 10-90% | Heteroplasmic — mixed population, clinically significant threshold varies by variant |
| 1-10% | Low-level heteroplasmy — may be age-related somatic mutation |
| <1% | Background noise |

## Runtime
~15-30 minutes per sample.

## Notes
- Run AFTER alignment (step 2). Independent of other analyses.
- MToolBox re-aligns reads to a mitochondrial reference (rCRS), so it may detect variants missed by the main pipeline's chrM calls.
- Heteroplasmy detection is MToolBox's main advantage over haplogrep3 — clinically relevant for mitochondrial disease assessment.
- The `--user root` flag is needed for write access to the output directory.
- Some mitochondrial variants are only pathogenic above certain heteroplasmy thresholds (e.g., m.3243A>G MELAS requires >60%).
