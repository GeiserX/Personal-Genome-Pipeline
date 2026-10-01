# Step 15: SV Quality Annotation with duphold

## What This Does
Adds depth-based quality scores to structural variant VCFs, enabling simple filtering of false positive SVs. duphold annotates each SV call with three scores derived from read-depth evidence around the breakpoints.

## Why
Manta (step 4) calls structural variants from paired-end and split-read evidence, but many calls are false positives. duphold adds depth-of-coverage annotations that allow filtering without losing true calls — it removes 63% of false positive deletions while retaining 99% of true ones.

## Tool
- **duphold** (Brent Pedersen)

## Docker Image
```
brentp/duphold:v0.2.3
```

## Annotations Added
duphold writes these as **FORMAT** fields (one value per sample), not INFO fields.

| Tag | Meaning | Interpretation |
|---|---|---|
| DHFC | Fold-change of depth inside the SV vs the rest of the chromosome it is on | General quality indicator |
| DHBFC | Fold-change of depth inside the SV vs genome bins with similar GC content | > 1.3 for duplications = true duplication (depth rises) |
| DHFFC | Fold-change of depth inside the SV vs its flanking regions | < 0.7 for deletions = true deletion (depth drops as expected) |

## Command
```bash
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  brentp/duphold:v0.2.3 \
  duphold \
  -v /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
  -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
  -f /genome/reference/Homo_sapiens_assembly38.fasta \
  -o /genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf
```

## Filtering Examples
The tags are FORMAT fields, so the expressions use `FMT/<tag>[0]` (the first sample).

```bash
# Keep only high-confidence deletions (DHFFC < 0.7)
bcftools view -i 'SVTYPE="DEL" && FMT/DHFFC[0] < 0.7' ${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf

# Keep only high-confidence duplications (DHBFC > 1.3)
bcftools view -i 'SVTYPE="DUP" && FMT/DHBFC[0] > 1.3' ${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf
```

## Runtime
~20 minutes per genome.

## Notes
- Run this AFTER Manta (step 4). Zero-cost quality improvement before AnnotSV (step 5).
- Requires the original BAM and reference FASTA — it re-calculates depth around each SV.
- Output (`duphold/${SAMPLE}_sv_duphold.vcf`, uncompressed) is the same VCF with three new FORMAT fields added. All downstream tools (AnnotSV, bcftools) work unchanged.
- Consider piping the duphold output into AnnotSV instead of the raw Manta VCF for cleaner results.
