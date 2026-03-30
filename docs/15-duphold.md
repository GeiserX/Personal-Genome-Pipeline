# Step 15: SV Quality Annotation with duphold

## What This Does
Adds depth-based quality scores to structural variant VCFs, enabling simple filtering of false positive SVs. duphold annotates each SV call with three scores derived from read-depth evidence around the breakpoints.

## Why
Manta (step 4) calls structural variants from paired-end and split-read evidence, but many calls are false positives. duphold adds depth-of-coverage annotations that allow filtering without losing true calls — it removes 63% of false positive deletions while retaining 99% of true ones.

## Tool
- **duphold** (Brent Pedersen)

## Docker Image
```
brentp/duphold:latest
```

## Annotations Added
| Tag | Meaning | Interpretation |
|---|---|---|
| DHFFC | Fold-change of depth inside SV vs flanking regions | < 0.7 for deletions = true deletion (depth drops as expected) |
| DHBFC | Fold-change of depth inside SV vs background chromosome depth | > 1.3 for duplications = true duplication (depth rises) |
| DHFC | Fold-change combining both flanking and background | General quality indicator |

## Command
```bash
docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  brentp/duphold:latest \
  duphold \
  -v /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
  -b /genome/${SAMPLE}/sorted.bam \
  -f /genome/reference/Homo_sapiens_assembly38.fasta \
  -o /genome/${SAMPLE}/manta/results/variants/diploidSV.duphold.vcf.gz
```

## Filtering Examples
```bash
# Keep only high-confidence deletions (DHFFC < 0.7)
bcftools view -i 'SVTYPE="DEL" && DHFFC < 0.7' diploidSV.duphold.vcf.gz

# Keep only high-confidence duplications (DHBFC > 1.3)
bcftools view -i 'SVTYPE="DUP" && DHBFC > 1.3' diploidSV.duphold.vcf.gz
```

## Runtime
~20 minutes per genome.

## Notes
- Run this AFTER Manta (step 4). Zero-cost quality improvement before AnnotSV (step 5).
- Requires the original BAM and reference FASTA — it re-calculates depth around each SV.
- Output is the same VCF with three new INFO fields added. All downstream tools (AnnotSV, bcftools) work unchanged.
- Consider piping the duphold output into AnnotSV instead of the raw Manta VCF for cleaner results.
