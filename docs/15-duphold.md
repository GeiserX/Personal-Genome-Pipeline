# Step 15: SV Quality Annotation with duphold

## What This Does
Adds depth-based quality scores to structural variant VCFs, enabling simple filtering of false positive SVs. duphold annotates each SV call with three scores derived from read-depth evidence around the breakpoints.

## Why
Manta (step 4) calls structural variants from paired-end and split-read evidence, but many calls are false positives. duphold adds depth-of-coverage annotations that allow filtering without losing true calls — it removes 63% of false positive deletions while retaining 99% of true ones.

## Tool
- **duphold** (Brent Pedersen)

## Docker Image
- `DUPHOLD_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Annotations Added
duphold writes these as **FORMAT** fields (one value per sample), not INFO fields.

| Tag | Meaning | Interpretation |
|---|---|---|
| DHFC | Fold-change of depth inside the SV vs the rest of the chromosome it is on | General quality indicator |
| DHBFC | Fold-change of depth inside the SV vs genome bins with similar GC content | > 1.3 for duplications supports a real duplication (depth rises) |
| DHFFC | Fold-change of depth inside the SV vs its flanking regions | < 0.7 for deletions supports a real deletion (depth drops as expected) |

## Command
```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
docker run --rm \
  --cpus 4 --memory 4g \
  -v ${GENOME_DIR}:/genome \
  "${DUPHOLD_IMAGE}" \
  duphold \
  -v /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
  -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
  -f "/genome/${REF_FASTA}" \
  -o /genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf
```

## The Depth Filter
The script then writes `duphold/${SAMPLE}_sv_filtered.vcf.gz` with its index: the same calls without the deletions and duplications the depth does not support. It is the filter of the Nextflow `DUPHOLD_FILTER`, and step 5 (AnnotSV) annotates this file. The tags are FORMAT fields, so the expression uses `FMT/<tag>[0]` (the first sample):

```bash
# Drop deletions with DHFFC >= 0.7 and duplications with DHBFC <= 1.3;
# other SV types, and records without a value, stay.
docker run --rm -v ${GENOME_DIR}:/genome "${BCFTOOLS_IMAGE}" \
  bcftools view \
  -e '(INFO/SVTYPE="DEL" && FMT/DHFFC[0] >= 0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0] <= 1.3)' \
  -Oz -o /genome/${SAMPLE}/duphold/${SAMPLE}_sv_filtered.vcf.gz \
  /genome/${SAMPLE}/duphold/${SAMPLE}_sv_duphold.vcf
docker run --rm -v ${GENOME_DIR}:/genome "${BCFTOOLS_IMAGE}" \
  bcftools index -t /genome/${SAMPLE}/duphold/${SAMPLE}_sv_filtered.vcf.gz
```

The script prints how many records it kept. How many SVs the filter removes from a real genome has not been measured here.

## Runtime
~20 minutes per genome.

## Notes
- Run this AFTER Manta (step 4) and before AnnotSV (step 5), which annotates the filtered file when it is there and Manta's calls otherwise.
- Requires the original BAM and reference FASTA — it re-calculates depth around each SV.
- Outputs: `duphold/${SAMPLE}_sv_duphold.vcf` (uncompressed), the same VCF with three new FORMAT fields, and `duphold/${SAMPLE}_sv_filtered.vcf.gz` with its `.tbi`, after the depth filter above.
