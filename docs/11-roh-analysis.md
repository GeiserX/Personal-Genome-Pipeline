# Step 11: Runs of Homozygosity (ROH) Analysis

## What This Does
Detects long stretches of homozygous genotypes (autozygous segments) in the genome. These arise when both copies of a chromosomal region are inherited from a common ancestor.

## Why
ROH analysis screens for consanguinity and uniparental disomy (UPD). Long ROH segments increase the risk of autosomal recessive disease by unmasking deleterious variants. ROH patterns also provide population-level ancestry information.

## Tool
- **bcftools roh** (samtools/bcftools)

## Docker Image
```
staphb/bcftools:1.21
```

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  staphb/bcftools:1.21 \
  bcftools roh \
    --AF-dflt 0.4 \
    /data/${SAMPLE}.vcf.gz \
    > ${GENOME_DIR}/${SAMPLE}/roh/${SAMPLE}_roh.txt

# Output: tab-delimited ROH segments with chromosomal coordinates and length
```

## Interpretation
| Total ROH | Interpretation |
|---|---|
| <100 MB | Normal outbred population |
| 100-300 MB | Possible distant consanguinity or population isolate |
| >300 MB | Suggests close consanguinity (e.g., second cousins or closer) |

| Individual ROH Segment | Interpretation |
|---|---|
| <1 MB | Common, population-level background |
| 1-10 MB | Distant shared ancestry |
| >10 MB | Recent identity-by-descent (IBD), possible UPD if single chromosome |

## Important Notes
- `--AF-dflt 0.4` sets a default allele frequency when population AF data is unavailable — suitable for single-sample WGS
- **Known false-positive regions** (centromeric/pericentromeric, always appear as ROH in WGS):
  - chr1: 125-143 MB
  - chr9: 42-60 MB
  - chr18: 15-20 MB
- These centromeric artifacts should be excluded when calculating total ROH burden
- A single large ROH (>20 MB) confined to one chromosome may indicate uniparental disomy — verify with SNP array or parental samples
- ROH within known disease gene regions warrants checking for homozygous pathogenic variants
