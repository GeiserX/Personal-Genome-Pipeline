# Step 11: Runs of Homozygosity (ROH) Analysis

## What This Does
Detects long stretches of homozygous genotypes (autozygous segments) in the genome. These arise when both copies of a chromosomal region are inherited from a common ancestor.

## Why
ROH analysis screens for consanguinity and uniparental disomy (UPD). Long ROH segments increase the risk of autosomal recessive disease by unmasking deleterious variants. ROH patterns also provide population-level ancestry information.

## Tool
- **bcftools roh** (samtools/bcftools)

## Docker Image
- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
export GENOME_DIR=/path/to/your/data
./scripts/11-roh-analysis.sh your_sample
```

What the script runs. Only PASS records (and records with no filter, as chip VCFs have) go into `bcftools roh`: DeepVariant's `RefCall` and other filtered records are not genotypes to count on. The Nextflow ROH module reads the same records with the same flags.

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  "${BCFTOOLS_IMAGE}" \
  bash -euo pipefail -c "bcftools view -f PASS,. -Ou /data/${SAMPLE}.vcf.gz \
    | bcftools roh --AF-dflt 0.4 -o /data/${SAMPLE}_roh.txt -"
```

For chip data (no `FORMAT/PL`), the script adds `-G30`.

## Output
- `${SAMPLE}/vcf/${SAMPLE}_roh.txt` — bcftools roh output: per-site states (`ST` lines) and segments (`RG` lines)
- `${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt` — the autosomal segments of **5 Mb or more** (chrom, start, end, length in bp and Mb). The script and the Nextflow module write the same file with the same threshold.

## Interpretation

### Total ROH and parental relationship

When parents are related, a child is autozygous (both copies from the same ancestor) over a fraction F of the genome, the inbreeding coefficient. The expected total of long ROH is F times the length of the autosomes, about 2,900 Mb. The table counts **autosomal segments of 5 Mb or more** (the ones the script prints), leaving out the centromeric artifacts listed below. Shorter segments come mostly from distant shared ancestry and population history, not from the parents' relationship.

| Parents' relationship | F | Expected total of ROH segments of 5 Mb or more |
|---|---|---|
| Not related | ~0 | None or a few segments |
| Second cousins | 1/64 | ~45 Mb |
| First cousins once removed | 1/32 | ~90 Mb |
| First cousins | 1/16 | ~180 Mb |
| Half siblings, uncle and niece, double first cousins | 1/8 | ~360 Mb |
| Parent and child, full siblings | 1/4 | ~720 Mb |

These are averages. Inheritance is random, so one person's total can be well above or below the value for their parents' relationship, and the ranges of neighbouring rows overlap. A total near one row is consistent with it; it does not prove it.

To compute the total from the output:

```bash
grep '^RG' "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt" \
  | awk '$3 ~ /^chr[0-9]+$/ && $6 >= 5000000 {sum += $6; n++} END {printf "%d segments, %.0f Mb\n", n, sum / 1e6}'
```

Subtract any segment that lies in one of the centromeric regions below.

### Single segments

| Individual ROH Segment | Interpretation |
|---|---|
| <1 Mb | Common, population-level background |
| 1-5 Mb | Distant shared ancestry; many of them is typical of population isolates |
| 5-10 Mb | Counted in the total above; a few can occur even when the parents are not related |
| >10 Mb | Recent shared ancestry; possible uniparental disomy if confined to one chromosome |

## Important Notes
- The script auto-detects chip data (no FORMAT/PL tag) and adds `-G30` for genotype-only mode
- `--AF-dflt 0.4` sets a default allele frequency when population AF data is unavailable — suitable for single-sample WGS
- Not done yet: a population allele-frequency file (`--AF-file`, for example from gnomAD) would replace that one default value with real frequencies per site and give better segment edges. Neither the script nor the module uses one today.
- **Known false-positive regions** (centromeric/pericentromeric, always appear as ROH in WGS):
  - chr1: 125-143 MB
  - chr9: 42-60 MB
  - chr18: 15-20 MB
- These centromeric artifacts should be excluded when calculating total ROH burden
- A single large ROH (>20 MB) confined to one chromosome may indicate uniparental disomy — verify with SNP array or parental samples
- ROH within known disease gene regions warrants checking for homozygous pathogenic variants
