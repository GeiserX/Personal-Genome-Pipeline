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
- `${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt` — the autosomal segments with **5 Mb or more outside the centromere and heterochromatin bands** (chrom, start, end, and `length_bp` and `length_mb` of the part outside the bands). The script and the Nextflow module write the same file with the same threshold and bands.

### Centromere and heterochromatin bands

Few markers map in GRCh38's centromeres (cytoband stain `acen`), its large heterochromatin blocks (`gvar`, such as 1q12, 9q12 and 16q11.2) and the acrocentric stalks (`stalk`), so `bcftools roh` reads long homozygous runs there that are artefacts. The script, the Nextflow module and the reports take the part of each segment that lies in one of these bands out of its length; a segment wholly inside one counts as nothing. The bands are the `acen`, `gvar` and `stalk` rows of UCSC's `cytoBand.hg38.txt` (the file `setup.sh` installs for step 10), neighbouring rows merged into one span around each centromere (two on chrY), written into `scripts/11-roh-analysis.sh`, `modules/local/roh/main.nf` and `bin/collect_summary.py` as the same line (a unit test checks they agree). On chr1 the span is 121.7-143.2 Mb, on chr9 42.2-61.5 Mb, on chr18 15.4-21.5 Mb.

The reports count autosomes only: Segments, Total and Largest keep chr1 to chr22, so they leave out chrX and chrY (a man's single X reads as one long run) and any other contig (chrM, chrEBV, chrUn_*), and the part of each segment in a band. They say how many chrX and chrY segments and how many megabases in bands were left out.

## Interpretation

### Total ROH and parental relationship

When parents are related, a child is autozygous (both copies from the same ancestor) over a fraction F of the genome, the inbreeding coefficient. The expected total of long ROH is F times the length of the autosomes, about 2,900 Mb. The table counts **autosomal segments of 5 Mb or more** (the ones the script prints), without their part in the centromere and heterochromatin bands above. Shorter segments come mostly from distant shared ancestry and population history, not from the parents' relationship.

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

That sums the raw segments; the summary file already has the band parts taken out, so `awk -F'\t' 'NR > 3 {sum += $4; n++} END {printf "%d segments, %.0f Mb\n", n, sum / 1e6}' "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"` gives the total the table means.

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
- No population allele frequencies are used, and the reports say so: read the total and the short segments as rough, the long segments are the signal. A population allele-frequency file (`--AF-file`, for example from gnomAD) would replace the one default value with real frequencies per site and give better segment edges; neither the script nor the module uses one today, and its effect on the total has not been measured.
- The centromere and heterochromatin bands (above) are taken out of every count; the step and the reports do it, so there is nothing to subtract by hand
- A single large ROH (>20 MB) confined to one chromosome may indicate uniparental disomy — verify with SNP array or parental samples
- ROH within known disease gene regions warrants checking for homozygous pathogenic variants
