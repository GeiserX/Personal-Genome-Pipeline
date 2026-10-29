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

When parents are related, a child is autozygous (both copies from the same ancestor) over a fraction F of the genome, the inbreeding coefficient. The expected total of long ROH is F times the length of the autosomes, about 2,900 Mb. The table counts **autosomal segments of 5 Mb or more** (the ones the script prints), without their part in the centromere and heterochromatin bands above. Those bands cover about 217 Mb of chr1 to chr22 (7.5% of 2,875 Mb), and the table below uses the full autosomes, so expect a reported total up to about 7.5% under its row. Shorter segments come mostly from distant shared ancestry and population history, not from the parents' relationship.

| Parents' relationship | F | Expected total of ROH segments of 5 Mb or more |
|---|---|---|
| Not related | ~0 | None or a few segments |
| Second cousins | 1/64 | ~45 Mb |
| First cousins once removed | 1/32 | ~90 Mb |
| First cousins | 1/16 | ~180 Mb |
| Half siblings, uncle and niece, double first cousins | 1/8 | ~360 Mb |
| Parent and child, full siblings | 1/4 | ~720 Mb |

These are averages. Inheritance is random, so one person's total can be well above or below the value for their parents' relationship, and the ranges of neighbouring rows overlap. A total near one row is consistent with it; it does not prove it.

Both reports print this as **F_ROH**, the measure of McQuillan et al. 2008 (PMID 18760389): the length of the autosomal segments of 5 Mb or more, without their band parts, divided by the length of chr1 to chr22 without the bands (2,875,001,522 bp less 217,000,000 bp, so 2,658,001,522 bp). Beside it they name the row of the table above whose F is nearest to F_ROH on a log scale: the boundaries are the geometric means of neighbouring rows (for example 1/16 and 1/8 meet at about 0.088), and an F_ROH under 1/128 (about 21 Mb) reads "no sign of related parents". When one chromosome carries most of the runs (next section), there is no band: that chromosome fills F_ROH by itself.

To compute the total the table means from the summary file, which takes the band parts out before it keeps the segments of 5 Mb or more, as the reports do:

```bash
awk -F'\t' 'NR > 3 {sum += $4; n++} END {printf "%d segments, %.0f Mb\n", n, sum / 1e6}' \
  "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh_summary.txt"
```

### One chromosome carrying most of the runs: possible uniparental disomy

In uniparental disomy (UPD) both copies of a chromosome come from one parent. When they are two copies of the same parental chromosome (isodisomy), that chromosome reads as one long run of homozygosity, or a few, while the rest of the genome has few long runs. Nakka et al. 2019 (PMID 31607426) found 675 cases of UPD this way, from runs of homozygosity and segments shared with relatives, among 4,400,363 23andMe research participants and 431,094 UK Biobank participants.

The reports flag a chromosome when its segments of 5 Mb or more (band parts left out) add up to at least 20 Mb **and** to more than half of the length of all such autosomal segments. Both numbers are our thresholds, not a published cut-off: 20 Mb on one chromosome is the long-standing sign of isodisomy (see the uniparental disomy note under [Important Notes](#important-notes)), and the share rule keeps out a child of related parents, whose runs spread over many chromosomes. The flag says "may": a run that long can also be chance in a child of related parents.

A flagged chromosome is looked up in `assets/imprinted_regions.tsv`, the chromosomes where UPD from one parent causes an imprinting disorder (overview: Eggermann et al. 2023, Nat Rev Dis Primers, PMID 37386011). A row with a region applies only when a run on that chromosome overlaps it; the others apply to the whole chromosome.

| Chromosome or region | UPD from the mother | UPD from the father | Source (PMID) |
|---|---|---|---|
| chr6 (6q24, PLAGL1, HYMAI) | none known | Transient neonatal diabetes mellitus | 20301706 |
| chr7 (GRB10 at 7p12, MEST at 7q32) | Silver-Russell syndrome | none known | 27585961 |
| 11p15, chr11:0-22,000,000 (H19/IGF2, KCNQ1OT1) | Silver-Russell syndrome | Beckwith-Wiedemann syndrome | 29377879, 27585961 |
| chr14 (14q32.2, DLK1-MEG3) | Temple syndrome | Kagami-Ogata syndrome | 24891339, 26377239 |
| 15q11-q13, chr15:20,500,000-33,400,000 (SNRPN, UBE3A) | Prader-Willi syndrome | Angelman syndrome | 22237428, 27615419 |
| chr20 (20q13.32, GNAS) | Mulchandani-Bhoj-Conlin syndrome (growth failure) | Pseudohypoparathyroidism type 1B | 26248010, 29959430 |

`bin/collect_summary.py` holds the same rows (a unit test checks that the two agree), because the Nextflow report task reads `bin/` but not `assets/`. A flagged chromosome with no row still matters: a recessive variant inside the run, if there is one, would be homozygous.

What this cannot tell, and the reports say so:

- **Which parent.** A single genome has no parent of origin; both columns are listed. The parents' genotypes (step 38, the family workflow) confirm a UPD and tell which. Where a row applies, a methylation test of its imprinted region tests for that imprinting disorder; on its own it does not prove a UPD, since an imprinting defect gives the same result.
- **Heterodisomy.** Two different copies from one parent leave no run at all, so it is not detected here. It needs the parents' genotypes.
- **Mosaic or segmental UPD** that covers only a short stretch, or a share of the cells (common in Beckwith-Wiedemann syndrome), may stay under the thresholds.

### Single segments

| Individual ROH Segment | Interpretation |
|---|---|
| <1 Mb | Common, population-level background |
| 1-5 Mb | Distant shared ancestry; many of them is typical of population isolates |
| 5-10 Mb | Counted in the total above; a few can occur even when the parents are not related |
| >10 Mb | Recent shared ancestry; possible uniparental disomy if confined to one chromosome (the reports flag one chromosome from 20 Mb, see above) |

## Important Notes
- The script auto-detects chip data (no FORMAT/PL tag) and adds `-G30` for genotype-only mode
- `--AF-dflt 0.4` sets a default allele frequency when population AF data is unavailable — suitable for single-sample WGS
- bcftools roh reads a site's frequency from INFO/AC and INFO/AN when the record has them, and uses 0.4 only where it does not. DeepVariant writes neither, so on its calls every site uses 0.4. A VCF from another caller may carry them: in a single-sample file they are the sample's own allele counts, not a population's.
- No population allele frequencies are used, and the reports say so: read the total and the short segments as rough, the long segments are the signal. A population allele-frequency file (`--AF-file`, for example from gnomAD) would replace the one default value with real frequencies per site and give better segment edges; neither the script nor the module uses one today, and its effect on the total has not been measured.
- The centromere and heterochromatin bands (above) are taken out of every count; the step and the reports do it, so there is nothing to subtract by hand
- Long runs concentrated on one chromosome may indicate uniparental disomy; the reports flag a chromosome whose segments of 5 Mb or more add up to at least 20 MB and to more than half of the length of all such autosomal segments (above). Confirm with the parents' genotypes; a methylation test of a listed imprinted region tests for its imprinting disorder
- ROH within known disease gene regions warrants checking for homozygous pathogenic variants
