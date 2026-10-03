# Step 6: ClinVar Pathogenic Variant Screening

## What This Does
Intersects your sample VCF against the ClinVar database of known pathogenic variants, identifying any positions where your genome carries a clinically reported disease variant.

## Why
ClinVar is the most widely used public database of clinically reported variants. This screen catches pathogenic SNPs and indels that have been submitted by clinical labs — carrier status, dominant disease risk, and pharmacogenomic flags. Note that ClinVar entries vary in evidence quality (see star ratings in [interpreting-results.md](interpreting-results.md)).

## Tool
- **bcftools norm**, **isec** and **annotate** — split and left-align both files, keep the sample's records whose allele is in ClinVar, and copy ClinVar's gene, significance and review status onto them

## Docker Image
- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sample VCF from DeepVariant (step 3)
- Reference FASTA and its `.fai`: `reference/Homo_sapiens_assembly38.fasta` (used to left-align indels)
- `clinvar_pathogenic_chr.vcf.gz` from reference setup (step 00) — chr-prefixed, filtered to Pathogenic/Likely_pathogenic only

## Command
```bash
export GENOME_DIR=/path/to/data
./scripts/06-clinvar-screen.sh <sample_name>

# For long-read Clair3 output:
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh <sample_name>
```

### What the Script Does

1. Lists the contigs each file holds records on (`bcftools index -s`) and keeps those the reference has too: ClinVar's contigs that are in the reference, and the sample's contigs among those. `bcftools norm` stops at the first record whose contig the reference lacks, such as a provider's unplaced scaffold or one of the `NT_` contigs the full ClinVar file has, and a record on a contig the other file lacks cannot match anyway. The step prints how many of the sample's records it leaves out and on which contigs. When no contig is shared it stops with an error instead of reporting zero hits, and says which file is named the other way (`1, 2, ...` against `chr1, chr2, ...`).
2. Builds `clinvar/clinvar_pathogenic_chr.norm.vcf.gz` beside the ClinVar file the first time, and again whenever the source file is newer: ClinVar's records on the reference's contigs, multiallelic records split (`bcftools norm -m -any`) and indels left-aligned against the reference. Every sample then reuses it. The Nextflow module reuses this file when it sits beside the `--clinvar` file and is newer; otherwise it normalises ClinVar inside the task.
3. Filters the sample VCF to PASS records. If the VCF has no PASS record at all (callers that leave FILTER as `.`), it uses `-f .,PASS` instead and prints a notice saying so. If no record is left, the step stops with an error.
4. Splits and left-aligns the sample the same way. `bcftools isec` matches only identical REF/ALT, so a pathogenic allele inside a multiallelic record (genotype `1/2`) is found only after this split.
5. Keeps the sample's records whose allele is in ClinVar (`bcftools isec -n=2 -w1`), copies ClinVar's `ID`, `GENEINFO`, `CLNSIG` and `CLNREVSTAT` onto them (`bcftools annotate --pair-logic exact`), and keeps only those whose genotype carries an ALT allele (`bcftools view -i 'GT~"[1-9]"'`, a non-zero allele index anywhere in the genotype). The match is by allele, so without this a `0/0` or `./.` record, the `0/0` half of a split multiallelic record and a reference-only row with ALT `.` would be listed as hits. A half call such as `./1` carries the ALT allele and is a hit; `GT="alt"` would drop it. The filter runs after the split, where it sees each allele's own genotype. No command line goes into the hits file's header (`--no-version`), so ClinVar's path on your machine stays out of it.
6. Prints the hits grouped by ClinVar's review stars (`bin/clinvar_hits.awk`): 4 practice guideline, 3 expert panel, 2 multiple submitters with no conflict, 1 a single submitter or conflicting classifications, 0 no assertion criteria. The screening file keeps every Pathogenic/Likely_pathogenic submission whatever its review status, so a zero-star hit counts as a hit; the stars say how much weight it deserves. The Nextflow module prints the same counts.

## Output

| File | Description |
|---|---|
| `clinvar/${SAMPLE}_clinvar_hits.vcf` | **The hits: your records that carry a ClinVar Pathogenic/Likely_pathogenic allele, with ClinVar's ID, gene (`GENEINFO`), significance (`CLNSIG`) and review status (`CLNREVSTAT`). Your genotype is in the sample column.** |
| `clinvar/${SAMPLE}_clinvar_hits.tsv` | The same hits, one row each: `chrom`, `pos`, `ref`, `alt`, `genotype`, `clinvar_id`, `geneinfo`, `clnsig`, `clnrevstat` |
| `clinvar/${SAMPLE}_pass.vcf.gz` | Filtered, split and left-aligned sample VCF (intermediate). It holds the records on the contigs shared with ClinVar and the reference only; versions before the shared-contig restriction kept every contig |
| `clinvar/clinvar_pathogenic_chr.norm.vcf.gz` (in `${GENOME_DIR}`) | Normalised ClinVar, shared by every sample |

Both reports (step 24 and `generate-report.sh`) read the hits file and show gene, genotype (het/hom), significance, review status and stars for each hit, the best-reviewed first, with the ClinVar file's release date. Older versions of this step wrote `clinvar/isec/0002.vcf`; that file holds only the sample's side of the intersection, with no gene or significance, and nothing reads it any more. Rerun step 6 to get the hits file.

`${SAMPLE}_clinvar_hits.tsv` is the table to read or load elsewhere. To list the hits from the VCF yourself:

```bash
source versions.env   # from the repository root
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools query -f '%CHROM:%POS %REF>%ALT [%GT] %INFO/GENEINFO %INFO/CLNSIG %INFO/CLNREVSTAT\n' \
  /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf
```

## Interpreting Results

This step screens against **Pathogenic and Likely_pathogenic variants only** — benign and VUS entries are excluded at the database level (see step 00 reference setup). Every hit is an allele that ClinVar classifies as Pathogenic or Likely_pathogenic, and the hit shows your genotype for it.

| Scenario | Meaning | Action |
|---|---|---|
| Heterozygous + autosomal recessive | Healthy carrier | Note for family planning only |
| Homozygous + autosomal recessive | Possibly affected — requires clinical confirmation | Investigate — confirm with clinical evaluation and ClinVar review status |
| Any genotype + autosomal dominant | Possibly affected — requires clinical confirmation | Investigate — check penetrance, ClinVar review status, and phenotype |
| Compound het (two variants, same gene) | Potentially affected (recessive) | Check if variants are on different alleles (phasing) |

## Limitations

- The reports show ClinVar's review status (`CLNREVSTAT`), but a status such as "criteria provided, single submitter" still rests on one lab. Always check the full ClinVar entry before acting on any result.
- Matching is by allele after both files are split and left-aligned, in the bash script and in the Nextflow module alike. A complex variant that the caller writes differently from ClinVar (for example an MNP against two SNVs) can still be missed.
- The screen sees only small variants in the VCF. Copy-number losses and gene deletions (for example SMN1 in spinal muscular atrophy) are invisible to it.
- Results are **research-grade**, not clinical diagnoses. Do not make medical decisions based solely on this output.

## Important Notes
- Most hits will be **heterozygous carriers of recessive conditions** — this is normal and expected
- A typical 30X WGS shows 0-10 pathogenic/likely pathogenic overlaps. The majority are benign carrier states for recessive conditions
- Focus review on: homozygous pathogenic, any autosomal dominant pathogenic, and compound heterozygous variants in the same gene
- The ClinVar pathogenic database must be chr-prefixed to match the BAM/VCF coordinate system (done in step 00)
- ClinVar is updated monthly — re-download periodically to catch newly classified variants
