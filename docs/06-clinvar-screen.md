# Step 6: ClinVar Pathogenic Variant Screening

## What This Does
Intersects your sample VCF against the ClinVar database of known pathogenic variants, identifying any positions where your genome carries a clinically reported disease variant.

## Why
ClinVar is the most widely used public database of clinically reported variants. This screen catches pathogenic SNPs and indels that have been submitted by clinical labs — carrier status, dominant disease risk, and pharmacogenomic flags. Note that ClinVar entries vary in evidence quality (see star ratings in [interpreting-results.md](interpreting-results.md)).

## Tool
- **bcftools norm**, **isec** and **annotate** — split and left-align both files, keep the sample's records whose allele is in ClinVar, and copy ClinVar's gene, significance and review status onto them
- **bin/clinvar_join_adjacent.awk** with **bcftools +fill-from-fasta** — build the combined allele of two neighbouring records, so a delins the caller wrote as two records is matched as one allele too

## Docker Image
- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites
- Sample VCF from DeepVariant (step 3)
- Reference FASTA and its `.fai`: `reference/GRCh38_no_alt_analysis_set.fasta` (used to left-align indels and to fill in the bases between the two parts of a combined allele)
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
5. Builds the combined alleles of neighbouring records. A caller can write one complex allele, a delins, as two records a few bases apart on the same haplotype, for example a deletion and an insertion. Matched one by one, the parts can hit two single-submitter ClinVar records and miss the better reviewed record ClinVar holds for the whole allele. `bcftools norm -m +any` joins only records at the same position, so it cannot build that allele. `bin/clinvar_join_adjacent.awk` pairs two records when the second starts at most 10 bases after the first ends (records that overlap are never paired) and their genotypes put both ALT alleles on one haplotype: phased with the same genotype in the same phase set (`FORMAT/PS`), or unphased with the same genotype (`0/1` and `1/0` count as the same het call; a call with a missing allele such as `./1` is never paired). The combined record runs from the first part's start to the second part's end; `bcftools +fill-from-fasta` takes the bases between the parts from the reference, and `bcftools norm` trims and left-aligns it the way ClinVar was. The gene is not checked here, because the sample VCF has no gene annotation yet: a combined allele only counts when ClinVar holds a record for it. The step prints how many combined alleles it built. Only pairs are built: an allele split into three or more records is matched by its pairs only.
6. Keeps the sample's records and combined alleles whose allele is in ClinVar (`bcftools isec -n=2 -w1`), copies ClinVar's `ID`, `GENEINFO`, `CLNSIG` and `CLNREVSTAT` onto them (`bcftools annotate --pair-logic exact`), and keeps only those whose genotype carries an ALT allele (`bcftools view -i 'GT~"[1-9]"'`, a non-zero allele index anywhere in the genotype). The match is by allele, so without this a `0/0` or `./.` record, the `0/0` half of a split multiallelic record and a reference-only row with ALT `.` would be listed as hits. A half call such as `./1` carries the ALT allele and is a hit; `GT="alt"` would drop it. The filter runs after the split, where it sees each allele's own genotype. No command line goes into the hits file's header (`--no-version`), so ClinVar's path on your machine stays out of it.
7. Marks the hits that belong together. A matched combined allele is a hit of its own and carries `INFO/PGP_JOINED_FROM`, its two parts as `CHROM:POS:REF:ALT`. Each of its parts that is a hit too keeps its row and gets `INFO/PGP_PART_OF`, the combined record it belongs to, so the reader sees one allele with the combined record's own review status, not two unrelated hits. Two hits in one gene that could be one allele by the rule of item 5 but matched no combined record each get `INFO/PGP_ONE_ALLELE`, the other hit: possibly one allele, check ClinVar for a combined record.
8. Prints the hits grouped by ClinVar's review stars (`bin/clinvar_hits.awk`): 4 practice guideline, 3 expert panel, 2 multiple submitters with no conflict, 1 a single submitter or conflicting classifications, 0 no assertion criteria, each with its combined-allele note when it has one. The screening file keeps every Pathogenic/Likely_pathogenic submission whatever its review status, so a zero-star hit counts as a hit; the stars say how much weight it deserves. The Nextflow module runs the same steps and prints the same counts.

## Output

| File | Description |
|---|---|
| `clinvar/${SAMPLE}_clinvar_hits.vcf` | **The hits: your records that carry a ClinVar Pathogenic/Likely_pathogenic allele, with ClinVar's ID, gene (`GENEINFO`), significance (`CLNSIG`) and review status (`CLNREVSTAT`). Your genotype is in the sample column.** A combined allele of two neighbouring records is a record of its own, with `PGP_JOINED_FROM`; its parts carry `PGP_PART_OF`, and two hits that may be one allele carry `PGP_ONE_ALLELE` |
| `clinvar/${SAMPLE}_clinvar_hits.tsv` | The same hits, one row each: `chrom`, `pos`, `ref`, `alt`, `genotype`, `clinvar_id`, `geneinfo`, `clnsig`, `clnrevstat` |
| `clinvar/${SAMPLE}_pass.vcf.gz` | Filtered, split and left-aligned sample VCF (intermediate). It holds the records on the contigs shared with ClinVar and the reference only; versions before the shared-contig restriction kept every contig |
| `clinvar/clinvar_pathogenic_chr.norm.vcf.gz` (in `${GENOME_DIR}`) | Normalised ClinVar, shared by every sample |

Both reports (step 24 and `generate-report.sh`) read the hits file and show gene, genotype (het/hom), significance, review status, stars and the combined-allele note ("combined allele of ...", "part of a combined allele: ...", or "possibly one allele with ...: check ClinVar for a combined record") for each hit, the best-reviewed first, with the ClinVar file's release date. Older versions of this step wrote `clinvar/isec/0002.vcf`; that file holds only the sample's side of the intersection, with no gene or significance, and nothing reads it any more. Rerun step 6 to get the hits file.

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
- Matching is by allele after both files are split and left-aligned, in the bash script and in the Nextflow module alike, plus the combined allele of two neighbouring records (item 5). A complex variant can still be missed when the caller splits it into three or more records, when its parts are more than 10 bases apart, when their genotypes differ or they sit on different phase sets, or when ClinVar's record normalises to another representation. Two hits in one gene within 10 bases whose genotypes allow one allele, and that matched no combined record, are flagged: check ClinVar for a combined record by hand.
- The screen sees only small variants in the VCF. Copy-number losses and gene deletions (for example SMN1 in spinal muscular atrophy) are invisible to it.
- Population frequency is not checked. A common allele with a P/LP record is often low-penetrance or misclassified, so look up a homozygous P/LP hit in gnomAD before reading much into it. The report's hit-count badge is never green for this reason.
- Results are **research-grade**, not clinical diagnoses. Do not make medical decisions based solely on this output.

## Important Notes
- Most hits will be **heterozygous carriers of recessive conditions** — this is normal and expected
- A typical 30X WGS shows 0-10 pathogenic/likely pathogenic overlaps. The majority are benign carrier states for recessive conditions
- Focus review on: homozygous pathogenic, any autosomal dominant pathogenic, and compound heterozygous variants in the same gene
- The ClinVar pathogenic database must be chr-prefixed to match the BAM/VCF coordinate system (done in step 00)
- ClinVar is updated monthly — re-download periodically to catch newly classified variants
