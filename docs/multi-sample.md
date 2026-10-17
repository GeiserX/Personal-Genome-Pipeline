# Multi-Sample Comparison Guide

Ran the pipeline on two or more people? This guide explains how to compare their results for carrier screening, family planning, and inherited disease investigation.

> **Privacy between the samples of one data directory.** The container of every step script gets the whole `GENOME_DIR` read-only at `/genome` (the references, the databases, and every sample's folder), and only its own sample's folder writable. So a container is not a confidentiality boundary between the samples of one data directory: any tool run for one person can read the other people's files. To keep two people's data apart, give each their own `GENOME_DIR`.

---

## Family analysis (step 38, opt-in)

When two or more related people are processed, step 38 analyses them together: one joint call, a check that the samples are the people the pedigree says, and for a child with a sequenced parent, which allele came from which parent. It runs only when you give it a pedigree; a single-sample run is not changed.

### The pedigree

A tab-separated file with a header line and one row per person ([`assets/pedigree.example.tsv`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/assets/pedigree.example.tsv)):

```text
sample	sex	father	mother	partner
child_a	female	father_a	mother_a	.
father_a	male	.	.	mother_a
mother_a	female	.	.	father_a
```

`sample` is the sample name the pipeline uses. `sex` is `male`, `female` or `unknown`. `father`, `mother` and `partner` name another person, or `.` (or empty) for none. A parent or partner named without a row of their own is a person who was not sequenced: a child's unsequenced father can still be named. The step stops, before any work, when a column name is unknown, a person has two rows, a father is female or a mother male, someone is their own parent or ancestor, two people name different partners, the samplesheet's sex contradicts the pedigree's, a gVCF holds another sample name than its own, or fewer than two named people have a gVCF.

### Running it

From the scripts, once each member has a gVCF from step 03 (`<sample>/vcf/<sample>.g.vcf.gz`) and somalier's sites are installed (`setup.sh --sample-qc-data`):

```bash
GENOME_DIR=/path/to/data ./scripts/38-family.sh my_family pedigree.tsv
```

The results go to `$GENOME_DIR/my_family/`, a folder of its own: the family name follows the sample-name rule and may not be a sample's name. Without the second argument the step reads `$GENOME_DIR/my_family/pedigree.tsv`.

In the Nextflow pipeline, give each member a samplesheet row and add `--pedigree pedigree.tsv --somalier_sites <sites.hg38.vcf.gz>`. A person processed earlier is a VCF row with the `gvcf` and `gvcf_index` columns; a row from FASTQ or a BAM is called in the same run and its gVCF is used. The family results go to `<outdir>/family/`.

### What it does

1. **Joint call.** GLnexus (`GLNEXUS_IMAGE`, its DeepVariantWGS preset) merges the members' DeepVariant gVCFs into one VCF, `family.vcf.gz`. A site where one person has a variant gets the others' genotypes from their reference blocks, so "not called" and "matches the reference" are told apart.
2. **Relationship check.** somalier extracts each member's genotypes at its sites from the gVCF and `somalier relate` measures the relatedness of every pair: about the share of the genome two people inherited from common ancestors (1 the same person, 0.5 parent and child or full siblings, 0.25 second degree, 0 unrelated). The pedigree gives the expected value. A pair must fall in the same degree as declared (bins at 0.707, 0.354, 0.177 and 0.088); a declared parent and child must also be opposite homozygotes at under 1% of the sites, which siblings are not; a pair declared unrelated or distant must not measure first degree. A pair that disagrees stops the step with exit 1, before any phasing, and `relationships.tsv` says which pair and why: two samples were swapped, or a pedigree line is wrong. A pair compared at fewer than 1,000 sites is reported as not checked.
3. **Phasing by inheritance.** For each child whose two parents are sequenced, WhatsHap (`WHATSHAP_IMAGE`, `whatshap phase --ped`) phases the family by inheritance alone; no reads are given. A site where the child and both parents are heterozygous stays unphased, because no genotype tells it.
4. **The report.** For each child with a sequenced parent, `<child>_inheritance.tsv.gz` gives, at every site of the joint VCF where the child and a sequenced parent have a genotype, the allele from the father and the allele from the mother when the genotypes tell (`how` is `genotype`), or `ambiguous` (both parents heterozygous like the child) or `inconsistent` (a genotype error, a deletion or a new mutation). On chrX and chrY outside the pseudoautosomal regions a son's X comes from his mother and his Y from his father; the child's sex must be known there. chrM is left out. With one parent sequenced, the allele the other parent passed on is known wherever the sequenced parent is homozygous or differs from the child: that is the half of the unsequenced parent's genome the child carries, and `family_report.md` gives its share of the sites. The half the child did not inherit is not seen.

   De novo candidates (`de_novo.tsv`, a child allele neither parent has) are called only when both parents are sequenced. A candidate is `pass` with GQ 20 or more and depth 10 or more in all three and, when the child is heterozygous, 25% to 75% of the child's reads on the new allele; otherwise `low`. Most candidates of a whole genome are genotype errors: confirm one that matters by another method.

   When step 23 (the clinical filter) ran for a child, every pair of its candidate variants the child carries heterozygous in one gene is listed in `compound_het.tsv`: `in trans` when one came from each parent (both copies of the gene are hit), `in cis` when both came from one parent, `unknown` when the genotypes do not tell.

Two members marked as partners also get the couple carrier screen of [step 39](39-couple.md): `scripts/39-couple.sh`, or in the Nextflow pipeline the `COUPLE_CARRIER` process of the family workflow (with `clinvar` in `--tools`).

### Outputs

The Nextflow pipeline writes them to `<outdir>/family/`, step 38 to `$GENOME_DIR/<family_name>/` (with the `pedigree.tsv` it read).

| File | Content |
|---|---|
| `family.vcf.gz` (+ `.tbi`) | the joint call |
| `members.tsv`, `family.ped`, `trios.ped` | the pedigree as read, somalier's pedigree, WhatsHap's trios |
| `somalier/family.pairs.tsv`, `family.samples.tsv`, `family.html` | somalier relate's own tables and plot |
| `relationships.tsv` | each pair: declared relationship, expected and measured relatedness, opposite homozygotes, sites, verdict and reason |
| `family.phased.vcf.gz` (+ `.tbi`) | WhatsHap's phasing (only when a child has both parents sequenced) |
| `<child>_inheritance.tsv.gz` | allele from each parent, site by site |
| `de_novo.tsv`, `compound_het.tsv` | as above |
| `family_report.md` | all of it in words, pair by pair and child by child |

### Limits

- The gVCFs must be DeepVariant's: the GLnexus preset reads DeepVariant's fields.
- No reads are used for phasing, so an all-heterozygous trio site stays unphased, and a duo is phased only where the genotypes tell.
- The relatedness is one number per pair for the whole genome. Shared segments (where two people share a stretch of DNA) are not reported.
- Runtime on whole genomes has not been measured yet; the CI test runs it on 2.45 Mb of the GIAB Ashkenazi trio.

---

## Carrier Screening for Partners

The most immediately useful multi-sample analysis: checking whether both partners carry pathogenic variants in the **same recessive gene**. If so, each child has a 25% chance of being affected.

[Step 39](39-couple.md) does this for two partners of a pedigree: the recessive genes where both carry a change, counting ClinVar hits, Delly deletions over an exon, SMN1 copy loss and FMR1 and FXN expansions, the X-linked findings, and what the screen does not see. The cross-check below is the same comparison by hand, for ClinVar hits only.

### Quick Cross-Check

```bash
source versions.env   # from the repository root
PARTNER_A="sample_a"
PARTNER_B="sample_b"

# Gene symbols of each partner's ClinVar hits (step 6 writes them with ClinVar's GENEINFO)
for SAMPLE in $PARTNER_A $PARTNER_B; do
  docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENEINFO\t%INFO/CLNSIG\t[%GT]\n' \
      /genome/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf \
    > /tmp/${SAMPLE}_clinvar_genes.txt
done

# GENEINFO is SYMBOL:GeneID, several joined by |
for SAMPLE in $PARTNER_A $PARTNER_B; do
  cut -f5 /tmp/${SAMPLE}_clinvar_genes.txt | tr '|' '\n' | cut -d: -f1 | grep -v '^\.$' | sort -u \
    > /tmp/${SAMPLE}_genes.txt
done

# Confirm each per-sample gene file is non-empty before reading the intersection
wc -l /tmp/${PARTNER_A}_genes.txt /tmp/${PARTNER_B}_genes.txt

# Find genes where BOTH partners have pathogenic hits
comm -12 /tmp/${PARTNER_A}_genes.txt /tmp/${PARTNER_B}_genes.txt
```

**Before reading the result, check the `wc -l` line.** A gene file with 0 lines means that partner has no hit with a gene name, or the hits file is missing or from an older run. An empty intersection then says nothing about carrier risk. Each file should list the genes of that partner's hits; rerun step 6 if it does not.

**If both gene files have lines and the intersection is empty:** No shared recessive carrier risk was found among the small variants in ClinVar. This is the most common result. It does not cover copy-number carriers such as SMA (see below).

**If genes appear in both lists:** Check whether:
1. Both variants are in the **same gene** (not just nearby genes)
2. Both are classified as **Pathogenic** or **Likely pathogenic** (not VUS)
3. The gene follows **autosomal recessive** inheritance (check OMIM or ClinVar)
4. Both partners are **heterozygous** (carriers), not homozygous

### Common Carrier Genes (Not Usually a Concern)

These genes frequently show up in ClinVar carrier screens. Being a carrier is very common and only relevant if your partner carries the same gene:

| Gene | Condition | Carrier Frequency |
|---|---|---|
| CFTR | Cystic fibrosis | 1 in 25 (European) |
| GJB2 | Hearing loss (DFNB1) | 1 in 30 |
| HFE | Hemochromatosis | 1 in 10 (H63D), 1 in 150 (C282Y) |
| MUTYH | Colorectal cancer risk | 1 in 50 |
| HEXA | Tay-Sachs disease | 1 in 30 (Ashkenazi), 1 in 300 (general) |

Spinal muscular atrophy is not in this table on purpose. Most SMA carriers have one copy of SMN1 instead of two: a copy-number loss, not a small variant. ClinVar screening of a VCF cannot see it, so this check says nothing about SMA carrier status.

The opt-in [step 35](35-paralogs.md) (Parascopy) estimates SMN1 and SMN2 copy number from each partner's BAM. What it can and cannot say for a couple:

- **Assessed:** each partner's SMN1 copy number, with a quality. One partner with one SMN1 copy (quality 20 or more, filter `PASS`) is a likely carrier; both partners with one copy is the pattern that gives a 1 in 4 chance of an affected child.
- **Not assessed:** a "2+0" carrier, two SMN1 copies on one chromosome and none on the other, reads as two copies like a non-carrier; small variants inside SMN1; the other paralog genes (GBA, CYP21A2, HBA1/HBA2 and so on).

A result that matters for family planning needs a clinical SMN1 carrier test; ask for one that reports the 2+0 risk.

---

## Comparing Pharmacogenomics

PharmCAT results can differ dramatically between partners. Compare the HTML reports side by side for genes that affect commonly prescribed medications. HLA-A, HLA-B and CYP2D6 come from the BAM-based callers through step 36; CYP2D6 has a result only when pypgx and Cyrius (opt-in) agree, so a partner whose CYP2D6 reads `indeterminate` has no comparable result, not a normal one:

| Gene | One Partner is Rapid, Other is Poor? | Clinical Impact |
|---|---|---|
| CYP2C19 | PPIs, SSRIs, clopidogrel | Different SSRI dosing needed |
| CYP2D6 | Codeine, tramadol, psych meds | Codeine dangerous for poor metabolizers |
| NAT2 | Isoniazid, hydralazine | Slow acetylators get more side effects at standard doses |
| CYP2C9 | Warfarin, NSAIDs | Warfarin dosing differs |

---

## Parent-Child Analysis

If you have WGS data for a parent and child, you can investigate:

### Inherited vs De Novo Variants

A de novo variant is one that appeared for the first time in the child (not present in either parent). These are rare (~50-100 per genome) and occasionally clinically significant.

```bash
source versions.env   # from the repository root
PARENT="parent_name"
CHILD="child_name"

# Find variants in the child that are NOT in the parent
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools isec -C \
    /genome/${CHILD}/vcf/${CHILD}.vcf.gz \
    /genome/${PARENT}/vcf/${PARENT}.vcf.gz \
    -p /genome/${CHILD}/vcf/de_novo_candidates/
```

**Note:** This is a rough screen. With both parents sequenced, [step 38](#family-analysis-step-38-opt-in) calls de novo candidates from a joint call instead. True de novo detection requires **trio analysis** (both parents + child) and careful filtering for sequencing errors. Many variants flagged by `bcftools isec` will be false positives (present in the parent but missed by the variant caller due to low coverage at that position).

### Carrier Inheritance Tracing

If the child is a carrier for a recessive condition, you can check which parent contributed the variant:

```bash
source versions.env   # from the repository root
GENE_REGION="chr13:20189473-20189473"  # Example: GJB2 position

for SAMPLE in $PARENT $CHILD; do
  echo "--- ${SAMPLE} ---"
  docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
    bcftools view -r "$GENE_REGION" /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz
done
```

---

## Structural Variant Comparison

SVs called by multiple callers in one person have lower false-positive rates. SVs shared between family members add further confidence:

```bash
source versions.env   # from the repository root
# Compare Manta SVs between two samples
# (Simple overlap check using bedtools-style comparison)
for SAMPLE in $PARTNER_A $PARTNER_B; do
  docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' \
      /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
    > /tmp/${SAMPLE}_svs.bed
done

echo "Partner A SVs: $(wc -l < /tmp/${PARTNER_A}_svs.bed)"
echo "Partner B SVs: $(wc -l < /tmp/${PARTNER_B}_svs.bed)"

# Exact position matches (most stringent)
comm -12 <(sort /tmp/${PARTNER_A}_svs.bed) <(sort /tmp/${PARTNER_B}_svs.bed) | wc -l
echo "Shared SVs (exact match)"
```

**Expected:** Partners (unrelated) share very few SVs at exact positions. Parent-child pairs share ~50% of SVs.

---

## Mitochondrial Haplogroup Comparison

| Relationship | Expected mtDNA Result |
|---|---|
| Partners | Different haplogroups (unless same maternal lineage) |
| Siblings | Identical haplogroup (same mother) |
| Mother-child | Identical haplogroup |
| Father-child | Different haplogroups (mtDNA is maternal only) |

If siblings have different mitochondrial haplogroups, it may indicate different biological mothers (adoption, etc.) or a very rare paternal mtDNA inheritance event.

The Y-chromosome haplogroup ([step 37](37-y-haplogroup.md), male samples) is the paternal counterpart: a father and his sons, and brothers, share it.

---

## Telomere Length Comparison

Telomere length (step 10) is most informative when compared between samples of **similar age** sequenced on the **same platform**:

- **Partners of similar age:** Should be roughly similar. Significant differences (>20%) may reflect lifestyle, stress, or genetic factors.
- **Parent-child:** Parent will generally have shorter telomeres. The difference correlates loosely with age gap.
- **Siblings:** Similar telomere lengths expected. Large differences may be worth investigating with a clinical telomere assay.

**Caution:** TelomereHunter's `tel_content` is a rough estimate. Only use it for relative comparisons between samples processed identically. Do not compare with values from other labs, papers, or sequencing platforms.

---

## What This Pipeline Cannot Do (Yet)

- **Read-backed or statistical phasing** (which variants are on which chromosome copy in one person): step 38 phases a child only by inheritance, from the parents' genotypes; a single person needs read-based phasing (long reads help) or statistical phasing tools like Eagle2/SHAPEIT
- **Shared segments between relatives** (identity by descent, stretch by stretch): step 38 gives one relatedness number per pair
- **Comparable polygenic risk scores**: step 25 computes raw scores per person, but turning them into percentiles that can be compared between people needs a matched reference population (see [step 25](25-prs.md))
- **Ancestry estimates**: one sample at a time. Step 26 projects each sample onto a 1000 Genomes panel and [step 26b](26b-ancestry-fine.md) estimates its admixture fractions and paints its chromosomes; neither compares two samples with each other, so read their outputs side by side

These are planned for future pipeline versions.
