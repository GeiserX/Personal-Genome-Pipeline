# Multi-Sample Comparison Guide

Ran the pipeline on two or more people? This guide explains how to compare their results for carrier screening, family planning, and inherited disease investigation.

---

## Carrier Screening for Partners

The most immediately useful multi-sample analysis: checking whether both partners carry pathogenic variants in the **same recessive gene**. If so, each child has a 25% chance of being affected.

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

Spinal muscular atrophy is not in this table on purpose. Most SMA carriers have one copy of SMN1 instead of two: a copy-number loss, not a small variant. ClinVar screening of a VCF cannot see it, so this check says nothing about SMA carrier status. That needs a copy-number test of SMN1 (a clinical carrier test, or a dedicated SMN1/SMN2 caller).

---

## Comparing Pharmacogenomics

PharmCAT results can differ dramatically between partners. Compare the HTML reports side by side for genes that affect commonly prescribed medications:

| Gene | One Partner is Rapid, Other is Poor? | Clinical Impact |
|---|---|---|
| CYP2C19 | PPIs, SSRIs, clopidogrel | Different SSRI dosing needed |
| CYP2D6 | Codeine, tramadol, psych meds | Codeine dangerous for poor metabolizers |
| NAT2 | Isoniazid, caffeine | Different caffeine sensitivity |
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

**Note:** This is a rough screen. True de novo detection requires **trio analysis** (both parents + child) and careful filtering for sequencing errors. Many variants flagged by `bcftools isec` will be false positives (present in the parent but missed by the variant caller due to low coverage at that position).

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

---

## Telomere Length Comparison

Telomere length (step 10) is most informative when compared between samples of **similar age** sequenced on the **same platform**:

- **Partners of similar age:** Should be roughly similar. Significant differences (>20%) may reflect lifestyle, stress, or genetic factors.
- **Parent-child:** Parent will generally have shorter telomeres. The difference correlates loosely with age gap.
- **Siblings:** Similar telomere lengths expected. Large differences may be worth investigating with a clinical telomere assay.

**Caution:** TelomereHunter's `tel_content` is a rough estimate. Only use it for relative comparisons between samples processed identically. Do not compare with values from other labs, papers, or sequencing platforms.

---

## What This Pipeline Cannot Do (Yet)

- **True trio analysis** (proband + both parents) with de novo calling: Requires tools like GATK's `--pedigree` mode or DeNovoGear
- **Phasing** (determining which variants are on which chromosome copy): Requires statistical phasing tools like Eagle2/SHAPEIT or long-read data
- **Comparable polygenic risk scores**: step 25 computes raw scores per person, but turning them into percentiles that can be compared between people needs a matched reference population (see [step 25](25-prs.md))
- **Ancestry estimates**: step 26 prepares a SNP set shared with 1000 Genomes, but a usable estimate needs joint PCA or admixture analysis against a reference panel, which is not implemented (see [step 26](26-ancestry.md))

These are planned for future pipeline versions.
