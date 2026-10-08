# Step 23: Clinical Variant Filter

## What This Does

Extracts the small subset of clinically interesting variants from your VEP-annotated VCF. Instead of manually searching through 4-5 million variants, this step produces a focused list of a few hundred variants that are rare and functionally impactful, plus the known ClinVar pathogenic ones.

## Why

The biggest challenge after running VEP annotation is: "I have millions of variants, what do I look at?" This step solves that by applying conservative filters to surface the variants most likely to be medically relevant.

## Tool

bcftools + `bcftools +split-vep` plugin (parses VEP CSQ fields structurally — no grep)

## Docker Image

- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Input

- The step 30 output `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz` when it is newer than the VEP output, else the VEP output of step 13 (`${SAMPLE}_vep.vcf.gz`, or `${SAMPLE}_vep.vcf`, which is compressed first). A derived file older than its source is ignored with a notice, so a re-run of step 13 is never hidden behind an old copy.
- Optional: the step 6 hits `${GENOME_DIR}/${SAMPLE}/clinvar/${SAMPLE}_clinvar_hits.vcf` and the gnomAD v4.1 constraint table `${GENOME_DIR}/annotations/gnomad_v4.1_constraint.tsv`.

## Command

```bash
./scripts/23-clinical-filter.sh your_name
```

## What Gets Filtered

Every tier starts from the PASS records. "Rare" means VEP's `MAX_AF` (the highest allele frequency in any 1000 Genomes, gnomAD exome or gnomAD genome population) is below 1% or missing. Without `MAX_AF` the step uses `gnomADe_AF` and `gnomADg_AF`, both below 1% or missing. A variant common in gnomAD genomes but absent from the exomes is therefore not rare. When the VEP output has none of these fields (VEP run without `--everything`, `--max_af` or `--af_gnomadg`), no tier is filtered by frequency, every MODERATE variant is kept, and the step prints a notice saying so.

The gene, impact and consequence of a variant are those of its most severe consequence (`bcftools +split-vep -s worst`).

### Rare HIGH impact
Stop-gained, frameshift, splice donor/acceptor, start-lost.

### Rare MODERATE impact
Missense variants and in-frame insertions/deletions.

### ClinVar pathogenic/likely pathogenic, at any frequency
- Preferred source: the step 6 hits file, built from the ClinVar file in `clinvar/` that `setup.sh` refreshes. The tier holds the records at those positions.
- Step 6 matches on a split, left-aligned copy of the sample, and this tier selects the VEP records at the same CHROM and POS. An SNV always matches. An indel matches only when the caller already wrote it left-aligned, as DeepVariant does; an indel whose position moves on left-alignment is missing from this tier, though it stays in the step 6 hits and in both reports' ClinVar section.
- Next, when step 6 has not run: `ClinVar_CLNSIG`, the same ClinVar file step 13 annotated with `--custom` (pathogenic or likely pathogenic, not conflicting). It follows a ClinVar refresh once step 13 runs again.
- Last: VEP's `CLIN_SIG`, from the VEP cache release, so a ClinVar refresh never reaches it. The step says which source it used.
- A common pathogenic allele (for example HFE p.C282Y) stays: this tier has no frequency filter.

### Rare high CADD (step 30)
CADD PHRED >= 20 for variants that are not HIGH or MODERATE.

### Rare cryptic splice (step 30)
A SpliceAI delta score >= 0.2 for any gene of the value. SpliceAI writes one entry per gene, joined by commas; every entry is tested.

### Rare deleterious missense (step 30)
REVEL >= 0.644 (ClinGen's PP3 Supporting threshold) or AlphaMissense >= 0.564 (AlphaMissense's own likely_pathogenic class boundary, not an ACMG evidence level).

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_clinical.vcf.gz` | All tiers merged |
| `${SAMPLE}_clinical_summary.tsv` | One row per variant: `CHROM`, `POS`, `REF`, `ALT`, `GT`, `IMPACT`, `GENE`, `Consequence`, `MAX_AF`, `CADD_PHRED`, `REVEL`, `AM_CLASS`, and with the constraint table `LOEUF`, `pLI`, `mis_z` |
| `${SAMPLE}_high_impact.vcf.gz` | Rare HIGH impact |
| `${SAMPLE}_rare_moderate.vcf.gz` | Rare MODERATE impact |
| `${SAMPLE}_clinvar_pathogenic.vcf.gz` | ClinVar P/LP (when step 6 hits, `ClinVar_CLNSIG` or `CLIN_SIG` exist) |
| `${SAMPLE}_cadd_high.vcf.gz`, `${SAMPLE}_spliceai_high.vcf.gz`, `${SAMPLE}_missense_deleterious.vcf.gz` | The step 30 tiers, when their scores exist |

`GENE` is the `SYMBOL` of the worst consequence, `.` for an intergenic one. A score or frequency the input does not carry is written as `.`.

The constraint columns come from `bin/constraint_join.awk`, the loader step 31 and the Nextflow slivar module run too: only canonical transcripts count, the Ensembl row wins over the RefSeq one, and `mis_z` is gnomAD v4.1's `mis.z_score`. When the table is present and rows carry gene symbols but not one matches it, the step fails instead of writing `.` everywhere.

The Nextflow `CLINICAL_FILTER` module applies the same tiers. It takes the ClinVar tier from `ClinVar_CLNSIG` (the `--clinvar` file, which the VEP module adds with `--custom`), else from VEP's `CLIN_SIG`, and adds no constraint columns.

### Secondary-findings genes (ACMG SF v3.3)

Both reports list, beside the clinical filter's counts, the variants in the 84 genes of the ACMG SF v3.3 list (Lee et al., Genet Med 2025; the genes in which the ACMG recommends reporting pathogenic variants found by chance, such as BRCA1, LDLR, MYBPC3 or TTN): the step 6 ClinVar hits in those genes, and the clinical filter's records with a HIGH-impact consequence in them. The list is versioned in `bin/collect_summary.py` (`ACMG_SF_VERSION`). It is a list to review with a clinician, not a set of findings: the list's own per-gene rules (HFE homozygous p.C282Y only; BTD, CYP27A1 and others only with two variants) are not applied, and a rare HIGH-impact variant is not a ClinVar classification. Step 17 (CPSR) reports the ACMG secondary findings with its own classification.

## Runtime

~5-10 minutes (I/O-bound, reading the large VEP VCF)

## How to Use the Output

### Quick look at the summary

```bash
# View the most important variants
column -t ${GENOME_DIR}/${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv | head -20
```

### Cross-reference with ClinVar

```bash
source versions.env   # from the repository root
# Find which clinical variants are also in ClinVar
docker run --rm -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" \
  bcftools isec -n=2 -w1 \
    /genome/${SAMPLE}/clinical/${SAMPLE}_clinical.vcf.gz \
    /genome/clinvar/clinvar.vcf.gz \
    -Oz -o /genome/${SAMPLE}/clinical/${SAMPLE}_clinical_clinvar.vcf.gz
```

### Load in a genome browser

The `_clinical.vcf.gz` file is small enough to load in [IGV Web](https://igv.org/app/) or [gene.iobio](https://gene.iobio.io/) for visual inspection.

## Limitations

- This is a **computational filter**, not a clinical interpretation
- Some pathogenic variants are LOW impact (e.g., synonymous variants affecting splicing, regulatory variants) and will be missed by this filter
- Frequency filtering depends on VEP having written `MAX_AF` or the gnomAD fields
- Always cross-reference findings with ClinVar and consult a professional for clinical decisions

## Notes

- No additional Docker images required — uses the same bcftools image as other steps
- Uses `bcftools +split-vep` to parse VEP's pipe-delimited CSQ annotation structurally (not grep)
- The VEP VCF is compressed and indexed automatically if needed
- PASS filter is applied to exclude low-quality variant calls
- Available CSQ subfields and INFO scores are detected from the header; a tier whose data is missing is skipped and the step says so
