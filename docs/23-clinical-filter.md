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

Every tier starts from the PASS records. "Rare" means VEP's `MAX_AF` (the highest allele frequency in any 1000 Genomes, gnomAD exome or gnomAD genome population) is below 1%. Without `MAX_AF` the step uses `gnomADe_AF` and `gnomADg_AF`: rare when the values it has are below 1%. A variant common in gnomAD genomes but absent from the exomes is therefore not rare. When the VEP output has none of these fields (VEP run without `--everything`, `--max_af` or `--af_gnomadg`), no tier is filtered by frequency, every MODERATE variant is kept, and the step prints a notice saying so.

A variant with no frequency value at all is **frequency unknown**, not rare. VEP gives a variant a frequency only when it finds the same allele in gnomAD, written the same way. A multi-base indel written right-shifted, or inside a multiallelic record, finds nothing and gets no value even when it is common. Step 13 normalises its input for that reason, but a VCF annotated elsewhere may not be. HIGH and MODERATE variants of unknown frequency go to their own tier (`_freq_unknown.vcf.gz`) and their own count, never into the rare HIGH or MODERATE counts. The CADD, SpliceAI and REVEL/AlphaMissense tiers keep both rare and unknown-frequency variants; the summary's `FREQ_TIER` column says which.

The gene, impact and consequence of a variant come from one transcript: its MANE Select transcript (VEP's `MANE_SELECT`, written with `--everything`), else the canonical one (`CANONICAL=YES`), else the most severe over every transcript. A variant with a frameshift on a minor transcript and an intronic consequence on the MANE Select one is ranked intronic, so it is not in the rare HIGH tier. The summary keeps the worst consequence beside the ranked one. `bin/csq_rank.awk` applies this rule, and the Nextflow module and the reports use the same one. Step 31 keeps the worst consequence for its tiers and shows the MANE Select one beside it.

### Genotype check

Each variant's call is checked against the sample's `GQ`, `DP` and `AD`. A call is **weak** when GQ is below 30 or DP below 10 (step 6's floor), when a heterozygote's alt reads are outside 20% to 80% of the depth, or when a homozygous alt call's alt reads are below 80%. A 1/1 call at GQ 11 with 12 reference and 5 alt reads is weak: the reads support a heterozygote or no call, not a homozygote. Weak calls stay in their tier and are marked (`CALL_CHECK`, `CALL_WHY`); the step prints how many there are, and the reports list them apart from the counted findings. The floors are settable: `CLINICAL_MIN_GQ`, `CLINICAL_MIN_DP`, `CLINICAL_HET_AB_MIN`, `CLINICAL_HET_AB_MAX`, `CLINICAL_HOM_AB_MIN` for the script, and `ext.min_gq`, `ext.min_dp`, `ext.het_ab_min`, `ext.het_ab_max`, `ext.hom_ab_min` on the `CLINICAL_FILTER` process. The clinical VCF records the floors it used (`##PGP_CALL_FLOOR`), and the reports read them from there.

### Rare HIGH impact
Stop-gained, frameshift, splice donor/acceptor, start-lost, on the ranked transcript.

### Rare MODERATE impact
Missense variants and in-frame insertions/deletions, on the ranked transcript.

### HIGH or MODERATE impact, frequency unknown
The same impacts with no frequency value. Their own tier, not rare (see above).

### ClinVar pathogenic/likely pathogenic, at any frequency
- Preferred source: the step 6 hits file, built from the ClinVar file in `clinvar/` that `setup.sh` refreshes. The tier holds the records at those positions.
- Step 6 matches on a split, left-aligned copy of the sample, and this tier selects the VEP records at the same CHROM and POS. An SNV always matches. An indel matches only when the caller already wrote it left-aligned, as DeepVariant does; an indel whose position moves on left-alignment is missing from this tier, though it stays in the step 6 hits and in both reports' ClinVar section.
- Next, when step 6 has not run: `ClinVar_CLNSIG`, the same ClinVar file step 13 annotated with `--custom` (pathogenic or likely pathogenic, not conflicting). It follows a ClinVar refresh once step 13 runs again.
- Last: VEP's `CLIN_SIG`, from the VEP cache release, so a ClinVar refresh never reaches it. The step says which source it used.
- A common pathogenic allele (for example HFE p.C282Y) stays: this tier has no frequency filter.

### Rare high CADD (step 30)
CADD PHRED >= 20 for variants that are not HIGH or MODERATE on the ranked transcript. Rare or of unknown frequency.

### Rare cryptic splice (step 30)
A SpliceAI delta score >= 0.2 for any gene of the value. SpliceAI writes one entry per gene, joined by commas; every entry is tested.

### Rare deleterious missense (step 30)
REVEL >= 0.644 (ClinGen's PP3 Supporting threshold) or AlphaMissense >= 0.564 (AlphaMissense's own likely_pathogenic class boundary, not an ACMG evidence level).

## Output

| File | Contents |
|---|---|
| `${SAMPLE}_clinical.vcf.gz` | All tiers merged |
| `${SAMPLE}_clinical_summary.tsv` | One row per variant: `CHROM`, `POS`, `REF`, `ALT`, `GT`, `IMPACT`, `GENE`, `Consequence`, `MAX_AF`, `CADD_PHRED`, `REVEL`, `AM_CLASS`, `TRANSCRIPT`, `WORST_IMPACT`, `WORST_Consequence`, `FREQ_TIER`, `GQ`, `DP`, `AD`, `CALL_CHECK`, `CALL_WHY`, and with the constraint table `LOEUF`, `pLI`, `mis_z` |
| `${SAMPLE}_high_impact.vcf.gz` | Rare HIGH impact |
| `${SAMPLE}_rare_moderate.vcf.gz` | Rare MODERATE impact |
| `${SAMPLE}_freq_unknown.vcf.gz` | HIGH or MODERATE impact of unknown frequency (when the VEP output has a frequency field) |
| `${SAMPLE}_clinvar_pathogenic.vcf.gz` | ClinVar P/LP (when step 6 hits, `ClinVar_CLNSIG` or `CLIN_SIG` exist) |
| `${SAMPLE}_cadd_high.vcf.gz`, `${SAMPLE}_spliceai_high.vcf.gz`, `${SAMPLE}_missense_deleterious.vcf.gz` | The step 30 tiers, when their scores exist |

`IMPACT`, `GENE` and `Consequence` are those of the ranked transcript; `TRANSCRIPT` says which set it came from (`MANE`, `CANONICAL` or `ANY`), and `WORST_IMPACT` and `WORST_Consequence` give the most severe consequence over every transcript. `GENE` is `.` for an intergenic variant. `FREQ_TIER` is `rare`, `unknown`, `common` (a ClinVar record) or `none` (no frequency field). `CALL_CHECK` is `ok`, `weak` or `unknown` (no GQ, DP or AD to judge), and `CALL_WHY` names what made a call weak, for example `GQ11,AB0.29`. A score, frequency or FORMAT field the input does not carry is written as `.`. The records of every tier carry the same marks as `PGP_*` INFO tags in the VCFs.

The constraint columns come from `bin/constraint_join.awk`, the loader step 31 and the Nextflow slivar module run too: only canonical transcripts count, the Ensembl row wins over the RefSeq one, and `mis_z` is gnomAD v4.1's `mis.z_score`. When the table is present and rows carry gene symbols but not one matches it, the step fails instead of writing `.` everywhere.

The Nextflow `CLINICAL_FILTER` module applies the same tiers, ranking, frequency tier and genotype check (it runs `bin/csq_rank.awk` too) and writes the same summary columns. It takes the ClinVar tier from `ClinVar_CLNSIG` (the `--clinvar` file, which the VEP module adds with `--custom`), else from VEP's `CLIN_SIG`, and adds no constraint columns.

### Secondary-findings genes (ACMG SF v3.3)

Both reports list, beside the clinical filter's counts, the variants in the 84 genes of the ACMG SF v3.3 list (Lee et al., Genet Med 2025; the genes in which the ACMG recommends reporting pathogenic variants found by chance, such as BRCA1, LDLR, MYBPC3 or TTN): the step 6 ClinVar hits in those genes, and the clinical filter's records with a HIGH-impact consequence on the ranked transcript in them. Of those, the reports count as "rare HIGH impact" only the rare ones with a sound call; a weak call or an unknown frequency is listed with its mark and counted apart, and a common one is left out. The list is versioned in `bin/collect_summary.py` (`ACMG_SF_VERSION`). It is a list to review with a clinician, not a set of findings: the list's own per-gene rules (HFE homozygous p.C282Y only; BTD, CYP27A1 and others only with two variants) are not applied, and a rare HIGH-impact variant is not a ClinVar classification. Step 17 (CPSR) reports the ACMG secondary findings with its own classification.

## Runtime

2.5 to 6.5 minutes, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)) (I/O-bound, reading the large VEP VCF)

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
- Frequency filtering depends on VEP having written `MAX_AF` or the gnomAD fields, and on VEP finding the allele: an unmatched one is "frequency unknown"
- Always cross-reference findings with ClinVar and consult a professional for clinical decisions

## Notes

- No additional Docker images required — uses the same bcftools image as other steps
- Uses `bcftools +split-vep` to parse VEP's pipe-delimited CSQ annotation structurally (not grep)
- The VEP VCF is compressed and indexed automatically if needed
- PASS filter is applied to exclude low-quality variant calls
- Available CSQ subfields and INFO scores are detected from the header; a tier whose data is missing is skipped and the step says so
