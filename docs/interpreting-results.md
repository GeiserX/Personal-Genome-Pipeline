# Interpreting Your Results

You've run the pipeline. Now you have directories full of VCFs, TSVs, and HTML reports. This guide explains what to look at first and what it all means — no bioinformatics degree required.

## Before You Look: What the Pipeline Can Tell You

Decide what you want to know before you open the reports. A default run looks further than most people expect:

- **Secondary findings are on by default.** Step 17 runs CPSR with `--secondary_findings`. Besides the cancer genes, CPSR then reports pathogenic and likely pathogenic variants in the genes of the ACMG SF list of secondary findings. That list includes genes for inherited heart conditions that can cause sudden death (cardiomyopathies, arrhythmias), familial hypercholesterolaemia and some metabolic diseases. A finding there can be serious and actionable even though you never asked about it.
- **APOE and Alzheimer's disease.** Step 25 computes a score for late-onset Alzheimer's disease (PGS000334). Its two largest weights are the APOE variants rs429358 and rs7412, which define the ε2, ε3 and ε4 alleles, and the summary report lists that score with the others. The ε4 allele raises the risk; it does not say who will get the disease. Many people choose not to learn their APOE status. Decide before you open the step 25 output.
- **Your relatives.** You share half your DNA with each parent, child and sibling. A pathogenic variant in a dominant gene means each of them has a 50% chance of carrying it too, so a result about you is also information about them. Comparing two genomes (see [multi-sample](multi-sample.md)) can also reveal unexpected family relationships.
- **Insurance.** The rules depend on where you live. In the United States, GINA stops health insurers and employers from using genetic information, but it does not cover life, disability or long-term-care insurance. Elsewhere the rules differ. Some insurers ask whether you have had a genetic test, and a result entered in your medical record can count. Check what applies to you before you act on a finding.
- **Switching secondary findings off.** Step 17 has no setting for it. Run the command on the [step 17 page](17-cpsr.md#command) without the `--secondary_findings` line, or delete that line from `scripts/17-cpsr.sh` (Nextflow: `modules/local/cpsr/main.nf`) before the run. CPSR then reports the cancer panel only.

---

## Before You Panic: What Every Genome Looks Like

If this is your first time looking at your own genomic data, the numbers can be alarming. Here is what a **completely normal, healthy person's genome** looks like:

| Finding | Normal Range | Why It Seems Scary |
|---|---|---|
| Total variants | 4.5-5.5 million | Sounds like millions of "mutations" — but >99.9% are normal human variation |
| ClinVar pathogenic hits (step 6) | 0-10 | Step 6 screens against Pathogenic/Likely_pathogenic only. Most hits are recessive carriers — you need TWO copies to be affected |
| HIGH impact variants (VEP) | 100-150 | Most are heterozygous in non-essential genes. For recessive genes, one copy is typically tolerated (but see haploinsufficiency). |
| Structural variants | 5,000-10,000 | Most are in non-coding regions. Your parents had them too. |
| Heteroplasmic mitochondrial variants | 20-40 | Low-level heteroplasmy (<5%) is universal and age-related |
| VUS (Variants of Uncertain Significance) | 20-200+ | "Uncertain" means **not enough data yet**, not "probably bad" |

### The VUS Trap

The single biggest source of unnecessary anxiety in personal genomics is **VUS — Variants of Uncertain Significance**. These are variants where:

- There is not enough scientific evidence to classify them as either pathogenic or benign
- The vast majority will eventually be reclassified as **benign** as more data accumulates
- They are **not actionable** — no clinical decision should be made based on a VUS
- CPSR may report dozens or hundreds of VUS. This is normal and expected.

**Rule of thumb:** If a variant is classified as VUS, treat it the same as if it were not tested. Do not change screening or management based on a VUS. Check back in 1-2 years when ClinVar may have reclassified it.

**One nuance:** If you have a strong family history of a condition AND a VUS appears in the relevant high-penetrance gene (e.g., BRCA1/2, TP53, MLH1/MSH2), it may be worth mentioning to a genetic counselor — not to act on the VUS, but because the family history itself may warrant enhanced screening regardless of the variant's classification.

### ClinVar Star Ratings

Not all ClinVar classifications are equally reliable. Each entry has a **review status** indicated by stars:

| Stars | Review Status | Reliability |
|---|---|---|
| 0 | No assertion criteria | Low — submitter did not explain their reasoning |
| 1 | Single submitter, criteria provided | Moderate — one lab's interpretation |
| 2 | Two or more submitters, no conflict | Good — multiple labs agree |
| 3 | Expert panel reviewed | High — reviewed by specialists |
| 4 | Practice guideline | Highest — established clinical standard |

**Focus on 2+ star entries.** Single-submitter (1-star) pathogenic calls are sometimes reclassified. If you find a scary-looking pathogenic variant with 0-1 stars, check the ClinVar entry directly at [ncbi.nlm.nih.gov/clinvar](https://www.ncbi.nlm.nih.gov/clinvar/) — look at the "Review status" and "Last evaluated" date.

### Carrier Status Is Not Disease

The most common "pathogenic" finding in any genome is **heterozygous carrier status for recessive conditions**. This means:

- You have ONE copy of a variant that causes disease when BOTH copies are affected
- For most recessive conditions, heterozygous carriers are **not clinically affected** (though some carrier states confer subtle phenotypic effects — e.g., sickle cell trait, HFE carriers and iron loading)
- The primary relevance is for **family planning**: if your partner carries the same gene, each child has a 25% chance of being affected
- **Note — MUTYH**: Biallelic (homozygous or compound het) MUTYH carriers have a well-established high colorectal cancer risk. For **monoallelic** (single-copy) carriers, the evidence is more nuanced: some meta-analyses show a modest risk increase, but a counseling framework for moderate-penetrance CRC genes (Genetics in Medicine) notes that risk estimates for monoallelic MUTYH are conflicting and that screening recommendations (e.g., earlier colonoscopy) were historically tied to carriers with a CRC family history. Discuss with a genetic counselor, especially if you have a family history of CRC
- Examples: GJB2 (hearing loss), CFTR (cystic fibrosis), HFE (hemochromatosis)

---

## Start Here: The Three Most Important Outputs

### 1. ClinVar Screen (Step 6)

**What it tells you:** Known pathogenic variants in your genome, as classified by ClinVar (NCBI's public database of clinically significant variants).

**Where to look:** `${SAMPLE}/clinvar/`

**How to read it:**
- Each line in the output VCF is a variant in your genome that matches a known ClinVar entry
- Step 6 intersects against the **pathogenic-only** ClinVar subset (Pathogenic + Likely_pathogenic). Every hit in this output is at a position ClinVar classifies as disease-associated — benign/VUS entries are excluded at the database level
- The `CLNSIG` field confirms the classification

**What to expect:**
- 0-10 pathogenic/likely pathogenic hits is typical for a 30X WGS
- Most pathogenic hits are **carrier status** (heterozygous) for recessive conditions — you carry one copy but aren't affected
- A heterozygous pathogenic variant in a recessive gene (like GJB2 for hearing loss) means you're a **carrier**, not affected
- A homozygous pathogenic variant, or a heterozygous variant in a dominant gene, needs attention

**When to worry:**
- Pathogenic variant in a **dominant** gene (one copy is enough to cause disease)
- **Two** pathogenic variants in the same recessive gene (one from each parent)
- Any variant in cancer predisposition genes (BRCA1, BRCA2, MLH1, MSH2, etc.)

### 2. PharmCAT Report (Step 7)

**What it tells you:** How your genes affect drug metabolism. These results are clinically relevant and should be shared with your prescribing physician.

**Where to look:** `${SAMPLE}/pharmcat/` (Nextflow) or `${SAMPLE}/vcf/` (bash scripts) — PharmCAT writes its reports there. Open the HTML report in a browser.

**Key genes to check:**

| Gene | Affects | Common Impact |
|---|---|---|
| CYP2C19 | PPIs, clopidogrel, SSRIs, voriconazole | Rapid metabolizers burn through drugs too fast |
| CYP2D6 | Codeine, tramadol, tamoxifen, many psych meds | Poor metabolizers get toxic buildup |
| CYP2C9 | Warfarin, NSAIDs, phenytoin | Dose adjustment needed |
| DPYD | 5-fluorouracil (cancer drug) | Poor metabolizers can die from standard doses |
| SLCO1B1 | Statins (simvastatin, atorvastatin) | Increased myopathy risk |
| NAT2 | Isoniazid (TB), caffeine | Slow acetylators have more side effects |
| UGT1A1 | Irinotecan, atazanavir | Two \*28 alleles are associated with Gilbert syndrome (elevated bilirubin) |

**What to do:** Share the PharmCAT report with your prescribing physician or pharmacist. PharmCAT is a research tool — its authors explicitly note that missing positions, unphased input, and undetected structural variation (especially CYP2D6) can affect genotype and phenotype calls. The report is a valuable starting point for pharmacogenomic-guided prescribing, but clinical confirmation may be warranted before making medication changes, especially for high-risk drugs (DPYD, CYP2D6-dependent opioids).

### 3. CPSR Report (Step 17)

**What it tells you:** Cancer predisposition screening using CPSR's curated cancer gene panels (panel 0 covers 500+ genes). Step 17 also turns on CPSR's secondary findings: the ACMG SF list, which includes cardiac and metabolic genes outside cancer (see [Before you look](#before-you-look-what-the-pipeline-can-tell-you)).

**Where to look:** `${SAMPLE}/cpsr/` — open the HTML report in a browser.

**How to read it:**
- CPSR puts every variant in its genes into one of the five ACMG/AMP classes: **Pathogenic**, **Likely pathogenic**, **VUS** (uncertain significance), **Likely benign** and **Benign**. A variant that already has a ClinVar entry can carry ClinVar's class instead of the one CPSR computed.
- Only **Pathogenic** and **Likely pathogenic** call for clinical attention. The report lists them first, and the secondary findings in their own section.

---

## Structural Variants (Steps 4, 5, 15, 18, 19, 22)

### What Are Structural Variants?

Unlike SNPs (single letter changes), structural variants are large rearrangements:
- **Deletions (DEL):** A chunk of DNA is missing
- **Duplications (DUP):** A chunk is copied extra times
- **Inversions (INV):** A chunk is flipped backwards
- **Translocations (BND):** A chunk moved to a different chromosome
- **Insertions (INS):** New DNA inserted

### A Note About BND (Breakend) Calls

If you run Manta or Delly, you will see many **BND** calls — often hundreds or thousands. BND indicates a "breakend" where one end of a read pair maps to a different chromosome or a distant location. This sounds alarming ("translocation!") but:

- **Most BND calls are artifacts** of repetitive regions, segmental duplications, or mobile elements
- A typical genome has 1,000-3,000 BND calls from Manta and 5,000+ from Delly
- **Fewer than 5 are likely real** inter-chromosomal translocations in a healthy genome
- BND calls require **multi-caller support** (called by both Manta and Delly at overlapping breakpoints) to be considered credible
- Unless a BND disrupts a known disease gene AND is confirmed by a second caller, it can be safely ignored

### How Many Is Normal?

A typical human genome has:
- ~5,000-10,000 structural variants total
- Most are in non-coding regions and harmless
- ~5-20 may affect genes
- 0-2 may be clinically significant

### Which Callers to Trust?

If you ran multiple SV callers:
- **Called by 2+ callers (Manta + Delly, or Manta + CNVpytor):** Lower false-positive rate
- **Called by 1 caller only:** Lower confidence, may be false positive
- **duphold DHFFC < 0.7 for deletions:** High confidence (depth drops as expected)
- **duphold DHBFC > 1.3 for duplications:** High confidence (depth rises as expected)

### SV Consensus (Step 22)

`${SAMPLE}/sv_merged/${SAMPLE}_sv_consensus.vcf.gz` keeps the SVs that two or more callers found (Manta, Delly, CNVpytor and, when they ran, GRIDSS, Sniffles2 and TIDDIT), grouped by chromosome, SV type and the 1 kb window their start position falls in. The end breakpoint is not compared, and two calls a few bases apart on either side of a window edge are not grouped. Expect a few hundred records. It is the short list to read first, but it drops real SVs that only one caller found; see [step 22](22-survivor-merge.md#limitations).

### AnnotSV Output

The AnnotSV TSV (step 5) adds clinical annotations to each SV:
- `ACMG_class`: 1 (benign) to 5 (pathogenic)
- `Overlapped_CDS_percent`: How much of a gene is affected
- Focus on SVs with `ACMG_class` 4 or 5 that overlap known disease genes

---

## STR Expansions (Step 9)

### What Are Repeat Expansions?

Some regions of DNA have short sequences repeated many times (e.g., CAG CAG CAG...). When the number of repeats exceeds a threshold, it can cause disease.

### How to Read ExpansionHunter Output

The output VCF lists each tested locus with the number of repeats found. Key loci:

| Locus | Gene | Normal | Intermediate / Premutation | Pathogenic | Disease |
|---|---|---|---|---|---|
| HTT | HTT | <=26 | 27-35 (mutable normal); 36-39 (reduced penetrance) | >=40 (full penetrance) | Huntington's disease |
| FMR1 | FMR1 | <45 | 45-54 (intermediate); 55-200 (premutation) | >200 | Fragile X syndrome |
| ATXN1 | ATXN1 | <33 | — | >39 | Spinocerebellar ataxia 1 |
| C9orf72 | C9orf72 | <24 | 24-30 (gray zone, lab cutoffs vary) | Typically hundreds-thousands | ALS / Frontotemporal dementia |
| DMPK | DMPK | <35 | 35-49 (premutation) | >=50 | Myotonic dystrophy type 1 |

**HTT 36-39 repeats (reduced penetrance):** Individuals in this range may or may not develop Huntington's disease. The risk increases with repeat length but is not certain. GeneReviews classifies >=40 as full penetrance and 36-39 as reduced penetrance. Alleles of 27-35 ("mutable normal") do not cause disease but may expand in offspring.

**C9orf72 gray zone:** The exact pathogenic threshold for C9orf72 is not established. Laboratory cutoffs are discordant (JNNP 2021 review). Clearly pathogenic expansions are typically hundreds to thousands of repeats. Short-read WGS has limited ability to size very large expansions accurately.

**FMR1 intermediate zone (45-54 repeats):** Not affected, but repeats may expand in offspring. Carriers should receive genetic counseling. Premutation (55-200) carries risk of FXTAS (males >50) and FXPOI.

**"ALL CLEAR"** means no locus exceeded its clearly pathogenic threshold. Intermediate-range results should be discussed with a genetic counselor.

### Stranger Annotation (Step 9b)

Step 9b adds a `STR_STATUS` field to each locus in the VCF so you do not need to look up thresholds manually:

| STR_STATUS | Meaning |
|---|---|
| `normal` | Repeat count is within the established normal range |
| `pre_mutation` | Elevated repeat count; not currently disease-causing but carries risk of expansion in offspring or late-onset carrier effects (e.g. FXTAS, FXPOI for FMR1) |
| `full_mutation` | Repeat count exceeds the pathogenic threshold for this locus |

The `Disease`, `OMIM`, `Inheritance`, `NormalMax`, and `PathologicMin` INFO fields give the clinical context for each locus. Use a VCF viewer or `bcftools query` to extract them:

```bash
bcftools query -f '%INFO/STR_STATUS\t%INFO/Disease\t%INFO/NormalMax\t%INFO/PathologicMin\n' \
    expansion_hunter/<sample>_eh_stranger.vcf
```

`full_mutation` at any locus warrants follow-up with a clinical geneticist. Short-read WGS has limited sizing accuracy for very large expansions (>150 repeats), so confirmatory testing may be recommended.

---

## Telomere Length (Step 10)

### What It Means

Telomeres are protective caps at the ends of chromosomes that shorten with cell division. TelomereHunter measures `tel_content` — the normalized telomere read count — as a proxy for relative telomere content.

### How to Interpret

- **Higher `tel_content` = more telomeric reads** (generally correlates with longer telomeres)
- **Lower `tel_content` = fewer telomeric reads** (generally correlates with shorter telomeres)
- There is no universal "normal" range — compare between samples of similar age, sequenced on the same platform
- Typical `tel_content` for 30X WGS: 300-800 (varies by sequencing platform and coverage)

### Limitations

- **Not a clinical telomere length measurement.** TelomereHunter was developed for cancer genome analysis, not as a validated healthy-population aging assay
- **Not a "biological age" readout.** While telomere length correlates with aging at a population level, it is a rough estimate of aging rate and is not established as a clinically important standalone risk marker for individuals (see Vaiserman & Krasnienkov, "Telomere Length as a Marker of Biological Age," 2021)
- Short-read WGS systematically underestimates telomere length compared to dedicated assays (TRF, FlowFISH)
- Useful only for **relative comparisons** between samples run on the same platform — not absolute measurements and not individual health predictions

---

## ROH Analysis (Step 11)

### What It Means

Runs of Homozygosity (ROH) are long stretches where both copies of your DNA are identical. Everyone has some ROH, but extensive ROH can indicate:
- Parental relatedness (consanguinity)
- Uniparental disomy
- Population bottleneck effects

### How to Read

Add up the autosomal segments of 5 Mb or more and compare the total with the table in [step 11](11-roh-analysis.md#total-roh-and-parental-relationship), which gives the expected total for each parental relationship (about 45 Mb for second cousins, 180 Mb for first cousins). In short:

- **No segments of 5 Mb or more, or only a few:** no sign that your parents are related
- **Many small segments (1-5 Mb):** population-level background, typical of population isolates (Ashkenazi, Finnish, etc.)
- **A single segment over 10 Mb on one chromosome, with little elsewhere:** possible uniparental disomy rather than related parents

### Centromeric Artifacts

Some apparent ROH near centromeres are artifacts of low-coverage sequencing in repetitive regions. These can be ignored.

---

## Mitochondrial Haplogroup (Step 12)

### What It Means

Your mitochondrial haplogroup traces your maternal ancestry lineage. It's determined by the specific set of variants in your mitochondrial DNA (inherited only from your mother).

### Common European Haplogroups

| Haplogroup | Origin | Notes |
|---|---|---|
| H | Western Europe | Most common in Europe (~40%) |
| U | Northern/Eastern Europe | Second most common (~15%) |
| T | Near East / Mediterranean | ~10% of Europeans |
| K | Near East | ~6% of Europeans, Ashkenazi ~30% |
| J | Near East | ~8% of Europeans |
| V | Iberian Peninsula / Scandinavia | ~4% |
| I | Near East / Europe | ~3% |

### Medical Relevance

Mitochondrial haplogroups have weak associations with some diseases (Parkinson's, diabetes, longevity), but these are population-level statistics, not individual predictions. The main clinical value is in step 20 (GATK Mutect2 mitochondrial mode), which detects disease-causing mitochondrial variants and heteroplasmy.

---

## VEP Annotation (Step 13)

### What It Adds

VEP annotates every variant with:
- **Consequence type:** missense, nonsense, synonymous, splice site, etc.
- **SIFT score:** Predicts if amino acid change is tolerated (>0.05) or damaging (<0.05)
- **PolyPhen score:** Predicts if change is benign (<0.15), possibly damaging (0.15-0.85), or probably damaging (>0.85)
- **gnomAD frequencies:** How common this variant is in gnomAD's exomes (`gnomADe_AF`) and genomes (`gnomADg_AF`)

### gnomAD Frequency: Your Best Sanity Check

The single most useful annotation VEP adds is the **gnomAD allele frequency** — how common a variant is in the general population. Step 13 runs VEP with `--everything`, which turns on both `--af_gnomade` and `--af_gnomadg`, so each variant gets two frequencies: `gnomADe_AF` from gnomAD's exomes and `gnomADg_AF` from its genomes. A non-coding variant outside exome capture regions has no exome frequency but can still have a genome frequency, so check both before calling a variant absent from gnomAD. The recipes below do that.

**Key principle:** A variant that is common in healthy people is almost certainly benign, regardless of what any prediction tool says.

| gnomAD AF | Interpretation | Action |
|---|---|---|
| > 5% (0.05) | Common polymorphism | Benign. Ignore. |
| 1-5% | Low-frequency variant | Almost certainly benign |
| 0.1-1% | Uncommon | Probably benign, but check ClinVar |
| 0.01-0.1% | Rare | Worth investigating if in a disease gene |
| < 0.01% | Very rare | Potentially significant. Check ClinVar + literature |
| Absent | Novel or ultra-rare | Could be significant OR a sequencing artifact. Verify with a second method |

**If a variant is "pathogenic" in ClinVar but has gnomAD AF > 1%:** The ClinVar entry may be outdated or wrong. Truly pathogenic variants for severe diseases are almost always rare (< 0.1%) because natural selection removes them from the population.

### Filtering Strategy

Step 23 (the [clinical filter](#clinical-filter-step-23)) already applies the filters most people want. For your own queries, use the [Quick Variant Filtering Recipes](#quick-variant-filtering-recipes) below. Do not `grep` the VEP VCF: VEP writes its values by position inside the pipe-separated `CSQ` field, never as `name=value`, so a grep for `gnomAD_AF` or `SYMBOL=` matches nothing, and a grep for `HIGH` or `1/1` can match text in other fields. `bcftools +split-vep` reads the fields by name.

### What "HIGH Impact" Means

VEP classifies variant impact as:

| Impact | Types | Interpretation |
|---|---|---|
| HIGH | Stop gained, frameshift, splice donor/acceptor | Likely breaks the protein |
| MODERATE | Missense, in-frame insertion/deletion | Changes the protein, may or may not matter |
| LOW | Synonymous, splice region | Probably no functional effect |
| MODIFIER | Intronic, intergenic, UTR | Usually non-functional |

**Everyone has ~100-150 HIGH impact variants.** Most are in one copy (heterozygous) of non-essential genes. Don't panic at the number.

---

## Pathogenicity Scores (Step 30 — vcfanno)

If you ran step 30, your VCF now includes quantitative pathogenicity scores beyond VEP's qualitative SIFT/PolyPhen predictions. These scores are widely used as computational evidence in variant interpretation (see ClinGen's PP3/BP4 calibration framework).

### CADD (Combined Annotation Dependent Depletion)

Scores **all** variant types (coding, non-coding, splice, regulatory). Uses a PHRED-like scale where higher = more deleterious.

| CADD PHRED | Interpretation | Context |
|---|---|---|
| < 10 | Likely benign | Bottom 90% of genome variation |
| 10-20 | Uncertain | Top 10%, but most are still benign |
| 20-25 | Potentially deleterious | Top 1% — investigate if in a disease gene |
| 25-30 | Likely deleterious | Top 0.3% — strong candidate for pathogenicity |
| > 30 | Highly deleterious | Top 0.1% — likely damaging if in a constrained gene |

**When to use CADD:** Best for non-coding and splice-region variants where SIFT/PolyPhen don't apply. For missense variants, REVEL and AlphaMissense are more specific.

### REVEL (Rare Exome Variant Ensemble Learner)

Scores **missense variants only**. Combines 13 individual tools into a single 0-1 score. Recommended by ClinGen for ACMG PP3/BP4 evidence.

The levels below are the REVEL row of Table 2 in [Pejaver et al. 2022](https://doi.org/10.1016/j.ajhg.2022.10.013), the ClinGen calibration of PP3/BP4:

| REVEL Score | ClinGen Evidence Level | Interpretation |
|---|---|---|
| <= 0.003 | BP4_Very Strong | Very strong evidence of benign |
| > 0.003 to <= 0.016 | BP4_Strong | Strong evidence of benign |
| > 0.016 to <= 0.183 | BP4_Moderate | Moderate evidence of benign |
| > 0.183 to <= 0.290 | BP4_Supporting | Supporting evidence of benign |
| > 0.290 to < 0.644 | No evidence | Uncertain significance |
| >= 0.644 to < 0.773 | PP3_Supporting | Supporting evidence of pathogenicity |
| >= 0.773 to < 0.932 | PP3_Moderate | Moderate evidence of pathogenicity |
| >= 0.932 | PP3_Strong | Strong evidence of pathogenicity |

There is no PP3 Very Strong level for REVEL: no score reached it in the calibration.

**When to use REVEL:** First-line score for evaluating missense variants. If REVEL >= 0.644, investigate the variant seriously.

### AlphaMissense

DeepMind's protein-structure-informed **missense** classifier. Uses AlphaFold2 protein structure predictions to assess amino acid substitution impact.

| am_pathogenicity | am_class | Interpretation |
|---|---|---|
| < 0.34 | likely_benign | Predicted benign by protein structure analysis |
| 0.34-0.564 | ambiguous | Uncertain — use other evidence |
| > 0.564 | likely_pathogenic | Predicted damaging based on protein structure |

These cut-offs are AlphaMissense's own class boundaries, set by its authors. They are not ACMG/ClinGen evidence levels.

**When to use AlphaMissense:** Complements REVEL as a second opinion. If they disagree, investigate further. Agreement does not add up to stronger PP3 evidence: ClinGen's calibration recommends one tool, used genome-wide, for PP3/BP4, so picking whichever tool scores a variant highest would bias the result.

### SpliceAI

Deep learning model predicting **splice-altering** variants. Scores four types of splice disruption: acceptor gain (AG), acceptor loss (AL), donor gain (DG), donor loss (DL).

| Max Delta Score | Interpretation |
|---|---|
| < 0.2 | Unlikely to affect splicing |
| 0.2-0.5 | May affect splicing — investigate |
| 0.5-0.8 | Likely affects splicing |
| > 0.8 | Strong evidence of splice disruption |

**When to use SpliceAI:** VEP already flags canonical splice site variants (GT/AG dinucleotides). SpliceAI catches **cryptic** splice variants — intronic or exonic variants that create new splice sites or disrupt existing ones through more subtle mechanisms.

### gnomAD Gene Constraint (Step 23 summary)

These are per-gene metrics (not per-variant) added to the clinical filter summary TSV. They measure how intolerant a gene is to different types of mutations.

| Metric | Threshold | Meaning |
|---|---|---|
| LOEUF < 0.35 | Constrained for loss-of-function | Gene is intolerant to LoF mutations — a HIGH impact variant here is more likely to cause disease |
| pLI >= 0.9 | Loss-of-function intolerant | Same as LOEUF but older metric. LOEUF is preferred. |
| mis_Z > 3.09 | Constrained for missense | Gene is intolerant to missense mutations — a REVEL-high missense here is more concerning |

**Combining scores:** A rare variant (gnomAD AF < 0.01%) with CADD > 25, in a constrained gene (LOEUF < 0.35), with ClinVar pathogenic classification, is a high-confidence pathogenic finding. Any one of these alone is insufficient.

---

### Quick Variant Filtering Recipes

Copy-paste these commands to extract the most clinically relevant variants. They read step 13's output, `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf`, with `bcftools +split-vep` from the pipeline's pinned bcftools image, which picks each value out of the `CSQ` field by name (the same way step 23 does). `-s worst` keeps the most severe consequence of each variant; `-d` keeps every transcript.

```bash
# Run from the repository root, with GENOME_DIR and SAMPLE set
source versions.env
bcf() { docker run --rm -i -v "${GENOME_DIR}:/genome" -w /genome "${BCFTOOLS_IMAGE}" bcftools "$@"; }
VEP_VCF="${SAMPLE}/vep/${SAMPLE}_vep.vcf"   # relative to GENOME_DIR

# Rare: below 0.1% in gnomAD exomes and genomes, or absent from them
RARE='(gnomADe_AF="." || gnomADe_AF<0.001) && (gnomADg_AF="." || gnomADg_AF<0.001)'
FIELDS='%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Consequence\t%gnomADe_AF\t%gnomADg_AF\t[%GT]\n'

# 1. Homozygous high-impact variants (stop gained, frameshift, splice donor/acceptor)
bcf view -i 'GT="AA"' "$VEP_VCF" |
  bcf +split-vep - -s worst -i 'IMPACT="HIGH"' -f "$FIELDS"

# 2. Rare high-impact variants
bcf +split-vep "$VEP_VCF" -s worst -i "IMPACT=\"HIGH\" && ${RARE}" -f "$FIELDS"

# 3. Compound heterozygous candidates: genes with two or more rare heterozygous
#    HIGH or MODERATE variants. The phase is unknown: both can sit on the same
#    copy of the gene, so this is a list to curate, not a finding.
bcf view -i 'GT="het"' "$VEP_VCF" |
  bcf +split-vep - -s worst -i "(IMPACT=\"HIGH\" || IMPACT=\"MODERATE\") && ${RARE}" -f '%SYMBOL\n' |
  sort | uniq -c | awk '$1 >= 2' | sort -rn

# 4. HIGH or MODERATE variants in genes you name, matched on the exact symbol.
#    The list here is a few cancer genes for illustration; CPSR (step 17) already
#    reports the full ACMG secondary-findings list. For pharmacogenes, use
#    GENES="CYP2D6 CYP2C19 CYP2C9 DPYD UGT1A1 SLCO1B1 TPMT NUDT15".
GENES="BRCA1 BRCA2 MLH1 MSH2 MSH6 PMS2 APC MUTYH TP53"
bcf +split-vep "$VEP_VCF" -d -i 'IMPACT="HIGH" || IMPACT="MODERATE"' -f "$FIELDS" |
  awk -F'\t' -v genes="$GENES" 'BEGIN { n = split(genes, g, " "); for (i = 1; i <= n; i++) want[g[i]] = 1 } $5 in want' |
  sort -u

# 5. Rare missense variants that SIFT and PolyPhen both call damaging
bcf +split-vep "$VEP_VCF" -s worst \
  -i "Consequence~\"missense_variant\" && SIFT~\"^deleterious[(]\" && PolyPhen~\"^probably_damaging\" && ${RARE}" \
  -f "$FIELDS"

# 6. Missense variants absent from both gnomAD sets (novel or ultra-rare)
bcf +split-vep "$VEP_VCF" -s worst \
  -i 'Consequence~"missense_variant" && gnomADe_AF="." && gnomADg_AF="."' -f "$FIELDS"
```

If the VCF came from a VEP run without gnomAD frequencies, `+split-vep` stops with `the tag "gnomADe_AF" is not defined`: drop the `RARE` term and the two gnomAD columns from `FIELDS`. PharmCAT (step 7) misses some alleles, and CYP2D6 needs the BAM-based callers (steps 21 and 32), so recipe 4 with the pharmacogene list is a cross-check, not a replacement.

**Important:** These are starting points, not definitive screens. Any interesting finding should be cross-referenced with ClinVar and ideally confirmed by a second method (Sanger sequencing or a clinical lab).

---

## CNVpytor Results (Step 18)

CNVpytor detects **copy number variants** using read depth analysis — complementary to Manta's paired-end/split-read approach.

**Where to look:** `${SAMPLE}/cnvpytor/${SAMPLE}_cnvs.txt`

**Format:** Each line has: type, region, size, normalized_RD, e-value1, e-value2, e-value3, e-value4, q0

**What to expect:**
- 3,000-4,000 total CNVs (mostly deletions)
- 1,500-2,000 significant (e-value < 0.01)
- Calls at chromosome starts (chr1:1-10000) are telomeric artifacts — ignore them

**Filtering:**
```bash
# Significant CNVs only (e-value < 0.01)
awk '$5 < 0.01' ${SAMPLE}_cnvs.txt

# Large deletions (>100kb, potentially clinically relevant)
awk '$1 == "deletion" && $3 > 100000 && $5 < 0.01' ${SAMPLE}_cnvs.txt

# Large duplications
awk '$1 == "duplication" && $3 > 100000 && $5 < 0.01' ${SAMPLE}_cnvs.txt
```

**Multi-caller overlap:** CNVs found by both Manta AND CNVpytor have lower false-positive rates. Cross-reference by checking if the same genomic region appears in both output files.

---

## Delly Results (Step 19)

Delly is a third structural variant caller, detecting deletions, duplications, inversions, and translocations.

**Where to look:** `${SAMPLE}/delly/${SAMPLE}_sv.vcf.gz`

**Quick summary:**
```bash
# Count SVs by type
bcftools query -f '%INFO/SVTYPE\n' ${SAMPLE}_sv.vcf.gz | sort | uniq -c

# Filter PASS variants only
bcftools view -f PASS ${SAMPLE}_sv.vcf.gz | grep -cv '^#'
```

**What to expect:**
- 5,000-15,000 total SV calls
- Most are small deletions (<1kb)
- PASS filter reduces count significantly

**Multi-caller overlap:** SVs detected by Manta + Delly + CNVpytor (or any 2 of 3) have substantially lower false-positive rates. Single-caller calls, especially large ones, should be viewed with caution.

---

## Mitochondrial Variants (Step 20)

GATK Mutect2 in mitochondrial mode detects variants with heteroplasmy fractions — the proportion of your mitochondria carrying each variant.

**Where to look:** `${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz`

**Key field:** `AF` (allele fraction) indicates heteroplasmy level:

| AF Level | Meaning |
|---|---|
| >0.95 | Homoplasmic — fixed in all mitochondria (haplogroup-defining) |
| 0.10-0.95 | Heteroplasmic — clinically significant range |
| 0.03-0.10 | Low-level heteroplasmy — often age-related somatic |
| <0.03 | Near detection limit |

**What to expect:**
- 50-70 PASS variants total
- 25-35 homoplasmic (haplogroup variants)
- 25-35 low-level heteroplasmic (mostly <5%)
- Poly-C tract variants at positions 302-310 are sequencing artifacts

**When to investigate further:**
- Heteroplasmic variant at a known disease position (check [MitoMap](https://www.mitomap.org/))
- m.3243A>G (MELAS) or m.8344A>G (MERRF) at detectable heteroplasmy levels
- Any position in MT-ATP6, MT-ND genes with AF >0.10

**Important caveats about heteroplasmy thresholds:**
- There is **no single absolute heteroplasmy threshold** that determines clinical significance. Thresholds vary by variant and by tissue (ClinGen/MSeqDR mtDNA interpretation specifications)
- **Blood and saliva underrepresent heteroplasmy** for many mitochondrial diseases. WGS from saliva, a cheek swab or blood may show lower heteroplasmy levels than affected tissues (muscle, nerve). m.3243A>G in particular shows different clinical phenotypes at very different heteroplasmy levels across tissues
- The AF values from this pipeline reflect the tissue your sample came from. A low or absent heteroplasmy level there does **not** rule out clinically significant heteroplasmy in other tissues
- For any detected pathogenic mtDNA variant, discuss with a specialist who can order tissue-specific testing if warranted

**Cross-reference:** Compare with step 12 (haplogrep3) — your homoplasmic variants should match your assigned haplogroup.

---

## Somatic Variants (Step 29) [EXPERIMENTAL]

Mutect2 in tumor-only mode looks for somatic mutations -- variants acquired during your lifetime rather than inherited. Consumer WGS is usually made from saliva or a cheek swab, a mix of cheek cells and white blood cells; the main category of interest is **clonal hematopoiesis (CHIP)**, which lives in the blood cells.

**Where to look:** `${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz`

**Key field:** `AF` (allele fraction) indicates the clone size:

| AF Range | Likely Source |
|---|---|
| 0.45-0.55 | Heterozygous germline (false positive) |
| ~1.0 | Homozygous germline (false positive) |
| 0.10-0.40 | Could be a large somatic clone, mosaic, or noisy germline |
| below 0.10 | At 30X this is one to three reads: mostly noise, see below |

**What to expect:**
- Thousands of PASS calls in a healthy individual -- the vast majority are germline false positives
- At 30X a variant in 2% of the reads has less than one supporting read on average, and one in 10% about three. Only variants at roughly 10% allele fraction or more can be told apart from noise. That is an approximate rule about reads, not a clone size: a heterozygous variant at 10% allele fraction sits in about 20% of the sampled cells, and copy number and the mix of cell types in the sample shift it. Most CHIP (defined from 2% allele fraction) is invisible here, and a clean result does not rule it out
- Without a matched normal sample, germline variants that are rare in gnomAD will often pass all filters

**CHIP genes to check:** DNMT3A, TET2, ASXL1, TP53, JAK2, SF3B1, SRSF2, PPM1D, CBL. CHIP prevalence increases with age and is associated with elevated cardiovascular risk and risk of hematologic malignancies.

**Important caveats:**
- This step has a much higher false positive rate than any other step in the pipeline
- Do not interpret PASS variants as confirmed somatic without cross-referencing with the germline VCF (step 3) and gnomAD frequencies (step 13)
- If a variant is also called at ~50% AF by DeepVariant in step 3, it is almost certainly germline

---

## Reads, Alignment, Variant Calling and Coverage (Steps 1b, 2, 3, 16, 16b)

These steps check that the data is good enough for everything else. Look at them first if a later result looks strange.

- **fastp (step 1b):** `${SAMPLE}/fastq_trimmed/${SAMPLE}_fastp.html`. A few percent of reads trimmed or dropped is normal; a large share points at a library or upload problem.
- **Alignment and variant calling (steps 2 and 3):** the [example output](#variant-calling-step-3) below shows what a 30X genome looks like: about 4.5 to 5.5 million variants and a Ti/Tv ratio of 2.0 to 2.1.
- **Coverage (step 16b, mosdepth):** `${SAMPLE}/mosdepth/${SAMPLE}.mosdepth.summary.txt`. The `total` row should show a mean near the depth you paid for (about 30). Below 15, small-variant calls lose accuracy. See [step 16b](16b-mosdepth.md#interpreting-results).
- **Sex check (step 16, indexcov):** `${SAMPLE}/indexcov/` estimates the copy number of chrX and chrY from the BAM index. About 1 and 1 is XY, about 2 and 0 is XX. A result that does not match the sex you gave the pipeline means a sample swap, a mislabelled file or a real sex-chromosome difference; [step 16](16-indexcov.md#interpretation) lists the patterns. Re-check the input before reading anything else.

## HLA Typing (Step 8)

**Where to look:** `${SAMPLE}/hla_t1k/${SAMPLE}_hla_genotype.tsv`, two alleles per HLA gene.

The main use is drug safety: a few HLA alleles predict severe reactions to specific drugs, for example HLA-B\*57:01 with abacavir and HLA-B\*58:01 with allopurinol ([step 8](08-hla-typing.md#key-hla-alleles-for-drug-safety) has the list). Typing from short-read WGS is approximate. HLA-A and HLA-B also reach PharmCAT (step 36), so the PharmCAT report and the CPIC recommendations (step 27) give the drug guidance for them; the CPIC file marks them `[outside call]`. If one of those alleles appears, or is missing and you are about to take the drug, ask for a clinical HLA test; do not rely on this output for transplant matching. With `KIR=true` (opt-in) the KIR genes are typed too, in `${SAMPLE}/kir_t1k/`; no step interprets them.

## Clinical Filter (Step 23)

**Where to look:** `${SAMPLE}/clinical/${SAMPLE}_clinical_summary.tsv`, one row per variant with gene, impact, scores and gnomAD constraint.

Step 23 keeps the HIGH-impact variants, the MODERATE ones below 1% in gnomAD exomes, the ClinVar pathogenic ones and, when step 30 ran, the variants that CADD, SpliceAI, REVEL or AlphaMissense flag. Expect a few hundred rows. Most are heterozygous variants in genes that tolerate one broken copy. Read the rows in constrained genes (low LOEUF, see [gene constraint](#gnomad-gene-constraint-step-23-summary)) first, then the homozygous ones. The HIGH-impact rows have no frequency filter, so check gnomAD for each before you worry. [Step 23](23-clinical-filter.md#what-gets-filtered) gives the exact rules.

## Variant Prioritization (Step 31)

**Where to look:** `${SAMPLE}/slivar/${SAMPLE}_slivar_summary.tsv` and `${SAMPLE}/slivar/${SAMPLE}_compound_hets.tsv`.

slivar sorts the rare, damaging variants into three groups (rare HIGH, rare MODERATE with damaging scores, ClinVar pathogenic) and lists genes where you carry two such variants. Those compound-het candidates are not phased: from one genome the pipeline cannot tell whether the two variants sit on different copies of the gene (which can cause recessive disease) or on the same copy (which usually does not). Expect a thousand or more candidate pairs; nearly all are noise. A pair matters only in a gene that fits your health history, and confirming it needs a parent's DNA or long reads. See [step 31](31-slivar.md#interpretation).

## More Pharmacogenomics: Cyrius, CPIC, pypgx and the Consensus (Steps 21, 27, 32, 36)

- **CPIC lookup (step 27):** `${SAMPLE}/cpic/${SAMPLE}_cpic_recommendations.txt` turns PharmCAT's calls into the drugs with CPIC guidance. Only genes where you are not a normal metabolizer get drug entries. Genes PharmCAT could not call are listed separately at the end; their absence from the drug list does not mean normal function. Its section "Calls From Other Tools" says which of HLA-A, HLA-B and CYP2D6 reached PharmCAT from the BAM-based callers, and why the others did not.
- **PGx consensus (step 36):** `${SAMPLE}/pgx_consensus/${SAMPLE}_pgx_consensus.tsv` shows, for HLA-A, HLA-B and CYP2D6, what each caller said and what was passed to PharmCAT. CYP2D6 is passed only when pypgx and Cyrius give the same diplotype and the depth at CYP2D6 passed its check; otherwise it reads `indeterminate`.
- **pypgx (step 32):** `${SAMPLE}/pypgx/${SAMPLE}_pypgx_summary.tsv` calls 23 genes, four of them (CYP2D6, CYP2A6, GSTM1, GSTT1) from the BAM, so it sees gene deletions and duplications PharmCAT cannot. `${SAMPLE}_pharmcat_comparison.tsv`, written into the same folder by step 27 (CPIC lookup), shows where the two tools agree.
- **Cyrius (step 21, opt-in, non-commercial licence):** `${SAMPLE}/cyrius/${SAMPLE}_cyp2d6.tsv` gives a second CYP2D6 call from the BAM. Without it, CYP2D6 stays `indeterminate`.

CYP2D6 is the hard gene: a nearby pseudogene and frequent copy-number changes confuse short reads. The pipeline gives CYP2D6 drug guidance only when two callers agree; even then, take any result that would change a prescription to a pharmacist or a certified pharmacogenomics test first.

## SMN1 and SMN2 Copy Number (Step 35, opt-in)

**Where to look:** `${SAMPLE}/paralogs/${SAMPLE}_smn_copy_number.tsv`, one row per stretch of the SMN1/SMN2 locus with the copy number of the two genes together (agCN) and of each (psCN), each with a quality.

One SMN1 copy suggests SMA carrier status; two do not rule it out (two copies on one chromosome, none on the other). Trust a value only with quality 20 or more and filter `PASS`, and confirm anything that matters with a clinical SMN1 test. [Step 35](35-paralogs.md#interpreting-results) explains the columns.

## Polygenic Risk Scores (Step 25)

**Where to look:** `${SAMPLE}/prs/${SAMPLE}_prs_summary.tsv`, one raw score per condition.

These raw sums are not percentiles, probabilities or comparable between conditions, and the pipeline ships no reference population to turn them into percentiles. They are also biased low or high because the VCF leaves out the sites where you match the reference ([step 25](25-prs.md#interpreting-results) explains why). Treat them as exploratory. The Alzheimer's score includes APOE: read [Before you look](#before-you-look-what-the-pipeline-can-tell-you) first.

## Reports (Steps 24 and 28)

- **HTML report (step 24):** `${SAMPLE}/${SAMPLE}_report.html` collects the headline results of the other steps on one page, including the ClinVar hits table. It contains health findings: share it only as you would share a medical record.
- **MultiQC (step 28):** `${SAMPLE}/multiqc/multiqc_report.html` puts the QC of fastp, samtools, mosdepth and the other tools in one page. It is about data quality, not about your health.

---

## What to Do Next

1. **Share your PharmCAT report** with your prescribing physician or pharmacist
2. **Review ClinVar pathogenic hits** — check if any are in dominant genes or if you're homozygous for recessive genes
3. **Read the CPSR HTML report** — it's designed for clinical interpretation and will highlight anything that needs attention
4. **If you find something concerning:** Don't panic. Discuss with a genetic counselor. Many "pathogenic" variants have incomplete penetrance (not everyone with the variant gets the disease)
5. **For carrier status findings:** Relevant mainly for family planning. If both partners carry the same recessive condition, each child has a 25% chance of being affected

---

## Re-running with Updated Databases

Genomic databases are updated continuously. Variants classified as VUS today may be reclassified next year. Periodic re-analysis is one of the most valuable things you can do.

### What to Update and When

| Database | Update Frequency | Pipeline Steps Affected | How to Update |
|---|---|---|---|
| ClinVar | Weekly | Step 6 (ClinVar screen) | Re-download from NCBI FTP (see [00-reference-setup.md](00-reference-setup.md)) |
| VEP cache | Every 6 months | Step 13 (VEP annotation) | Download new release from Ensembl FTP |
| PCGR/CPSR data | Every 6-12 months | Step 17 (CPSR) | Download new bundle from PCGR GitHub releases |
| PharmCAT | Every few months | Step 7 (pharmacogenomics) | Bump `PHARMCAT_IMAGE` in `versions.env`, then pull it |

### Recommended Re-analysis Schedule

- **Every 6 months:** Re-run steps 6 (ClinVar) and 17 (CPSR) with updated databases. These are the fastest steps (~35 minutes total) and the most likely to have new classifications.
- **Every 12 months:** Re-run step 13 (VEP) with updated cache for new gnomAD frequencies and consequence predictions.
- **After major database releases:** ClinVar periodically reclassifies large batches of variants. Follow [@ClinVarUpdates](https://twitter.com/ClinVarUpdates) or check the NCBI blog for announcements.

### What You Do NOT Need to Re-run

- Steps 2-3 (alignment + variant calling): Your variants don't change. Only re-run if a major DeepVariant version is released with improved accuracy.
- Steps 4, 18, 19 (SV callers): Structural variant calling is compute-intensive and results don't change with database updates.
- Step 10 (telomere): Telomere content doesn't change with database updates.

---

## Example Outputs: What Correct Results Look Like

Invented examples in the real output formats, so you know what to expect. Every number and every call below is made up; none comes from a real genome.

### Variant Calling (Step 3)

```
bcftools stats output:
SN  0  number of samples:     1
SN  0  number of records:     5500000
SN  0  number of SNPs:        4200000
SN  0  number of indels:      1300000
SN  0  number of multiallelic sites:  45000

# PASS variants only: 4,650,000-4,700,000
# Ti/Tv ratio: 2.05-2.10 (if < 1.8, something is wrong)
```

### ClinVar Screen (Step 6)

```
Pathogenic/Likely Pathogenic hits: 3

  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
  chr<N>:<pos> <REF>><ALT> (rs<id>)   — <GENE> carrier (<condition>, recessive)
```

A handful of heterozygous (0/1) hits in recessive genes is the usual result and means carrier status only.

### PharmCAT (Step 7)

The HTML report will show a table like:

```
Gene        Diplotype           Phenotype              Affected Drugs
CYP2C19    *x/*y               <phenotype>             PPIs, SSRIs, clopidogrel
CYP2C9     *x/*y               <phenotype>             Warfarin, NSAIDs
NAT2       *x/*y               <phenotype>             Isoniazid, caffeine
DPYD       *x/*y               <phenotype>             5-FU
SLCO1B1    *x/*y               <phenotype>             Statins
```

The phenotype column reads Poor, Intermediate, Normal, Rapid or Ultrarapid Metabolizer (Normal or Decreased Function for transporters such as SLCO1B1).

Typically 18-21 of 23 genes will have confident calls. CYP2D6 may be "Inconclusive" from short-read WGS (known limitation).

### ExpansionHunter (Step 9)

```json
{
  "LocusResults": {
    "HTT": { "Genotype": "<repeats>/<repeats>" },
    "FMR1": { "Genotype": "<repeats>" },
    "C9orf72": { "Genotype": "<repeats>/<repeats>" },
    "ATXN1": { "Genotype": "<repeats>/<repeats>" },
    "DMPK": { "Genotype": "<repeats>/<repeats>" }
  }
}
```

Each number is the repeat count on one allele. Counts below every locus threshold = ALL CLEAR; step 9b (Stranger) marks the ones that are not.

### CPSR (Step 17)

The HTML report counts the variants in each class:

```
Pathogenic:                  0 variants
Likely pathogenic:           0 variants
VUS:                        22 variants
Likely benign:             <n> variants
Benign:                    <n> variants
Secondary findings:          0 variants
```

No Pathogenic or Likely pathogenic variant, in the cancer panel or in the secondary findings, means nothing actionable was found in those genes. The VUS count varies widely (20-200+) and is not cause for concern.

### ROH (Step 11)

```
Autosomal ROH >5MB (potential consanguinity signal):

NOTE: Centromeric ROH (chr1:125-143MB, chr9:42-60MB, chr18:15-20MB) are technical artifacts, not real.
```

No line under the heading means no segment of 5 Mb or more: no evidence of related parents. If there are some, add them up and compare the total with the table in [step 11](11-roh-analysis.md#total-roh-and-parental-relationship).

### Telomere Length (Step 10)

```
tel_content: <value>
```

No universal "normal" range — compare between samples of the same age, sequenced on the same platform.

---

## Investigating Specific Variants

Found something interesting? These free tools help you dig deeper:

### Visual Inspection

- **[IGV Web](https://igv.org/app/)** — Load your BAM file (or a region of it) to visually inspect read-level evidence for a variant. Essential for confirming structural variants and checking for sequencing artifacts.
- **[gene.iobio](https://gene.iobio.io/)** — Clinically-driven variant interrogation tool. Load your VCF and BAM, search for specific genes, and see coverage, variant calls, and population frequency in one view.
- **[UCSC Genome Browser](https://genome.ucsc.edu/)** — Search for any genomic coordinate to see the surrounding genes, conservation, regulatory elements, and known variants.

### Database Lookups

- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — Search by rsID, gene name, or genomic position. Check the review status (stars) and submission history.
- **[gnomAD](https://gnomad.broadinstitute.org/)** — Search any variant to see its population frequency across 800,000+ individuals (v4). If common in gnomAD, almost certainly benign.
- **[OMIM](https://www.omim.org/)** — The definitive catalog of genetic disorders. Search by gene name to understand what conditions it causes and the inheritance pattern.
- **[GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/)** — Expert-written disease descriptions for genetic conditions. The single best resource for understanding a specific genetic disease.

---

## Annotation Tool Disagreement

An important caveat: **different annotation tools may classify the same variant differently**.

VEP (used in step 13), SnpEff, and ANNOVAR are the three most common variant annotation tools. Studies have shown that they disagree on consequence predictions for ~5-10% of variants, particularly at splice sites and multi-transcript genes.

**What this means for you:**
- If VEP says a variant is "HIGH impact" but ClinVar says it is benign, **trust ClinVar** (human-reviewed evidence > computational prediction)
- If VEP and ClinVar agree on pathogenicity, this strengthens the interpretation
- If you find a potentially significant variant using VEP that is NOT in ClinVar, search gnomAD for its population frequency before drawing conclusions
- For the most important findings, consider running a second annotation tool as validation

The pipeline uses VEP because it is the most widely used and well-maintained tool, with direct gnomAD frequency integration. But no single tool is perfect.

---

## What This Pipeline Does Not Assess

A clean result from these steps says nothing about the following. Each needs a different test or a different kind of data.

- **Copy-number changes in genes with a near-identical copy**, except SMN1/SMN2 with the opt-in step 35. Most alpha-thalassaemia (HBA1/HBA2 deletions), GBA1 and CYP21A2 changes sit in regions where short reads cannot tell the gene from its paralog. The VCF-based screens do not see them, and without step 35 neither is SMA carrier status. Step 35 estimates SMN1 copy number but cannot see a "2+0" carrier; SMA carrier status still needs a clinical SMN1 test.
- **Mobile-element insertions.** Manta reports insertions but does not classify them as Alu, LINE-1 or SVA insertions, and no step looks for them.
- **Methylation and phasing from long reads.** The long-read branch stops at alignment, small variants and structural variants; see the [long-read guide](long-read-guide.md).
- **Mosaic copy-neutral loss of heterozygosity**, present in only some cells (common in blood with age). Step 11 sees runs of homozygosity that are in every cell; it cannot see a change carried by a fraction of them.
- **Repeat expansions outside ExpansionHunter's 31 loci**, and accurate sizing of very large expansions.

## Important Caveats

- **This is not a clinical diagnosis.** These tools use the same algorithms as clinical labs, but the pipeline has not been clinically validated.
- **False positives exist.** Short-read WGS has limitations in repetitive regions, homologous genes (CYP2D6, HLA), and structural variants.
- **False negatives exist.** Some pathogenic variants are in regions that short reads can't cover (deep intronic, repeat expansions beyond read length, large structural variants).
- **VUS (Variants of Uncertain Significance)** are not actionable. They may be reclassified in the future as more data accumulates.
- **ClinVar classifications can change.** A variant classified as pathogenic today may be reclassified as benign (or vice versa) as new evidence emerges. Re-run step 6 periodically with updated ClinVar databases.
