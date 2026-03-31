# Interpreting Your Results

You've run the pipeline. Now you have directories full of VCFs, TSVs, and HTML reports. This guide explains what to look at first and what it all means — no bioinformatics degree required.

## Before You Panic: What Every Genome Looks Like

If this is your first time looking at your own genomic data, the numbers can be alarming. Here is what a **completely normal, healthy person's genome** looks like:

| Finding | Normal Range | Why It Seems Scary |
|---|---|---|
| Total variants | 4.5-5.5 million | Sounds like millions of "mutations" — but >99.9% are normal human variation |
| ClinVar pathogenic hits | 0-10 | Everyone carries a few recessive disease variants. You need TWO copies (one from each parent) to be affected |
| HIGH impact variants (VEP) | 100-150 | Most are heterozygous in non-essential genes. Having one broken copy is fine. |
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
- You are **not affected** and will never develop the condition
- The only relevance is for **family planning**: if your partner carries the same gene, each child has a 25% chance of being affected
- **Exception — MUTYH**: heterozygous carriers have a modestly elevated colorectal cancer risk (~2x population risk). MUTYH biallelic carriers have a much higher risk, but even single-carrier status warrants earlier colonoscopy screening (discuss with your doctor)
- Examples: GJB2 (hearing loss), CFTR (cystic fibrosis), HFE (hemochromatosis)

---

## Start Here: The Three Most Important Outputs

### 1. ClinVar Screen (Step 6)

**What it tells you:** Known pathogenic variants in your genome, as classified by ClinVar (NCBI's public database of clinically significant variants).

**Where to look:** `${SAMPLE}/clinvar/`

**How to read it:**
- Each line in the output VCF is a variant in your genome that matches a known ClinVar entry
- The `CLNSIG` field tells you the clinical significance: `Pathogenic`, `Likely_pathogenic`, `Benign`, `Likely_benign`, `Uncertain_significance`
- **Focus on `Pathogenic` and `Likely_pathogenic` only**

**What to expect:**
- Everyone has 50-100+ ClinVar matches, most are benign or uncertain
- 0-5 pathogenic/likely pathogenic hits is typical
- Most pathogenic hits are **carrier status** (heterozygous) for recessive conditions — you carry one copy but aren't affected
- A heterozygous pathogenic variant in a recessive gene (like GJB2 for hearing loss) means you're a **carrier**, not affected
- A homozygous pathogenic variant, or a heterozygous variant in a dominant gene, needs attention

**When to worry:**
- Pathogenic variant in a **dominant** gene (one copy is enough to cause disease)
- **Two** pathogenic variants in the same recessive gene (one from each parent)
- Any variant in cancer predisposition genes (BRCA1, BRCA2, MLH1, MSH2, etc.)

### 2. PharmCAT Report (Step 7)

**What it tells you:** How your genes affect drug metabolism. This is immediately actionable — it can change which medications your doctor prescribes.

**Where to look:** `${SAMPLE}/vcf/` — PharmCAT writes its reports alongside the VCF. Open the HTML report in a browser.

**Key genes to check:**
| Gene | Affects | Common Impact |
|---|---|---|
| CYP2C19 | PPIs, clopidogrel, SSRIs, voriconazole | Rapid metabolizers burn through drugs too fast |
| CYP2D6 | Codeine, tramadol, tamoxifen, many psych meds | Poor metabolizers get toxic buildup |
| CYP2C9 | Warfarin, NSAIDs, phenytoin | Dose adjustment needed |
| DPYD | 5-fluorouracil (cancer drug) | Poor metabolizers can die from standard doses |
| SLCO1B1 | Statins (simvastatin, atorvastatin) | Increased myopathy risk |
| NAT2 | Isoniazid (TB), caffeine | Slow acetylators have more side effects |
| UGT1A1 | Irinotecan, atazanavir | *28/*28 = Gilbert syndrome (elevated bilirubin) |

**What to do:** Print the PharmCAT report and give it to your doctor. It's the single most actionable output of this entire pipeline.

### 3. CPSR Report (Step 17)

**What it tells you:** Cancer predisposition screening using CPSR's curated cancer gene panels (panel 0 covers ~200+ genes). This is broader than the 81-gene ACMG SF v3.2 list and focused specifically on cancer predisposition.

**Where to look:** `${SAMPLE}/cpsr/` — open the HTML report in a browser.

**How to read it:**
- Variants are classified into tiers:
  - **Tier 1:** Pathogenic / Likely pathogenic — needs clinical attention
  - **Tier 2:** Variant of Uncertain Significance (VUS) with some evidence
  - **Tier 3:** VUS with limited evidence
  - **Tier 4:** Likely benign / Benign
- Focus on **Tier 1** variants only for clinical action

---

## Structural Variants (Steps 4, 5, 15, 18, 19)

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
- **Called by 2+ callers (Manta + Delly, or Manta + CNVnator):** Lower false-positive rate
- **Called by 1 caller only:** Lower confidence, may be false positive
- **duphold DHFFC < 0.7 for deletions:** High confidence (depth drops as expected)
- **duphold DHBFC > 1.3 for duplications:** High confidence (depth rises as expected)

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

| Locus | Gene | Normal | Premutation | Full Expansion | Disease |
|---|---|---|---|---|---|
| HTT | HTT | <27 | 27-35 | >36 | Huntington's disease |
| FMR1 | FMR1 | <45 | 55-200 | >200 | Fragile X syndrome |
| ATXN1 | ATXN1 | <33 | — | >39 | Spinocerebellar ataxia 1 |
| C9orf72 | C9orf72 | <24 | — | >30 | ALS / Frontotemporal dementia |
| DMPK | DMPK | <35 | — | >50 | Myotonic dystrophy type 1 |

**FMR1 intermediate zone (45-54 repeats):** Not affected, but repeats may expand in offspring. Carriers should receive genetic counseling. Premutation (55-200) carries risk of FXTAS (males >50) and FXPOI.

**"ALL CLEAR"** means no locus exceeded its disease threshold.

---

## Telomere Length (Step 10)

### What It Means

Telomeres are protective caps at the ends of chromosomes that shorten with age. TelomereHunter measures `tel_content` — the normalized telomere read count.

### How to Interpret

- **Higher = longer telomeres** (younger biological age)
- **Lower = shorter telomeres** (older biological age)
- There is no universal "normal" range — compare between samples of similar age
- Typical `tel_content` for 30X WGS: 300-800 (varies by sequencing platform and coverage)

### Limitations

- This is a rough estimate, not a clinical telomere length measurement
- Short-read WGS underestimates telomere length compared to dedicated assays (TRF, FlowFISH)
- Useful for relative comparisons between samples run on the same platform, not absolute measurements

---

## ROH Analysis (Step 11)

### What It Means

Runs of Homozygosity (ROH) are long stretches where both copies of your DNA are identical. Everyone has some ROH, but extensive ROH can indicate:
- Parental relatedness (consanguinity)
- Uniparental disomy
- Population bottleneck effects

### How to Read

- **Total ROH > 300 Mb:** Suggests parental relatedness (first-cousin equivalent)
- **Total ROH > 100 Mb but < 300 Mb:** May indicate distant relatedness
- **Total ROH < 100 Mb with all segments < 10 Mb:** Normal for outbred populations
- **Individual ROH segments > 10 Mb:** Recent inbreeding event
- **Many small ROH segments (1-5 Mb):** Population-level background (Ashkenazi, Finnish, etc.)

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

Mitochondrial haplogroups have weak associations with some diseases (Parkinson's, diabetes, longevity), but these are population-level statistics, not individual predictions. The main clinical value is in step 20 (MToolBox), which detects disease-causing mitochondrial variants and heteroplasmy.

---

## VEP Annotation (Step 13)

### What It Adds

VEP annotates every variant with:
- **Consequence type:** missense, nonsense, synonymous, splice site, etc.
- **SIFT score:** Predicts if amino acid change is tolerated (>0.05) or damaging (<0.05)
- **PolyPhen score:** Predicts if change is benign (<0.15), possibly damaging (0.15-0.85), or probably damaging (>0.85)
- **gnomAD frequency:** How common this variant is in the general population

### gnomAD Frequency: Your Best Sanity Check

The single most useful annotation VEP adds is the **gnomAD allele frequency** — how common a variant is in the general population. gnomAD v4 contains ~807,000 exomes and ~76,000 genomes; if your VEP cache uses an older release, the counts will differ but the principle is the same.

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

For finding potentially significant variants in the annotated VCF:

```bash
# High-impact variants (loss of function: stop-gain, frameshift, splice donor/acceptor)
grep "HIGH" ${SAMPLE}_vep.vcf | grep -v "^#"

# Rare missense variants predicted damaging by both SIFT and PolyPhen
grep "missense_variant" ${SAMPLE}_vep.vcf | grep "deleterious" | grep "probably_damaging"

# Variants absent from gnomAD (novel/ultra-rare)
grep "missense_variant" ${SAMPLE}_vep.vcf | grep -v "gnomAD_AF"
```

### What "HIGH Impact" Means

VEP classifies variant impact as:
| Impact | Types | Interpretation |
|---|---|---|
| HIGH | Stop gained, frameshift, splice donor/acceptor | Likely breaks the protein |
| MODERATE | Missense, in-frame insertion/deletion | Changes the protein, may or may not matter |
| LOW | Synonymous, splice region | Probably no functional effect |
| MODIFIER | Intronic, intergenic, UTR | Usually non-functional |

**Everyone has ~100-150 HIGH impact variants.** Most are in one copy (heterozygous) of non-essential genes. Don't panic at the number.

### Quick Variant Filtering Recipes

Copy-paste these commands to extract the most clinically relevant variants. All assume your VEP-annotated VCF is at `${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf`.

```bash
VEP_VCF="${GENOME_DIR}/${SAMPLE}/vep/${SAMPLE}_vep.vcf"

# 1. Homozygous loss-of-function variants (most likely to cause disease)
grep -v "^#" "$VEP_VCF" | grep "HIGH" | grep "1/1" | head -20

# 2. Rare HIGH-impact variants (gnomAD AF < 0.1%)
#    These are the variants most likely to be clinically significant
grep -v "^#" "$VEP_VCF" | grep "HIGH" | grep -v "gnomADe_AF=0\.[0-9]" | head -20

# 3. Compound heterozygous candidates: genes with 2+ heterozygous variants
#    (potential autosomal recessive — needs manual curation)
grep -v "^#" "$VEP_VCF" | grep "0/1" | grep -oP 'SYMBOL=[^;|]+' | \
  sort | uniq -c | sort -rn | awk '$1 >= 2' | head -20

# 4. Known ACMG actionable genes (59 genes recommended for return of results)
#    Quick check if any HIGH/MODERATE variants land in these genes
ACMG_GENES="BRCA1|BRCA2|MLH1|MSH2|MSH6|PMS2|APC|MUTYH|TP53|RB1|MEN1|RET|VHL|SDHB|SDHD|TSC1|TSC2|WT1|NF2|PTEN|STK11|BMPR1A|SMAD4|CDH1|PALB2|CHEK2|ATM|NBN|BARD1|RAD51C|RAD51D|BRIP1"
grep -v "^#" "$VEP_VCF" | grep -E "HIGH|MODERATE" | grep -E "$ACMG_GENES" | head -20

# 5. PharmCAT-relevant variants not caught by step 7
#    (PharmCAT misses some alleles — check CYP2D6, DPYD, UGT1A1 manually)
grep -v "^#" "$VEP_VCF" | grep -E "CYP2D6|CYP2C19|CYP2C9|DPYD|UGT1A1|SLCO1B1|TPMT|NUDT15" | head -20
```

**Important:** These are starting points, not definitive screens. Any interesting finding should be cross-referenced with ClinVar and ideally confirmed by a second method (Sanger sequencing or a clinical lab).

---

## CNVnator Results (Step 18)

CNVnator detects **copy number variants** using read depth analysis — complementary to Manta's paired-end/split-read approach.

**Where to look:** `${SAMPLE}/cnvnator/${SAMPLE}_cnvs.txt`

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

**Multi-caller overlap:** CNVs found by both Manta AND CNVnator have lower false-positive rates. Cross-reference by checking if the same genomic region appears in both output files.

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

**Multi-caller overlap:** SVs detected by Manta + Delly + CNVnator (or any 2 of 3) have substantially lower false-positive rates. Single-caller calls, especially large ones, should be viewed with caution.

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
- m.3243A>G (MELAS) at >10% heteroplasmy
- m.8344A>G (MERRF) at >10% heteroplasmy
- Any position in MT-ATP6, MT-ND genes with AF >0.10

**Cross-reference:** Compare with step 12 (haplogrep3) — your homoplasmic variants should match your assigned haplogroup.

---

## What to Do Next

1. **Print your PharmCAT report** and bring it to your next doctor visit
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
| PharmCAT | Every few months | Step 7 (pharmacogenomics) | Pull new Docker image (`docker pull pgkb/pharmcat:latest`) |

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

Sanitized examples from a real 30X WGS run, so you know what to expect.

### Variant Calling (Step 3)

```
bcftools stats output:
SN  0  number of samples:     1
SN  0  number of records:     5560412
SN  0  number of SNPs:        4198753
SN  0  number of indels:      1361659
SN  0  number of multiallelic sites:  45231

# PASS variants only: 4,650,000-4,700,000
# Ti/Tv ratio: 2.05-2.10 (if < 1.8, something is wrong)
```

### ClinVar Screen (Step 6)

```
Pathogenic/Likely Pathogenic hits: 4

  chr2:47637270 T>C (rs80338939)      — GJB2 carrier (hearing loss, recessive)
  chr1:45331175 G>A (rs36053993)      — MUTYH carrier (CRC risk, recessive)
  chr10:124774641 C>T (rs28936670)    — ACADSB carrier (metabolic, recessive)
```

All heterozygous (0/1) = carrier status only. This is a completely normal result.

### PharmCAT (Step 7)

The HTML report will show a table like:

```
Gene        Diplotype           Phenotype              Affected Drugs
CYP2C19    *1/*17              Rapid Metabolizer       PPIs, SSRIs, clopidogrel
CYP2C9     *1/*1               Normal Metabolizer      Warfarin, NSAIDs
NAT2       *5/*6               Slow Acetylator         Isoniazid, caffeine
DPYD       *1/*1               Normal Metabolizer      5-FU (safe at standard dose)
SLCO1B1    *1a/*1a             Normal Function         Statins (standard dosing)
```

Typically 18-21 of 23 genes will have confident calls. CYP2D6 may be "Inconclusive" from short-read WGS (known limitation).

### ExpansionHunter (Step 9)

```json
{
  "LocusResults": {
    "HTT": { "Genotype": "17/19" },
    "FMR1": { "Genotype": "29" },
    "C9orf72": { "Genotype": "2/3" },
    "ATXN1": { "Genotype": "29/29" },
    "DMPK": { "Genotype": "12/13" }
  }
}
```

All values well below disease thresholds = ALL CLEAR.

### CPSR (Step 17)

The HTML report tier summary:

```
Tier 1 (Pathogenic/Likely pathogenic):    0 variants
Tier 2 (VUS with evidence):              3 variants
Tier 3 (VUS limited evidence):           21 variants
Tier 4 (Likely benign/Benign):           ~53,000 variants
```

Zero Tier 1 = ALL CLEAR for cancer predisposition. The VUS count varies widely (20-200+) and is not cause for concern.

### ROH (Step 11)

```
# Autosomal ROH > 5 MB: 0
# Total autosomal ROH: 47.3 MB (all segments < 3 MB)
# Conclusion: No evidence of parental relatedness
```

Normal outbred individual. If total ROH > 100 MB or any segment > 10 MB, investigate further.

### Telomere Length (Step 10)

```
tel_content: 553.52
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
- If VEP and ClinVar agree on pathogenicity, confidence is high
- If you find a potentially significant variant using VEP that is NOT in ClinVar, search gnomAD for its population frequency before drawing conclusions
- For the most important findings, consider running a second annotation tool as validation

The pipeline uses VEP because it is the most widely used and well-maintained tool, with direct gnomAD frequency integration. But no single tool is perfect.

---

## Important Caveats

- **This is not a clinical diagnosis.** These tools use the same algorithms as clinical labs, but the pipeline has not been clinically validated.
- **False positives exist.** Short-read WGS has limitations in repetitive regions, homologous genes (CYP2D6, HLA), and structural variants.
- **False negatives exist.** Some pathogenic variants are in regions that short reads can't cover (deep intronic, repeat expansions beyond read length, large structural variants).
- **VUS (Variants of Uncertain Significance)** are not actionable. They may be reclassified in the future as more data accumulates.
- **ClinVar classifications can change.** A variant classified as pathogenic today may be reclassified as benign (or vice versa) as new evidence emerges. Re-run step 6 periodically with updated ClinVar databases.
