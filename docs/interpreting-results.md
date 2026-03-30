# Interpreting Your Results

You've run the pipeline. Now you have directories full of VCFs, TSVs, and HTML reports. This guide explains what to look at first and what it all means — no bioinformatics degree required.

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

**Where to look:** `${SAMPLE}/pharmcat/` — open the HTML report in a browser.

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

**What it tells you:** Cancer predisposition screening across ACMG SF v3.2 genes (81 genes associated with hereditary cancer syndromes).

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

### How Many Is Normal?

A typical human genome has:
- ~5,000-10,000 structural variants total
- Most are in non-coding regions and harmless
- ~5-20 may affect genes
- 0-2 may be clinically significant

### Which Callers to Trust?

If you ran multiple SV callers:
- **Called by 2+ callers (Manta + Delly, or Manta + CNVnator):** High confidence
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

---

## What to Do Next

1. **Print your PharmCAT report** and bring it to your next doctor visit
2. **Review ClinVar pathogenic hits** — check if any are in dominant genes or if you're homozygous for recessive genes
3. **Read the CPSR HTML report** — it's designed for clinical interpretation and will highlight anything that needs attention
4. **If you find something concerning:** Don't panic. Discuss with a genetic counselor. Many "pathogenic" variants have incomplete penetrance (not everyone with the variant gets the disease)
5. **For carrier status findings:** Relevant mainly for family planning. If both partners carry the same recessive condition, each child has a 25% chance of being affected

---

## Important Caveats

- **This is not a clinical diagnosis.** These tools use the same algorithms as clinical labs, but the pipeline has not been clinically validated.
- **False positives exist.** Short-read WGS has limitations in repetitive regions, homologous genes (CYP2D6, HLA), and structural variants.
- **False negatives exist.** Some pathogenic variants are in regions that short reads can't cover (deep intronic, repeat expansions beyond read length, large structural variants).
- **VUS (Variants of Uncertain Significance)** are not actionable. They may be reclassified in the future as more data accumulates.
- **ClinVar classifications can change.** A variant classified as pathogenic today may be reclassified as benign (or vice versa) as new evidence emerges. Re-run step 6 periodically with updated ClinVar databases.
