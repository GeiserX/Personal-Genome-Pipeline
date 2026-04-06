# Whole Exome Sequencing (WES) Guide

If your data comes from whole exome sequencing instead of whole genome sequencing, most of this pipeline still works -- but some steps need different settings, and others should be skipped entirely. This guide explains what WES is, which pipeline steps apply, and how to run them.

---

## What WES Is (and How It Differs from WGS)

Whole exome sequencing targets the **exons** of your genome -- the ~22,000 protein-coding gene regions that make up roughly 2% of total DNA. A capture kit (set of probes) binds to these regions, pulling them out of a DNA sample for sequencing. Everything else is discarded before sequencing begins.

| | WGS (30X) | WES (100-150X) |
|---|---|---|
| **Genome covered** | ~100% | ~2% (exons only) |
| **Raw data size** | 60-90 GB FASTQ | 5-10 GB FASTQ |
| **BAM size** | 80-120 GB | 5-15 GB |
| **Variants called** | ~4-5 million | ~30,000-50,000 |
| **Mean depth (on-target)** | 30X | 100-150X |
| **Cost** | $200-1,000 | $100-500 |
| **Detects coding SNPs/indels** | Yes | Yes (higher depth) |
| **Detects intronic/intergenic variants** | Yes | No |
| **Detects structural variants** | Yes (large SVs) | Limited |
| **Detects repeat expansions** | Yes | Only within captured regions |
| **Telomere length** | Yes | No |
| **CNV calling (genome-wide)** | Yes | No |

**Bottom line:** WES is optimized for finding coding variants -- the mutations most likely to directly affect protein function. It has higher per-base depth in exonic regions than typical WGS, which can mean better sensitivity for rare coding variants. But it is completely blind to non-coding regions, structural variants, telomeres, and most repeat expansions.

---

## What You Need: The Capture BED File

WES data is only meaningful in the regions targeted by the capture kit. Outside those regions, coverage drops to near zero and any reads that map there are off-target noise.

You must provide a **BED file** defining the capture regions. This file tells the pipeline where your data actually has coverage. Without it, coverage statistics will be misleading (showing ~2X average across the whole genome instead of ~100-150X on-target) and some tools will waste time analyzing off-target regions.

### How to Set It

Pass the BED file path via the `CAPTURE_BED` environment variable:

```bash
export CAPTURE_BED="${GENOME_DIR}/reference/capture_regions.bed"
```

### Where to Get Your BED File

Contact your sequencing provider or download from the capture kit vendor's website:

| Capture Kit | Vendor | Exonic Targets | BED Download |
|---|---|---|---|
| SureSelect Human All Exon V8 | Agilent | ~35 Mb | [Agilent SureDesign](https://earray.chem.agilent.com/suredesign/) (free account) |
| SureSelect Clinical Research Exome V3 | Agilent | ~54 Mb | Same portal |
| Twist Human Core Exome | Twist Bioscience | ~33 Mb | [Twist support](https://www.twistbioscience.com/resources/data-files) |
| xGen Exome Hyb Panel v2 | IDT | ~34 Mb | [IDT download page](https://www.idtdna.com/pages/products/next-generation-sequencing/hybridization-capture/lockdown-panels/xgen-exome-panel) |
| KAPA HyperExome | Roche | ~43 Mb | [Roche Sequencing Solutions](https://sequencing.roche.com/) |
| Nextera Rapid Capture Exome | Illumina | ~37 Mb | [Illumina support](https://support.illumina.com/) |

> **Build matters.** Make sure the BED file matches your reference genome build (GRCh38/hg38). If you only have an hg19/GRCh37 BED, convert it with UCSC liftOver or Picard LiftoverVcf before using it.

> **If you cannot find your BED file:** Ask your sequencing provider directly. The lab that performed your sequencing knows exactly which capture kit was used and can provide the target regions file. Without it, coverage QC and some variant-calling adjustments cannot be applied correctly.

---

## Per-Step WES Compatibility

### Works Unchanged

These steps operate on VCF data and do not depend on genome-wide BAM coverage. They work identically for WES and WGS.

| Step | Name | Why It Works |
|---|---|---|
| **6** | ClinVar Screen | Checks your variants (whatever you have) against ClinVar |
| **7** | PharmCAT | Star allele calling from VCF genotypes |
| **8** | HLA Typing (T1K) | HLA genes are included in all major exome capture kits |
| **11** | ROH Analysis | Runs on VCF genotype data. Detects large ROH (>5 Mb) but misses smaller ones due to sparse non-coding coverage |
| **12** | Mito Haplogroup | Works if your capture kit includes mitochondrial targets (most clinical kits do) |
| **13** | VEP Annotation | Annotates whatever variants you supply |
| **17** | CPSR | Cancer predisposition screening on your VCF. All ACMG SF v3.2 genes are coding |
| **25** | PRS | Matches scoring file variants against your VCF. Expect ~5-15% variant matching (vs ~28% from WGS) because most PRS scoring files include non-coding variants |
| **27** | CPIC Recommendations | Operates on PharmCAT output |

### Needs Adjustment

These steps require different parameters when processing WES data.

| Step | Name | What Changes | Why |
|---|---|---|---|
| **3** | DeepVariant | `MODEL_TYPE=WES ./scripts/03-deepvariant.sh sample` | DeepVariant has separate trained models for WGS and WES. The WES model is calibrated for the higher depth, sharper coverage boundaries, and different error profiles of exome data. Using the WGS model on WES data inflates false positive rates at capture region boundaries. |
| **4** | Manta (SVs) | Add `--exome` flag to `configManta.py` | Without `--exome`, Manta expects genome-wide coverage and misinterprets off-target regions (zero coverage) as evidence of large deletions. The `--exome` flag restricts analysis to the capture regions and adjusts depth expectations. |
| **9** | ExpansionHunter | Limited to captured loci only | ExpansionHunter estimates repeat lengths from reads spanning the repeat region. For WES, only loci that fall within capture regions have sufficient coverage. Most clinically important repeat expansions (Huntington's HTT, Fragile X FMR1, ALS C9orf72) are in UTRs or introns and are **not captured** by standard exome kits. Expect very few usable results. |
| **16b** | mosdepth | `--by ${CAPTURE_BED}` instead of genome-wide | Calculates on-target coverage (mean depth within capture regions) rather than genome-wide coverage. Without the BED file, mosdepth reports ~2X mean coverage across the whole genome, which is technically correct but completely misleading for WES. |

### Skip These Steps

These steps require genome-wide data and produce incorrect or meaningless results with WES.

| Step | Name | Why It Fails with WES |
|---|---|---|
| **4b** | GRIDSS | Assembly-based SV caller that needs genome-wide read pairs. Off-target WES regions have no coverage, causing GRIDSS to generate thousands of false positive breakpoints at capture boundaries. |
| **10** | TelomereHunter | Estimates telomere length from telomeric repeat reads (TTAGGG). Exome capture kits do not target telomeres, so there are essentially zero telomeric reads. Output will show near-zero telomere content regardless of actual telomere length. |
| **16** | indexcov (Coverage QC) | Designed for whole-genome BAMs. Reports normalized coverage across 16kb windows, which is meaningless when 98% of the genome has zero coverage by design. Use mosdepth with `--by ${CAPTURE_BED}` instead (step 16b). |
| **18** | CNVnator | Detects CNVs from read-depth signal across the genome. WES has coverage in ~2% of the genome, making genome-wide read-depth analysis impossible. For exome CNV calling, dedicated tools like ExomeDepth or XHMM are needed (not currently in this pipeline). |
| **19** | Delly | Relies on discordant read pairs and split reads across the genome. WES read pairs that span capture boundaries look like structural variant evidence to Delly, causing massive false positive rates. |
| **20** | MToolBox (Mito variants) | GATK Mutect2 mitochondrial mode expects reads covering the full mtDNA circle. Most exome kits capture few or no mitochondrial regions, so heteroplasmy detection fails. If your kit does capture mtDNA, step 12 (haplogrep3) still works. |
| **26** | Ancestry PCA | PCA requires hundreds of thousands of evenly-distributed common variants across the genome. WES captures only coding variants, which are biased toward conserved regions and not representative of population-level genetic variation. Results will be unreliable. |

---

## Expected Variant Counts

If you have run WGS before, the variant counts from WES will look very different:

| Metric | WGS (30X) | WES (100-150X) |
|---|---|---|
| Total variants | ~4-5 million | ~30,000-50,000 |
| PASS variants | ~4.5 million | ~25,000-40,000 |
| Coding variants | ~20,000-25,000 | ~20,000-25,000 |
| ClinVar hits (pathogenic + LP) | 0-10 typical | 0-10 typical |
| PharmCAT gene calls | 20+ genes | 20+ genes |
| VEP HIGH impact | 100-200 | 80-150 |

The coding variant counts are nearly identical because WES specifically targets those regions. Where WES differs is the total count: you are missing ~99% of variants, almost all of which are in non-coding regions. For clinical interpretation of coding mutations, WES and WGS produce comparable results.

---

## Quick Start for WES Data

### 1. Gather Your Files

You need:
- **BAM or FASTQ** from your WES provider
- **Capture BED file** from your provider or the kit vendor (see table above)
- The standard pipeline reference data ([step 00](00-reference-setup.md))

Place your BED file in the reference directory:
```bash
cp /path/to/your/capture_regions.bed "${GENOME_DIR}/reference/capture_regions.bed"
```

### 2. Set Environment Variables

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_name
export DATA_TYPE=WES
export CAPTURE_BED="${GENOME_DIR}/reference/capture_regions.bed"
```

### 3. Align (if starting from FASTQ)

```bash
./scripts/02-alignment.sh $SAMPLE
```

Alignment works identically for WES and WGS data. minimap2 does not need different settings.

### 4. Call Variants

```bash
MODEL_TYPE=WES ./scripts/03-deepvariant.sh $SAMPLE
```

### 5. Run VCF-Based Analysis

These steps work unchanged:

```bash
./scripts/06-clinvar-screen.sh $SAMPLE     # ClinVar pathogenic variants
./scripts/07-pharmacogenomics.sh $SAMPLE   # PharmCAT drug-gene interactions
./scripts/13-vep-annotation.sh $SAMPLE     # Functional annotation
./scripts/17-cpsr.sh $SAMPLE               # Cancer predisposition screening
./scripts/27-cpic-lookup.sh $SAMPLE        # Drug recommendations
```

### 6. Run BAM-Based Steps (with Adjustments)

```bash
# Manta: add --exome manually until scripts support DATA_TYPE
# mosdepth: pass capture BED for on-target coverage stats
CAPTURE_BED=${GENOME_DIR}/reference/your_capture.bed ./scripts/16b-mosdepth.sh $SAMPLE
```

### 7. Skip These Steps

Do **not** run: TelomereHunter (10), indexcov (16), CNVnator (18), Delly (19), GRIDSS (4b), Mito analysis (20), Ancestry PCA (26).

---

## WES Coverage QC: What to Look For

With WGS, you check for ~30X mean genome-wide coverage. With WES, the key metrics are different:

| Metric | Good WES | Concerning | How to Check |
|---|---|---|---|
| **Mean on-target depth** | >80X | <50X | mosdepth with `--by` BED |
| **Fraction of targets at 20X** | >95% | <85% | mosdepth thresholds file |
| **Fraction of targets at 10X** | >98% | <90% | mosdepth thresholds file |
| **On-target rate** | >70% | <50% | On-target reads / total reads |
| **Duplication rate** | <20% | >30% | fastp or Picard MarkDuplicates |

Low on-target rates mean the capture was inefficient and many reads are wasted on non-target regions. This does not affect variant calling in captured regions but means you paid for sequencing that produced no usable data.

---

## WES Vendors and Capture Kits

Most clinical and research WES providers use one of a few standard capture kits. Knowing which kit was used helps you find the correct BED file:

| Provider | Typical Kit | Notes |
|---|---|---|
| Invitae, GeneDx (clinical) | Custom panels | Usually target specific gene lists, not the full exome |
| Ambry Genetics | Agilent SureSelect CRE | Clinical-grade exome |
| Blueprint Genetics | Twist or IDT | Research and clinical |
| Novogene, BGI (research) | Agilent SureSelect V6/V7/V8 | Common for research WES |
| Macrogen | IDT xGen or Agilent | Varies by order |
| Illumina (TruSeq/Nextera) | Nextera Rapid Capture or Illumina Exome Panel | Older kits; BED files on Illumina support site |

> **Clinical vs research WES:** Clinical exome labs often use curated gene panels that cover fewer regions than a full research exome kit. If your report says "clinical exome" or names specific gene panels, the target region may be smaller than a standard exome capture. The pipeline still works -- just expect fewer total variants.

---

## Limitations

1. **No structural variant calling.** WES cannot reliably detect deletions, duplications, inversions, or translocations. The exome SV callers that exist (ExomeDepth, XHMM, CoNIFER) detect exon-level CNVs from read depth, but they are not currently in this pipeline.

2. **No telomere or repeat expansion analysis.** These require reads from non-coding regions that WES does not capture.

3. **Reduced PRS accuracy.** Most polygenic risk scoring files include many non-coding variants. WES typically matches only 5-15% of a scoring file's variants (vs ~28% from WGS). Imputation from exome data does not help much because the target regions are too sparse for accurate imputation of flanking non-coding variants.

4. **ROH detection is less precise.** Runs of homozygosity analysis on WES data detects large ROH segments (>5 Mb) reasonably well because exome targets are spread across all chromosomes. Smaller ROH segments may be missed due to gaps between captured regions.

5. **Ancestry estimation does not work.** PCA-based ancestry requires uniformly-distributed variants across the genome, which exome data cannot provide.

6. **Deep intronic and regulatory variants are invisible.** Some clinically significant variants are in deep intronic regions (e.g., splicing regulators) or promoters/enhancers. WES will not detect these. If you have a clinical suspicion that WES missed, WGS is the appropriate follow-up.

---

## Should You Get WES or WGS?

If you are choosing between the two:

| Factor | WES Advantage | WGS Advantage |
|---|---|---|
| **Cost** | Cheaper ($100-500) | More comprehensive ($200-1,000) |
| **Storage** | Smaller files (5-15 GB BAM) | Much larger (80-120 GB BAM) |
| **Coding variant detection** | Higher per-base depth in exons | Same variants, slightly lower depth |
| **Non-coding variants** | None | Full coverage |
| **Structural variants** | Very limited | Comprehensive |
| **Pharmacogenomics** | Same (PGx genes are coding) | Same |
| **Repeat expansions** | Very limited | Full |
| **ClinVar screening** | ~90% of ClinVar pathogenic variants are coding | Catches the other ~10% too |
| **Future reanalysis value** | Limited to exons | Full genome for future discoveries |

**If you already have WES data:** This pipeline extracts substantial value from it. Run the compatible steps and you will get ClinVar screening, pharmacogenomics, cancer predisposition, and functional annotation that is comparable to WGS for coding regions.

**If you are ordering new sequencing:** WGS is generally recommended because the price difference has shrunk dramatically while WGS captures everything WES does plus non-coding regions, structural variants, and repeat expansions. The extra data storage is the main downside.
