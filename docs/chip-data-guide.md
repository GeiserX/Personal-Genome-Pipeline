# Using This Pipeline with Genotyping Array Data (23andMe, MyHeritage, AncestryDNA)

If you have raw data from a consumer genotyping service instead of whole genome sequencing, you can still run a meaningful subset of this pipeline. This guide explains what works, what doesn't, and how to get the most out of ~600K SNPs.

## What Genotyping Arrays Are (and Are Not)

| | WGS (30X) | Genotyping Array |
|---|---|---|
| **Positions tested** | ~3 billion (entire genome) | ~600,000-700,000 (preselected SNPs) |
| **Coverage** | 100% of sequenceable genome | ~0.02% |
| **Can detect** | SNPs, indels, SVs, CNVs, repeats, mito | Common SNPs only |
| **Cannot detect** | — | Rare variants, indels, SVs, STR expansions, CNVs |
| **Raw data format** | FASTQ (reads) or BAM (aligned reads) | Text file of genotypes (no reads) |
| **Typical file size** | 60-120 GB | 10-30 MB |
| **Cost** | $200-1,000 | $79-229 |

**Bottom line:** Arrays test positions that were preselected because they are common in human populations and useful for ancestry or trait prediction. They miss the rare variants that are most likely to be clinically significant. But for pharmacogenomics and polygenic risk, common variants are exactly what you need.

---

## Step 1: Download Your Raw Data

### 23andMe
1. Go to **Settings > 23andMe Data** (or the "Your DNA" section)
2. Click **Download Raw Data**
3. You get a `.txt` file (tab-separated, ~25 MB zipped)

### MyHeritage
1. Go to **DNA > Manage DNA Kits**
2. Click **Download Raw DNA data**
3. You get a `.csv` file (~15 MB zipped)

### AncestryDNA
1. Go to **Settings > DNA Membership Details**
2. Click **Download DNA Data**
3. You get a `.txt` file (tab-separated, ~15 MB zipped)

> **Important:** All three services use **GRCh37 (hg19)** coordinates. This pipeline requires **GRCh38**. The conversion steps below handle this.

---

## Step 2: Convert Raw Data to GRCh38 VCF

### Why This Is Non-Trivial

Consumer chip raw data files contain genotype calls (e.g., "AG" at a position) but do **not** tell you which allele is the reference allele on the human genome. To create a valid VCF, every position needs its REF allele looked up from a reference FASTA. Without this step, homozygous ALT calls get written as REF/REF and heterozygous calls can have ref/alt swapped — silently corrupting all downstream analysis.

The conversion requires three stages:
1. **Import** raw genotypes using plink2 (handles 23andMe/AncestryDNA format natively)
2. **Fix ref/alt** using plink2's `--ref-from-fa` with a GRCh37 reference FASTA
3. **Liftover** coordinates from GRCh37 to GRCh38

### Prerequisites

One-time downloads (~3.5 GB total):

```bash
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
mkdir -p "${GENOME_DIR}/liftover" "${GENOME_DIR}/reference_hg19"

# GRCh37 reference FASTA (needed for ref/alt resolution, ~3 GB)
wget -q -O "${GENOME_DIR}/reference_hg19/human_g1k_v37.fasta.gz" \
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
gunzip "${GENOME_DIR}/reference_hg19/human_g1k_v37.fasta.gz"

# Index the reference
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.21 \
  samtools faidx /genome/reference_hg19/human_g1k_v37.fasta

# GRCh37-to-GRCh38 liftover chain file (~500 KB)
wget -q -O "${GENOME_DIR}/liftover/hg19ToHg38.over.chain.gz" \
  "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
```

### Conversion Workflow

Place your raw data file at `${GENOME_DIR}/${SAMPLE}/raw/${SAMPLE}_raw.txt`, then:

```bash
SAMPLE=your_name
GENOME_DIR=/path/to/your/data
mkdir -p "${GENOME_DIR}/${SAMPLE}/vcf"

# --- Stage 1: Import raw genotypes and fix ref/alt ---
#
# plink2 --23file handles 23andMe and AncestryDNA tab-separated formats.
# --ref-from-fa resolves which allele is reference at each position.
# Without --ref-from-fa, the VCF will have wrong REF/ALT assignments.
#
# For MyHeritage CSV: convert to 23andMe-like TSV first (see note below).

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  pgscatalog/plink2:2.00a5.10 \
  plink2 \
    --23file "/genome/${SAMPLE}/raw/${SAMPLE}_raw.txt" "${SAMPLE}" \
    --ref-from-fa "/genome/reference_hg19/human_g1k_v37.fasta" \
    --export vcf bgz \
    --out "/genome/${SAMPLE}/raw/${SAMPLE}_hg19" \
    --chr 1-22 \
    --snps-only just-acgt \
    --output-chr chr26

# Index the hg19 VCF
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz"

# --- Stage 2: Liftover to GRCh38 ---

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar LiftoverVcf \
    I="/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz" \
    O="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    CHAIN=/genome/liftover/hg19ToHg38.over.chain.gz \
    R=/genome/reference/Homo_sapiens_assembly38.fasta \
    REJECT="/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" \
    WARN_ON_MISSING_CONTIG=true

# Index the final VCF
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

echo "Done. VCF at: ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
```

**Expected output:** A bgzipped, tabix-indexed VCF on GRCh38 coordinates with 400,000-600,000 autosomal SNPs. Some variants (~5-15%) will be rejected during liftover because they map to regions that changed between builds.

> **MyHeritage CSV format:** plink2's `--23file` expects tab-separated data with columns: rsID, chromosome, position, genotype. MyHeritage CSVs use a different column layout. Convert first: strip the header, reorder columns to match 23andMe format, and replace commas with tabs.

> **X/Y/MT chromosomes:** This workflow restricts to autosomes (`--chr 1-22`). Mitochondrial and sex chromosome SNPs from your chip are not converted. For mitochondrial haplogroup estimation from chip data, dedicated tools like [HaploGrep](https://haplogrep.i-med.ac.at/) accept raw 23andMe files directly.

### Optional: Imputation

Imputation can expand your 600K chip variants to ~40M by predicting untyped genotypes from population reference panels. This significantly improves PRS variant matching.

1. Prepare per-chromosome VCFs from the hg19 data
2. Upload to the [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/) — accepts single-sample submissions and outputs GRCh38 natively
3. Download the imputed VCF, filter to R2 > 0.3, and use as your pipeline input

> **Note on Michigan Imputation Server:** MIS may require multiple samples per job (see [step 14 docs](14-imputation-prep.md)). TOPMed is generally more accessible for single-sample chip data. Check each server's current policies before uploading.

---

## Which Pipeline Steps Work with Chip Data

### Works Well

| Step | Name | Why It Works | Notes |
|---|---|---|---|
| **6** | ClinVar screen | Checks your variants against known pathogenic entries | You'll only find pathogenic variants that happen to be on the chip. Most clinically significant rare variants will be missed. |
| **7** | PharmCAT | Pharmacogenomic star alleles from SNP genotypes | This is one of the **best uses** of chip data. Most PGx-relevant positions are common SNPs that arrays cover well. CYP2D6 will still be limited. |
| **11** | ROH analysis | Runs of homozygosity from SNP genotypes | Works, but resolution is lower (~600K markers vs 5M). Large ROH (>5 MB) are still detectable. |
| **25** | PRS | Polygenic risk scores from common variants | Works **very well** — PRS scoring files are derived from GWAS arrays that overlap heavily with consumer chips. Variant matching rate may actually be higher than with WGS. |
| **27** | CPIC lookup | Drug-gene recommendations | Works if step 7 (PharmCAT) succeeds. |

### Works with Limitations

| Step | Name | Limitation |
|---|---|---|
| **13** | VEP annotation | Runs, but annotating 600K variants is much less useful than annotating 5M. The rare, potentially significant variants are the ones arrays miss. |
| **17** | CPSR | Runs, but cancer predisposition screening on chip data has very low sensitivity. Most pathogenic variants in cancer genes are rare and not on the chip. A negative CPSR result from chip data does NOT rule out cancer predisposition. |

### Does Not Work

| Step | Name | Why |
|---|---|---|
| **2** | Alignment | No reads to align |
| **3** | DeepVariant | No BAM file |
| **4, 18, 19** | SV callers (Manta, CNVnator, Delly) | Need read-level evidence |
| **5, 15** | SV annotation (AnnotSV, duphold) | No SV calls to annotate |
| **8** | HLA typing (T1K) | Needs reads spanning HLA region |
| **9** | ExpansionHunter | Needs reads spanning repeat regions |
| **10** | TelomereHunter | Needs telomeric reads |
| **12** | Mito haplogroup | The conversion workflow produces autosomal-only VCF. For mt haplogroup from chip data, use [HaploGrep](https://haplogrep.i-med.ac.at/) directly with your raw file. |
| **16** | Coverage QC (indexcov) | No BAM to assess coverage |
| **20** | Mito variant calling (Mutect2) | Needs BAM |
| **21** | CYP2D6 (Cyrius) | Needs BAM |
| **22** | SV consensus merge | No SV calls |
| **23** | Clinical filter | Requires VEP-annotated VCF with gnomAD. Limited value on chip data. |
| **26** | Ancestry PCA | The current step 26 implementation requires >=2 samples for PCA and produces no output for a single sample. For ancestry from chip data, use the provider's built-in ancestry tools or upload to a service like [DNA Painter](https://dnapainter.com/). |

---

## What You Can Realistically Learn

### From chip data alone (no imputation)

| Analysis | Usefulness | Confidence |
|---|---|---|
| **Pharmacogenomics** | High | Good — most PGx SNPs are on the chip |
| **Carrier screening** | Low-moderate | Only finds carriers for variants on the chip |
| **Cancer predisposition** | Very low | Most pathogenic variants are rare and NOT on the chip |
| **Polygenic risk scores** | Moderate | Reasonable, but fewer matched variants than WGS or imputed data |
| **Ancestry** | N/A via pipeline | Use your provider's ancestry tools or HaploGrep for mt haplogroup |

### From chip data + imputation

| Analysis | Usefulness | Confidence |
|---|---|---|
| **PRS** | High | Comparable to WGS for well-imputed regions |
| **ClinVar screening** | Moderate | Imputed rare variants have lower confidence (check R2 scores) |
| **PharmCAT** | High | Same as chip alone (PGx SNPs are directly genotyped) |

### What you cannot learn regardless of imputation

- Structural variants (deletions, duplications, inversions)
- Repeat expansions (Huntington's, Fragile X, ALS)
- Novel rare variants not in the imputation reference panel
- Telomere length
- Mitochondrial heteroplasmy levels

---

## Recommended Workflow for Chip Data

```
Download raw data from 23andMe / MyHeritage / AncestryDNA
         |
         v
Convert to GRCh38 VCF (plink2 + liftover, see above)
         |
         +---> Step 7:  PharmCAT -----> Step 27: CPIC drug-gene report
         |
         +---> Step 25: PRS (polygenic risk scores)
         |
         +---> Step 6:  ClinVar screen (limited but useful)
         |
         +---> Step 11: ROH analysis
         |
    (optional)
         |
         v
Impute via TOPMed server
         |
         +---> Re-run step 25 with imputed VCF (better PRS)
```

**Minimum useful run** (takes ~15 minutes): Steps 7 + 27 (pharmacogenomics + drug recommendations). This is the single highest-value analysis you can do with chip data.

---

## Limitations to Keep in Mind

1. **A clean ClinVar screen does NOT mean you have no pathogenic variants.** The chip only tests ~600K positions out of 3 billion. The vast majority of known pathogenic variants are rare and not on the chip.

2. **PRS from chip data is reasonable but less precise** than PRS from WGS. The scoring files may reference variants that aren't on your chip and can't be imputed.

3. **PharmCAT coverage depends on chip version.** 23andMe v5 covers most CYP2C19, CYP2C9, and DPYD positions. Older chip versions (v3, v4) cover fewer PGx SNPs. PharmCAT will report "Not called" for genes without sufficient data.

4. **Imputed genotypes are predictions, not observations.** Imputation accuracy varies by ancestry and local linkage disequilibrium. Common variants (MAF > 5%) impute well. Rare variants (MAF < 1%) impute poorly and should not be used for clinical decisions.

5. **If you find something concerning, get WGS.** Chip data is a screening tool. Any significant finding should be confirmed with either WGS or targeted Sanger sequencing through a clinical lab.

---

## Cost Comparison: Getting More from Your Data

| Approach | Cost | Variants | Best For |
|---|---|---|---|
| Chip alone | $0 (already have data) | ~600K | PharmCAT, basic PRS |
| Chip + imputation | $0 (free servers) | ~40M (imputed) | Better PRS |
| Budget WGS (Novogene) | $200-400 | ~5M (observed) | Full pipeline |
| Standard WGS (Nebula, Dante) | $300-600 | ~5M (observed) | Full pipeline + data retention |

If you are serious about genomic analysis, WGS is worth the investment. But if you already have chip data and want to start exploring, the pharmacogenomics and PRS steps deliver real value today.
