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

The conversion requires two stages:
1. **Import + fix ref/alt** in a single step using `bcftools convert --tsv2vcf` with a GRCh37 reference FASTA
2. **Liftover** coordinates from GRCh37 to GRCh38

> **Why not plink?** An earlier version of this guide used plink 1.9 (`--23file`) + plink2 (`--ref-from-fa`). That pipeline silently corrupts homozygous ALT genotypes in single-sample data because plink's binary format cannot represent both alleles for monomorphic sites. `bcftools convert --tsv2vcf` reads the reference FASTA directly and handles all genotype classes correctly.

### Prerequisites

One-time downloads (~3.5 GB total, plus the GRCh38 reference from [step 00](00-reference-setup.md)):

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

> **GRCh38 reference required:** The liftover step (stage 3) needs `Homo_sapiens_assembly38.fasta` in `${GENOME_DIR}/reference/`. If you haven't set up the pipeline's reference data yet, follow [step 00 — reference setup](00-reference-setup.md) first.

### Conversion Workflow

All three vendor formats need to be converted to a tab-separated file with columns: `rsID  chromosome  position  genotype` (23andMe/AncestryDNA are already in this format). Then `bcftools convert --tsv2vcf` creates a proper VCF by looking up each position's reference allele from the FASTA.

A ready-to-use script is provided at `scripts/chip-to-vcf.sh`. You can also run the steps manually:

```bash
SAMPLE=your_name
GENOME_DIR=/path/to/your/data
mkdir -p "${GENOME_DIR}/${SAMPLE}/vcf"

# --- Pre-step: Convert MyHeritage CSV to TSV ---
# (Skip this for 23andMe/AncestryDNA — their files are already TSV)
#
# MyHeritage CSVs have quoted fields and a different header.
# Strip comments, headers, and quotes, then rearrange to TSV.

grep -v "^#" "${GENOME_DIR}/${SAMPLE}/raw/MyHeritage_raw_dna_data.csv" | \
  grep -v "^RSID" | \
  sed 's/"//g' | \
  awk -F',' '{print $1"\t"$2"\t"$3"\t"$4}' \
  > "${GENOME_DIR}/${SAMPLE}/raw/${SAMPLE}_raw.txt"

# --- Stage 1: Import genotypes + fix ref/alt (single step) ---
#
# bcftools convert --tsv2vcf reads the TSV and looks up the REF allele
# from the FASTA at each position. This correctly handles:
#   - Homozygous reference (e.g., AA where REF=A → GT 0/0)
#   - Heterozygous (e.g., AG where REF=A → GT 0/1)
#   - Homozygous ALT (e.g., AA where REF=T → GT 1/1)
# The -c flag maps columns: ID=rsid, CHROM=chromosome, POS=position, AA=genotype

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools convert --tsv2vcf "/genome/${SAMPLE}/raw/${SAMPLE}_raw.txt" \
    -f /genome/reference_hg19/human_g1k_v37.fasta \
    -s "${SAMPLE}" \
    -c ID,CHROM,POS,AA \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz"

# Add "chr" prefix (required by the liftover chain file)
printf '%s\n' $(seq 1 22) X Y MT | \
  awk '{print $1" chr"$1}' > "${GENOME_DIR}/reference_hg19/chr_rename.txt"

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools annotate \
    --rename-chrs /genome/reference_hg19/chr_rename.txt \
    "/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz"

# --- Stage 2: Liftover to GRCh38 ---

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar LiftoverVcf \
    I="/genome/${SAMPLE}/raw/${SAMPLE}_hg19_chr.vcf.gz" \
    O="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    CHAIN=/genome/liftover/hg19ToHg38.over.chain.gz \
    R=/genome/reference/Homo_sapiens_assembly38.fasta \
    REJECT="/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz" \
    WARN_ON_MISSING_CONTIG=true

# Index the final VCF
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t -f "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

echo "Done. VCF at: ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"
```

**Expected output:** A bgzipped, tabix-indexed VCF on GRCh38 coordinates with ~600,000 SNPs. Typical breakdown: ~430K hom-ref, ~107K het, ~66K hom-alt, ~1K missing. Some variants (~0.3%) are rejected during liftover; ~900 have swapped REF/ALT between builds.

> **X/Y/MT chromosomes:** All chromosomes are converted (MT is renamed to chrM to match GRCh38 convention). Chip arrays cover very few mtDNA positions, so for mitochondrial haplogroup estimation, dedicated tools like [HaploGrep](https://haplogrep.i-med.ac.at/) that accept raw 23andMe files directly will give better results.

> **23andMe / AncestryDNA:** These are already tab-separated. Skip the MyHeritage CSV conversion pre-step and place your file directly at `${GENOME_DIR}/${SAMPLE}/raw/${SAMPLE}_raw.txt`.

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
| **7** | PharmCAT | Pharmacogenomic star alleles from SNP genotypes | Calls many genes (CYP2B6, CYP4F2, DPYD, NUDT15, TPMT, SLCO1B1, UGT1A1) but **misses key genes** like CYP2C19 and VKORC1 on some chip versions due to missing positions. Expect 888+ missing PGx positions. Always compare with WGS results if available. |
| **11** | ROH analysis | Runs of homozygosity from SNP genotypes | Works, but requires the `-G30` flag (chip VCFs lack PL tags). Large ROH (>1 MB) are detectable. |
| **25** | PRS | Polygenic risk scores from common variants | Works with `no-mean-imputation` flag (single sample lacks allele frequencies). Matches ~12% of large scoring files (vs ~28% from WGS). Scores are not directly comparable to WGS scores. |
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
| **12** | Mito haplogroup | Chip arrays cover very few mtDNA positions. For mt haplogroup from chip data, use [HaploGrep](https://haplogrep.i-med.ac.at/) directly with your raw file — it accepts raw 23andMe/AncestryDNA format. |
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
| **Pharmacogenomics** | Moderate-High | Calls many genes correctly, but misses some critical ones (e.g., CYP2C19, VKORC1) depending on chip version. Some calls may differ from WGS (e.g., CYP3A5). |
| **Carrier screening** | Low-moderate | Only finds carriers for variants on the chip |
| **Cancer predisposition** | Very low | Most pathogenic variants are rare and NOT on the chip |
| **Polygenic risk scores** | Low-Moderate | Matches ~12% of large scoring files. Raw scores differ substantially from WGS and are not comparable without a matched reference cohort. |
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
Convert to GRCh38 VCF (bcftools + liftover, see above)
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

3. **PharmCAT coverage depends on chip version.** In our smoke test with MyHeritage GSA chip data, PharmCAT correctly called CYP2B6, CYP4F2, DPYD, NUDT15, TPMT, and UGT1A1 — but **failed to call CYP2C19** (25 missing positions) and **VKORC1** (1 missing position), and **miscalled CYP3A5** as \*1/\*1 instead of \*3/\*3 (4 missing positions). 23andMe v5 may cover more PGx positions. PharmCAT will report "Not called" for genes without sufficient data — this is preferable to a wrong call.

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
