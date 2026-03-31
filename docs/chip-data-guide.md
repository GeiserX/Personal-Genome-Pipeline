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

This conversion has three stages: raw text to plink format, liftover from GRCh37 to GRCh38, and export to VCF.

### Prerequisites

Download the GRCh37-to-GRCh38 liftover chain file (one-time, ~500 KB):

```bash
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
mkdir -p "${GENOME_DIR}/liftover"
wget -q -O "${GENOME_DIR}/liftover/hg19ToHg38.over.chain.gz" \
  "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
```

### Conversion Script

Place your raw data file at `${GENOME_DIR}/${SAMPLE}/raw/${SAMPLE}_raw.txt` (or `.csv` for MyHeritage), then run:

```bash
SAMPLE=your_name
GENOME_DIR=/path/to/your/data
RAW_FILE="${GENOME_DIR}/${SAMPLE}/raw/${SAMPLE}_raw.txt"
OUTDIR="${GENOME_DIR}/${SAMPLE}/vcf"
mkdir -p "$OUTDIR"

# --- Stage 1: Convert raw genotypes to plink2 format ---

# Detect format and convert to a simple 5-column TSV: chr, pos, id, ref, alt, genotype
# (23andMe and AncestryDNA are tab-separated; MyHeritage is CSV)
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  python:3.11-slim \
  python3 -c "
import csv, sys, gzip, os

sample = '${SAMPLE}'
raw = '/genome/${SAMPLE}/raw/${SAMPLE}_raw.txt'
# Try .csv for MyHeritage
if not os.path.exists(raw):
    raw = '/genome/${SAMPLE}/raw/${SAMPLE}_raw.csv'

out_path = f'/genome/${SAMPLE}/raw/${SAMPLE}_cleaned.tsv'
count = 0

with open(raw) as f:
    # Skip comment lines (23andMe/AncestryDNA use # comments)
    lines = [l for l in f if not l.startswith('#')]

with open(out_path, 'w') as out:
    reader = csv.reader(lines, delimiter='\t' if '\t' in lines[0] else ',')
    header = next(reader)  # skip header row
    for row in reader:
        if len(row) < 4:
            continue
        rsid, chrom, pos, genotype = row[0], row[1], row[2], row[3]
        # Skip indels (D/I), no-calls (--), and MT/X/Y for now
        if genotype in ('--', '00', 'DD', 'DI', 'II', 'D', 'I'):
            continue
        if chrom in ('MT', 'X', 'Y', 'XY', '0'):
            continue
        # Add chr prefix if missing
        if not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        out.write(f'{chrom}\t{pos}\t{rsid}\t{genotype}\n')
        count += 1

print(f'Converted {count} autosomal SNPs')
"

# --- Stage 2: Create a minimal VCF on GRCh37 coordinates ---

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  python:3.11-slim \
  python3 -c "
sample = '${SAMPLE}'
tsv = f'/genome/${SAMPLE}/raw/${SAMPLE}_cleaned.tsv'
vcf = f'/genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf'

with open(vcf, 'w') as out:
    out.write('##fileformat=VCFv4.2\n')
    out.write('##INFO=<ID=RS,Number=1,Type=String,Description=\"rsID\">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n')
    for i in range(1, 23):
        out.write(f'##contig=<ID=chr{i}>\n')
    out.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n')

    with open(tsv) as f:
        for line in f:
            chrom, pos, rsid, geno = line.strip().split('\t')
            if len(geno) != 2:
                continue
            a1, a2 = geno[0], geno[1]
            if a1 == a2:
                # Homozygous — we don't know which is ref without a lookup,
                # so write as-is. Liftover + bcftools norm will fix.
                out.write(f'{chrom}\t{pos}\t{rsid}\t{a1}\t.\t.\tPASS\t.\tGT\t0/0\n')
            else:
                out.write(f'{chrom}\t{pos}\t{rsid}\t{a1}\t{a2}\t.\tPASS\t.\tGT\t0/1\n')

print('VCF created')
"

# --- Stage 3: Liftover to GRCh38 ---

# Sort and compress the hg19 VCF
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bash -c "
    bcftools sort /genome/${SAMPLE}/raw/${SAMPLE}_hg19.vcf -Oz \
      -o /genome/${SAMPLE}/raw/${SAMPLE}_hg19_sorted.vcf.gz &&
    bcftools index -t /genome/${SAMPLE}/raw/${SAMPLE}_hg19_sorted.vcf.gz
  "

# Liftover using Picard
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  broadinstitute/picard:latest \
  java -jar /usr/picard/picard.jar LiftoverVcf \
    I=/genome/${SAMPLE}/raw/${SAMPLE}_hg19_sorted.vcf.gz \
    O=/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    CHAIN=/genome/liftover/hg19ToHg38.over.chain.gz \
    R=/genome/reference/Homo_sapiens_assembly38.fasta \
    REJECT=/genome/${SAMPLE}/raw/${SAMPLE}_liftover_rejected.vcf.gz \
    WARN_ON_MISSING_CONTIG=true

# Index the final VCF
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz"

echo "Done. VCF at: ${OUTDIR}/${SAMPLE}.vcf.gz"
```

**Expected output:** A bgzipped, tabix-indexed VCF on GRCh38 coordinates with 400,000-600,000 variants. Some variants (~5-15%) will be rejected during liftover because they map to regions that changed between builds.

### Alternative: Imputation First

For much better PRS and ancestry results, consider imputing your chip data before running the pipeline. Imputation fills in ~40 million positions from your 600K by using population reference panels.

1. Prepare chromosome-split VCFs from the hg19 data (step 14 in this pipeline)
2. Upload to the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/) or [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/)
3. Select the TOPMed or 1000G reference panel, GRCh38 output
4. Download the imputed VCF and use it as your pipeline input

Imputation typically takes a few hours on the server side. The imputed VCF will have ~40M variants with imputation quality scores (R2). Filter to R2 > 0.3 for downstream use.

> **Note:** Imputation servers may require a minimum number of samples or have usage policies. Single-sample imputation works on both servers listed above.

---

## Which Pipeline Steps Work with Chip Data

### Works Well

| Step | Name | Why It Works | Notes |
|---|---|---|---|
| **6** | ClinVar screen | Checks your variants against known pathogenic entries | You'll only find pathogenic variants that happen to be on the chip. Most clinically significant rare variants will be missed. |
| **7** | PharmCAT | Pharmacogenomic star alleles from SNP genotypes | This is one of the **best uses** of chip data. Most PGx-relevant positions are common SNPs that arrays cover well. CYP2D6 will still be limited. |
| **11** | ROH analysis | Runs of homozygosity from SNP genotypes | Works, but resolution is lower (~600K markers vs 5M). Large ROH (>5 MB) are still detectable. |
| **25** | PRS | Polygenic risk scores from common variants | Works **very well** — PRS scoring files are derived from GWAS arrays that overlap heavily with consumer chips. Variant matching rate may actually be higher than with WGS. |
| **26** | Ancestry PCA | Principal component analysis against 1000G | Works well for ancestry estimation. Common SNPs are what PCA uses. |
| **27** | CPIC lookup | Drug-gene recommendations | Works if step 7 (PharmCAT) succeeds. |

### Works with Limitations

| Step | Name | Limitation |
|---|---|---|
| **12** | Mito haplogroup | Only if your chip includes mitochondrial SNPs (23andMe does, ~3,000 mt positions). May assign a less specific haplogroup than WGS. |
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
| **16** | Coverage QC (indexcov) | No BAM to assess coverage |
| **20** | Mito variant calling (Mutect2) | Needs BAM |
| **21** | CYP2D6 (Cyrius) | Needs BAM |
| **22** | SV consensus merge | No SV calls |
| **23** | Clinical filter | Requires VEP-annotated VCF with gnomAD. Limited value on chip data. |

---

## What You Can Realistically Learn

### From chip data alone (no imputation)

| Analysis | Usefulness | Confidence |
|---|---|---|
| **Pharmacogenomics** | High | Good — most PGx SNPs are on the chip |
| **Ancestry / haplogroup** | High | Good for continental-level ancestry |
| **Carrier screening** | Low-moderate | Only finds carriers for variants on the chip |
| **Cancer predisposition** | Very low | Most pathogenic variants are rare and NOT on the chip |
| **Polygenic risk scores** | Moderate | Reasonable, but fewer matched variants than WGS or imputed data |

### From chip data + imputation

| Analysis | Usefulness | Confidence |
|---|---|---|
| **PRS** | High | Comparable to WGS for well-imputed regions |
| **Ancestry PCA** | High | Excellent with millions of imputed SNPs |
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
Convert to GRCh38 VCF (instructions above)
         |
         +---> Step 7:  PharmCAT -----> Step 27: CPIC drug-gene report
         |
         +---> Step 25: PRS (polygenic risk scores)
         |
         +---> Step 26: Ancestry PCA
         |
         +---> Step 6:  ClinVar screen (limited but useful)
         |
         +---> Step 11: ROH analysis
         |
    (optional)
         |
         v
Impute via Michigan/TOPMed server
         |
         +---> Re-run steps 25, 26 with imputed VCF (better PRS, better PCA)
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
| Chip + imputation | $0 (free servers) | ~40M (imputed) | Better PRS, ancestry |
| Budget WGS (Novogene) | $200-400 | ~5M (observed) | Full pipeline |
| Standard WGS (Nebula, Dante) | $300-600 | ~5M (observed) | Full pipeline + data retention |

If you are serious about genomic analysis, WGS is worth the investment. But if you already have chip data and want to start exploring, the pharmacogenomics and PRS steps deliver real value today.
