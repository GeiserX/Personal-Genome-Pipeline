# Long-Read Sequencing Guide

## What Is Long-Read Sequencing?

Standard WGS (from Illumina, Novogene, Dante Labs, etc.) produces **short reads** — millions of DNA fragments typically 150 base pairs long. Two companies offer a fundamentally different approach: reading much longer stretches of DNA in a single pass.

### Oxford Nanopore Technology (ONT)

ONT sequencers (MinION, PromethION, GridION) push a strand of DNA through a tiny protein pore. As each nucleotide passes through, it changes the electrical current. The device reads these current changes in real time to determine the sequence.

- **Read lengths:** 1,000 to 100,000+ bp (some reads exceed 1 million bp)
- **Accuracy:** ~99% with R10.4.1 chemistry and super-accuracy basecalling (Q20+)
- **Common data format:** FASTQ (basecalled from FAST5/POD5 raw signal)
- **Vendors/services:** ONT direct (MinION, PromethION), some clinical labs

### PacBio HiFi (High-Fidelity)

PacBio sequencers (Sequel II, Revio) use a different approach: a polymerase enzyme copies a circular DNA template while a camera records fluorescent signals from each incorporated nucleotide. By reading the same circle multiple times, errors cancel out.

- **Read lengths:** 10,000-25,000 bp
- **Accuracy:** >99.9% (Q30+, comparable to Illumina)
- **Common data format:** Unaligned BAM or FASTQ (HiFi/CCS reads)
- **Vendors/services:** PacBio direct, Revio service providers, some clinical labs

### Why Use Long Reads?

Long reads solve problems that short reads cannot:

| Capability | Short Reads (150bp) | Long Reads (10-100Kbp) |
|---|---|---|
| SNPs and small indels | Excellent | Excellent (HiFi) / Good (ONT) |
| Structural variants (>50bp) | Limited | Excellent |
| Repeat expansions | Indirect (ExpansionHunter) | Direct measurement |
| Phasing (which allele is on which chromosome) | Statistical | Read-based, much more accurate |
| Segmental duplications / tandem repeats | Mostly invisible | Resolvable |
| Methylation | Requires bisulfite conversion | Native (ONT) |

---

## Pipeline Compatibility

Not every step in this pipeline works with long-read data. Here is the full breakdown.

### Works As-Is (No Changes Needed)

These steps take a VCF or BAM and work identically regardless of read technology:

| Step | Tool | Notes |
|---|---|---|
| 5 | AnnotSV | Takes any SV VCF — point it at Sniffles2 output |
| 6 | ClinVar Screen | VCF-only, technology-independent |
| 7 | PharmCAT | VCF-only, technology-independent |
| 11 | ROH Analysis | VCF-only via plink2 |
| 12 | Mito Haplogroup | VCF-only via haplogrep3 |
| 13 | VEP Annotation | VCF-only |
| 14 | Imputation Prep | VCF-only |
| 17 | CPSR | VCF-only |
| 22 | SURVIVOR Merge | Takes any SV VCF set |
| 23 | Clinical Filter | VCF-only |
| 24 | HTML Report | Aggregates existing outputs |
| 25 | PRS | VCF-only via plink2 |
| 26 | Ancestry | VCF-only via plink2 |
| 27 | CPIC Lookup | Reads PharmCAT JSON output |

### Needs Long-Read Specific Scripts

| Step | Short-Read Tool | Long-Read Alternative | Script |
|---|---|---|---|
| 2 | minimap2 (`-x sr`) | minimap2 (`-x map-ont` / `-x map-hifi`) | `scripts/02b-alignment-longread.sh` |
| 3 | DeepVariant (WGS model) | **Clair3** or DeepVariant (ONT/PACBIO model) | `scripts/03e-clair3.sh` |
| 4 | Manta | **Sniffles2** | `scripts/04c-sniffles2.sh` |

### Does NOT Work with Long Reads

These tools are specifically designed for short-read data and will produce incorrect results or crash with long-read input:

| Step | Tool | Why It Fails | Alternative |
|---|---|---|---|
| 1b | fastp | Expects paired-end Illumina reads; no long-read mode | NanoPlot, LongQC, or chopper for long-read QC |
| 2a | BWA-MEM2 | Short-read aligner only; cannot handle reads >500bp | Use minimap2 with long-read preset |
| 3a | GATK HaplotypeCaller | Active-region model assumes short reads; poor with long reads | Clair3 or DeepVariant long-read mode |
| 3b | FreeBayes | Bayesian model tuned for short reads | Clair3 or DeepVariant long-read mode |
| 3c | Strelka2 | Scoring model trained on BWA-MEM short-read data | Clair3 or DeepVariant long-read mode |
| 3d | Octopus | Haplotype model designed for short reads | Clair3 or DeepVariant long-read mode |
| 4 | Manta | Illumina-specific insert size model | Sniffles2 |
| 4a | TIDDIT | Short-read coverage/discordance model | Sniffles2 |
| 4b | GRIDSS | Assembly-based, short-read specific | Sniffles2 |
| 9 | ExpansionHunter | Illumina short-read graph model; expects paired-end data | TRGT (PacBio), STRique (ONT), or direct long-read spanning |
| 15 | duphold | Re-genotypes SVs using short-read depth models | Not needed — Sniffles2 QUAL scores are reliable |
| 18 | CNVnator | Read-depth model calibrated for short reads | Sniffles2 detects CNVs natively |
| 19 | Delly | Paired-end and split-read model | Sniffles2 |
| 21 | Cyrius (CYP2D6) | Short-read depth-based star allele caller | Paraphase (long-read CYP2D6 resolver) |

---

## Step-by-Step: Long-Read Analysis

### Prerequisites

Same reference data as the short-read pipeline:

```bash
export GENOME_DIR=/path/to/your/data
# Reference genome and index should already exist from initial setup
# See docs/00-reference-setup.md
```

### Step 1: Prepare Your Data

Place your long-read data in the sample directory:

```bash
# Create sample directory
mkdir -p ${GENOME_DIR}/your_sample/fastq/

# Copy or symlink your long-read FASTQ
# ONT: typically one FASTQ per run (basecalled from FAST5/POD5)
cp /path/to/ont_reads.fastq.gz ${GENOME_DIR}/your_sample/fastq/your_sample.fastq.gz

# PacBio: typically one FASTQ or unaligned BAM per cell
cp /path/to/hifi_reads.fastq.gz ${GENOME_DIR}/your_sample/fastq/your_sample.fastq.gz
# OR for unaligned BAM:
cp /path/to/hifi_reads.bam ${GENOME_DIR}/your_sample/fastq/your_sample.bam
```

Note: long-read data is **single-end** (no R1/R2 pairs). The alignment script expects a single file, not paired files.

### Step 2: Align Reads

```bash
# For Oxford Nanopore
PLATFORM=ont ./scripts/02b-alignment-longread.sh your_sample

# For PacBio HiFi
PLATFORM=hifi ./scripts/02b-alignment-longread.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/aligned_longread/your_sample_sorted.bam`

The script auto-detects input files in this order:
1. `fastq/your_sample.fastq.gz`
2. `fastq/your_sample_lr.fastq.gz`
3. `fastq/your_sample.fq.gz`
4. `fastq/your_sample.bam` (unaligned BAM)

Override with: `INPUT=/path/to/reads.fastq.gz PLATFORM=ont ./scripts/02b-alignment-longread.sh your_sample`

### Step 3: Variant Calling (SNPs + Indels)

You have two options:

#### Option A: Clair3 (recommended for long reads)

```bash
# For ONT
PLATFORM=ont ./scripts/03e-clair3.sh your_sample

# For PacBio HiFi
PLATFORM=hifi ./scripts/03e-clair3.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/vcf_clair3/your_sample.vcf.gz`

Clair3 uses deep-learning models trained specifically on long-read error profiles. It auto-selects the correct model based on your PLATFORM setting.

#### Option B: DeepVariant (long-read mode)

DeepVariant 1.6.0+ also supports long reads with dedicated models:

```bash
SAMPLE=your_sample

docker run --rm \
  --cpus 8 --memory 32g \
  -v "${GENOME_DIR}:/genome" \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=ONT_R104 \
    --ref="/genome/reference/Homo_sapiens_assembly38.fasta" \
    --reads="/genome/${SAMPLE}/aligned_longread/${SAMPLE}_sorted.bam" \
    --output_vcf="/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    --num_shards=8

# For PacBio, use --model_type=PACBIO instead of ONT_R104
```

**Clair3 vs DeepVariant for long reads:**

| Aspect | Clair3 | DeepVariant |
|---|---|---|
| ONT accuracy | Excellent (purpose-built) | Good (newer models improving) |
| PacBio HiFi accuracy | Excellent | Excellent |
| Speed | Generally faster | Slower (more compute-heavy) |
| Indel calling | Very good | Very good |
| Community support | Active (HKU) | Active (Google) |
| Model updates | Frequent (matches new chemistries) | Periodic |

For PacBio HiFi data, both tools perform comparably. For ONT data, Clair3 tends to have an edge because its models are updated more frequently to match new ONT basecalling improvements.

### Step 4: Structural Variant Calling

```bash
# Default: uses aligned_longread/ directory
./scripts/04c-sniffles2.sh your_sample

# Or explicitly:
ALIGN_DIR=aligned_longread ./scripts/04c-sniffles2.sh your_sample
```

Output: `${GENOME_DIR}/your_sample/sv_sniffles/your_sample_sv.vcf.gz`

Sniffles2 excels at detecting all SV types from long reads. Typical output for a 30X long-read genome:
- 15,000-25,000 SVs (vs 7,000-9,000 from short-read Manta)
- Much better at insertions (short reads systematically miss these)
- Better breakpoint resolution (often single-nucleotide)
- Reliable in repetitive regions where short reads struggle

### Step 5: Continue with Standard Pipeline Steps

After alignment and variant calling, most downstream steps work normally. Point them at the long-read outputs:

```bash
SAMPLE=your_sample

# ClinVar screen (uses VCF from Clair3)
VCF_DIR=vcf_clair3 ./scripts/06-clinvar-screen.sh "$SAMPLE"

# PharmCAT (uses VCF from Clair3)
# Note: PharmCAT needs the VCF at the default location, so either:
# 1. Symlink: ln -s vcf_clair3/${SAMPLE}.vcf.gz vcf/${SAMPLE}.vcf.gz
# 2. Or copy the Clair3 VCF to the expected vcf/ directory
./scripts/07-pharmacogenomics.sh "$SAMPLE"

# VEP annotation
./scripts/13-vep-annotation.sh "$SAMPLE"

# Annotate Sniffles2 SVs
SV_VCF="${GENOME_DIR}/${SAMPLE}/sv_sniffles/${SAMPLE}_sv.vcf.gz" ./scripts/05-annotsv.sh "$SAMPLE"
```

---

## Docker Images Used

| Tool | Image | Size |
|---|---|---|
| minimap2 | `quay.io/biocontainers/minimap2:2.28--he4a0461_0` | ~30 MB |
| samtools | `staphb/samtools:1.20` | ~200 MB |
| Clair3 | `hkubal/clair3:v2.0.0` | ~3 GB (includes all models) |
| Sniffles2 | `quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0` | ~200 MB |
| DeepVariant | `google/deepvariant:1.6.0` | ~5 GB |

Pre-pull images before your first run:

```bash
docker pull quay.io/biocontainers/minimap2:2.28--he4a0461_0
docker pull staphb/samtools:1.20
docker pull hkubal/clair3:v2.0.0
docker pull quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0
```

---

## ONT-Specific Notes

### Basecalling Matters

ONT raw data (FAST5/POD5) must be basecalled before entering this pipeline. Use the latest Dorado basecaller with the **super accuracy (SUP)** model for best variant calling:

```bash
# This happens BEFORE the pipeline — on a GPU machine
dorado basecaller sup pod5_dir/ > reads.bam
dorado demux reads.bam  # if barcoded
```

The Clair3 model (`r1041_e82_400bps_sup_v500`) assumes R10.4.1 chemistry with SUP basecalling. Using a different chemistry or basecaller may reduce accuracy.

### Read Length Distribution

ONT generates a wide range of read lengths. Very short reads (<1kb) add noise. If your N50 is below 5kb, consider filtering:

```bash
# Filter reads shorter than 1kb (optional, before alignment)
# Use chopper or NanoFilt
```

### Methylation

ONT natively detects methylation (5mC, 6mA) during basecalling with Dorado. Methylation tags are stored in the BAM as MM/ML tags. This pipeline does not currently process methylation data, but the aligned BAM preserves these tags for future use.

---

## PacBio-Specific Notes

### HiFi vs CLR

This pipeline supports **HiFi (CCS)** reads only, not older CLR (continuous long reads). HiFi reads are generated by PacBio Sequel II, Sequel IIe, and Revio systems using circular consensus sequencing.

If you received CLR data, you'll need to generate CCS reads first using `ccs` or `pbccs` (PacBio tool, outside this pipeline).

### Kinetic Tags

PacBio HiFi BAMs may contain kinetic modification tags (fi, fp, ri, rp) for methylation detection. Like ONT methylation, these are preserved in alignment but not processed by this pipeline.

---

## Hybrid Analysis (Short + Long Reads)

If you have both short-read and long-read data for the same sample, you can run both pipelines and combine results:

1. **Short-read pipeline:** standard `02-alignment.sh` + `03-deepvariant.sh` + `04-manta.sh`
2. **Long-read pipeline:** `02b-alignment-longread.sh` + `03e-clair3.sh` + `04c-sniffles2.sh`
3. **Merge variant calls:**
   - SNPs/indels: use `benchmark-variants.sh` to compare, then take the intersection for high-confidence calls
   - SVs: use `22-survivor-merge.sh` with both Manta and Sniffles2 VCFs for consensus SVs

This is the gold standard approach when both data types are available.

---

## Troubleshooting

### "minimap2 runs out of memory"

Long-read alignment uses more memory than short-read because the index is loaded into RAM. Increase Docker memory limit:

```bash
# In the script, change --memory 16g to --memory 24g or higher
```

### "Clair3 model not found"

The Docker image `hkubal/clair3:v2.0.0` bundles models at `/opt/models/`. If you get a model-not-found error:

1. Verify the image is pulled: `docker images | grep clair3`
2. Check available models: `docker run --rm hkubal/clair3:v2.0.0 ls /opt/models/`
3. If your ONT chemistry is different (e.g., R9.4.1), look for a matching model in the container

### "Sniffles2 produces 0 SVs"

- Verify the BAM is from long-read alignment (not short-read): `samtools stats your.bam | grep "average length"`
- Average read length should be >1000bp for ONT or >5000bp for HiFi
- Check that the BAM is sorted and indexed

### "DeepVariant long-read model fails"

- Use `--model_type=ONT_R104` for ONT (not `ONT` alone)
- Use `--model_type=PACBIO` for PacBio HiFi
- DeepVariant 1.6.0 requires the exact model type string — check the DeepVariant release notes for your version

### "VCF downstream steps fail on Clair3 output"

Some downstream steps expect the VCF at the default `vcf/` path. Either:
1. Symlink: `ln -s ../vcf_clair3/${SAMPLE}.vcf.gz ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz`
2. Or set `VCF_DIR=vcf_clair3` if the step supports it

---

## Expected Output Comparison

Rough numbers for a 30X human genome:

| Metric | Short-Read (Illumina) | ONT (R10.4.1 SUP) | PacBio HiFi (Revio) |
|---|---|---|---|
| SNPs | ~4.5M | ~4.3M | ~4.5M |
| Indels | ~1.0M | ~0.8M | ~1.0M |
| PASS rate | ~93% | ~85-90% | ~93% |
| SVs (total) | ~7-9K (Manta) | ~20-25K (Sniffles2) | ~18-22K (Sniffles2) |
| SV insertions | ~500 (poor) | ~8-10K (excellent) | ~7-9K (excellent) |
| Alignment time | 1-2 hours | 1-3 hours | 1-2 hours |
| Variant calling time | 2-4 hours (DV) | 2-4 hours (Clair3) | 2-4 hours (Clair3) |

The biggest advantage of long reads is SV detection — especially insertions, which short reads systematically miss.
