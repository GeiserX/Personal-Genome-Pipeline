# Step 29: Somatic Variant Calling (Mutect2 Tumor-Only) [EXPERIMENTAL]

## What This Does

Detects **somatic mutations** -- variants acquired during your lifetime, not inherited from your parents. These include:

- **Clonal hematopoiesis of indeterminate potential (CHIP)**: Age-related mutations in blood cells (common after age 40, clinically relevant for cardiovascular risk and blood cancers)
- **Mosaic variants**: Mutations present in only a fraction of cells, arising from post-zygotic events
- **Other acquired mutations**: Environmental damage, replication errors, etc.

This is fundamentally different from the germline variant calling in step 3 (DeepVariant), which finds variants you were born with. Somatic variants have **low allele fractions** (often 1-30%) because they only exist in a subset of cells.

## Why [EXPERIMENTAL]

This step runs Mutect2 in **tumor-only mode** -- meaning there is no matched normal sample for comparison. In a clinical setting, somatic variant calling uses a tumor sample AND a matched normal (blood) to distinguish somatic mutations from germline variants. Without a matched normal:

- **High false positive rate**: Many germline variants (especially rare ones not in gnomAD) will be called as somatic
- **Limited sensitivity for low-AF variants**: Without a normal sample to establish the baseline, true low-frequency somatic events are harder to distinguish from noise
- **Germline contamination**: Common germline SNPs that happen to be slightly off 0.5/1.0 allele fraction due to sequencing noise can appear "somatic"

The gnomAD resource and Panel of Normals help reduce false positives significantly, but cannot eliminate them entirely.

**Use this step for exploratory analysis only. Do not make medical decisions based on tumor-only somatic calls without clinical validation.**

## When Is This Useful?

- **CHIP screening**: Looking for age-related clonal hematopoiesis variants (DNMT3A, TET2, ASXL1, TP53, etc.) from blood-derived WGS
- **Mosaicism**: Detecting mosaic variants that germline callers miss because they expect 50%/100% allele fractions
- **Research**: Exploring somatic mutation burden, mutational signatures, or clonal dynamics
- **Complement to germline**: Some pathogenic variants in cancer genes may be somatic rather than germline

## Tool

- **GATK Mutect2** (Broad Institute) in tumor-only mode (no `-I normal` or `-normal` flags)

## Docker Image

```
broadinstitute/gatk:4.6.1.0
```

Already used in step 20 (mitochondrial analysis). No additional download needed.

## Prerequisites

### Required

- Sorted BAM with index (from step 2)
- GRCh38 reference FASTA with `.fai` and `.dict` (see [00-reference-setup.md](00-reference-setup.md))

### Recommended Resources (Reduce False Positives)

These resources are optional but **strongly recommended**. Without them, the output will contain far more false positive somatic calls.

#### gnomAD AF-Only VCF (~6.5 GB)

Contains population allele frequencies from gnomAD. Mutect2 uses this to flag likely germline variants (high AF in the population = probably not somatic).

```bash
mkdir -p ${GENOME_DIR}/somatic

# Download via gsutil (if installed)
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
  ${GENOME_DIR}/somatic/
gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
  ${GENOME_DIR}/somatic/

# Alternative: direct HTTPS download
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/af-only-gnomad.hg38.vcf.gz.tbi
```

#### Panel of Normals (PoN, ~1 GB)

Built from 1000 Genomes data. Contains technical artifacts and recurrent sequencing errors seen across many normal samples. Mutect2 uses this to filter out variants that are likely noise rather than real somatic events.

```bash
# Download via gsutil (if installed)
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
  ${GENOME_DIR}/somatic/
gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
  ${GENOME_DIR}/somatic/

# Alternative: direct HTTPS download
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz
wget -c https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi \
  -O ${GENOME_DIR}/somatic/1000g_pon.hg38.vcf.gz.tbi
```

> **Note:** `gsutil` is part of the Google Cloud SDK. If you do not have it installed, use the `wget` alternative URLs above (same files, just accessed over HTTPS instead of the gs:// protocol).

## Command

```bash
export GENOME_DIR=/path/to/your/data
./scripts/29-mutect2-somatic.sh your_name
```

### Environment Variables

| Variable | Default | Description |
|---|---|---|
| `GENOME_DIR` | (required) | Path to your data directory |
| `THREADS` | 4 | CPU threads for Mutect2 |
| `INTERVALS` | (empty = full genome) | Restrict to a region, e.g. `chr22` or `chr17:7500000-7700000` (TP53 locus) |
| `ALIGN_DIR` | `aligned` | Use `aligned_bwamem2` for BWA-MEM2 alignments |

### Quick Test on a Single Chromosome

Running on the full genome takes 2-6 hours. To test quickly:

```bash
INTERVALS=chr22 ./scripts/29-mutect2-somatic.sh your_name
# ~15-30 minutes
```

## Output

All files are written to `${GENOME_DIR}/${SAMPLE}/somatic/`:

| File | Description |
|---|---|
| `${SAMPLE}_somatic_unfiltered.vcf.gz` | Raw Mutect2 calls (before filtering) |
| `${SAMPLE}_somatic_unfiltered.vcf.gz.stats` | Mutect2 internal statistics (used by FilterMutectCalls) |
| `${SAMPLE}_somatic_filtered.vcf.gz` | Filtered calls with PASS/FAIL annotations |
| `${SAMPLE}_somatic_filtered.vcf.gz.tbi` | Tabix index for the filtered VCF |

## Interpreting Results

### What the FILTER Field Means

After FilterMutectCalls, each variant gets a FILTER status:

| Filter | Meaning |
|---|---|
| `PASS` | Passed all filters -- candidate somatic variant |
| `germline` | Likely germline (high gnomAD AF or high allele fraction) |
| `normal_artifact` | Matches a Panel of Normals entry (technical artifact) |
| `weak_evidence` | Low quality scores / insufficient reads supporting the variant |
| `strand_bias` | Variant reads come overwhelmingly from one strand (artifact signal) |
| `contamination` | Possible sample contamination |
| `orientation` | Orientation bias artifact (common in FFPE samples, rare in blood WGS) |

### Allele Fraction (AF) Interpretation

In tumor-only mode from blood WGS:

| AF Range | Likely Source |
|---|---|
| 0.45-0.55 | Heterozygous germline (false positive) |
| ~1.0 | Homozygous germline (false positive) |
| 0.01-0.10 | Possible low-frequency somatic (CHIP candidate) |
| 0.10-0.40 | Could be somatic, mosaic, or germline with noise |

**Key insight**: In a healthy individual's blood WGS, the vast majority of PASS calls will be germline variants that escaped filtering. True somatic variants are rare events -- a healthy 40-year-old might have 0-20 genuine CHIP mutations detectable at 30X coverage.

### Finding CHIP Candidates

Clonal hematopoiesis variants are found in specific genes. After running, check for PASS variants in known CHIP genes:

```bash
# Extract PASS variants and look for known CHIP genes using VEP output (step 13)
bcftools view -f PASS ${GENOME_DIR}/${SAMPLE}/somatic/${SAMPLE}_somatic_filtered.vcf.gz \
  | grep -E "DNMT3A|TET2|ASXL1|TP53|JAK2|SF3B1|SRSF2|PPM1D|CBL|GNB1|IDH1|IDH2"
```

### Cross-Referencing with Other Steps

| Step | How It Helps |
|---|---|
| Step 3 (DeepVariant) | Compare: if DeepVariant calls it at ~50% AF, it is almost certainly germline |
| Step 6 (ClinVar) | Check if the somatic variant is in ClinVar as pathogenic |
| Step 13 (VEP) | Annotate somatic calls with functional impact and gnomAD frequencies |
| Step 17 (CPSR) | Cancer predisposition report covers germline cancer gene variants |

## Runtime

| Scope | Approximate Time | Memory |
|---|---|---|
| Full genome (no intervals) | 2-6 hours | 8 GB |
| Single chromosome (`INTERVALS=chr22`) | 15-30 minutes | 8 GB |
| Targeted region (e.g., TP53 locus) | <5 minutes | 8 GB |

## Tumor-Only vs. Tumor-Normal: Why It Matters

| Aspect | Tumor-Only (this script) | Tumor-Normal (clinical) |
|---|---|---|
| **Input** | One BAM (blood/tissue) | Two BAMs (tumor + matched normal) |
| **Germline filtering** | gnomAD AF + PoN only | Direct subtraction of normal genotype |
| **False positive rate** | High (thousands of calls in healthy tissue) | Low (tens to hundreds of true somatic calls) |
| **Rare germline variants** | Often called as somatic | Correctly filtered out |
| **Use case** | CHIP screening, mosaicism research | Cancer genomics, treatment selection |
| **Clinical validity** | Exploratory only | Clinically actionable |

If you have access to a matched normal sample (e.g., blood for a solid tumor), you can modify the script to add the normal BAM:

```bash
# This is NOT implemented in the script -- manual modification needed
gatk Mutect2 \
  -R reference.fasta \
  -I tumor.bam \
  -I normal.bam \
  -normal normal_sample_name \
  --germline-resource af-only-gnomad.hg38.vcf.gz \
  -pon 1000g_pon.hg38.vcf.gz \
  -O somatic.vcf.gz
```

## Notes

- The script is **idempotent**: if the filtered output VCF already exists, it skips execution. Delete the output file to force re-run.
- `--max-mnp-distance 0` prevents merging adjacent SNPs into multi-nucleotide polymorphisms (consistent with step 20).
- The `.stats` file generated by Mutect2 is automatically consumed by FilterMutectCalls. Do not delete it before filtering completes.
- The `ALIGN_DIR` variable lets you use BWA-MEM2 alignments (`ALIGN_DIR=aligned_bwamem2`) if available.
- GATK Docker image is ~2.2 GB but is shared with step 20 (mitochondrial analysis) -- no extra download if you already have it.
- For WGS data from consumer vendors (Nebula, Dante, etc.), the sequencing is from blood, so CHIP is the primary type of somatic variant you can detect. Tissue-specific somatic variants require sequencing the relevant tissue.
