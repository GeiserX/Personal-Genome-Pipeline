# Variant Caller Benchmarking

## Why Benchmark Callers?

No single variant caller is universally best. Each tool uses a different algorithm -- deep learning (DeepVariant), local haplotype assembly (GATK HaplotypeCaller), Bayesian haplotype-based calling (FreeBayes), or heuristic scoring (Strelka2) -- and each makes different tradeoffs between precision and recall. A variant found by one caller may be missed by another, and vice versa.

Benchmarking lets you:

1. **Quantify accuracy** against a known truth set (GIAB) -- precision, recall, and F1 for your data and hardware
2. **Measure concordance** between callers -- understand how much they agree and where they diverge
3. **Choose the right caller** for your use case -- highest accuracy for a clinical question, broadest sensitivity for research exploration, or fastest runtime for iterative analyses
4. **Build confidence** in your results -- variants called by two or more independent tools are far more likely to be real

For background on caller performance differences, see:
- [PLOS ONE variant caller comparison (2024)](https://doi.org/10.1371/journal.pone.0339891) -- systematic comparison of DeepVariant, GATK, and FreeBayes on WGS data
- [UMCCR BWA vs minimap2 (2021)](https://umccr.org/blog/bwa-mem-vs-minimap2/) -- alignment-level differences and their downstream effects on variant calling

---

## Available Alternative Caller Scripts

The pipeline ships with four variant callers. Each writes output to a separate directory so they can run in parallel without conflicts.

| Script | Caller | Output Directory | Description |
|---|---|---|---|
| `scripts/03-deepvariant.sh` | DeepVariant 1.6.0 | `vcf/` | Default caller. Deep learning model, highest precision and F1. |
| `scripts/03a-gatk-haplotypecaller.sh` | GATK HaplotypeCaller 4.6.1 | `vcf_gatk/` | Gold standard in clinical labs. Good precision/recall balance, GVCF support. |
| `scripts/03b-freebayes.sh` | FreeBayes 1.3.6 | `vcf_freebayes/` | Bayesian caller. Highest sensitivity, most false positives, single-threaded. |
| `scripts/03c-strelka2-germline.sh` | Strelka2 2.9.10 | `vcf_strelka2/` | Fast heuristic caller (SNVs + indels). Best with BWA-MEM2 alignments. |

All four scripts accept the same arguments:

```bash
export GENOME_DIR=/path/to/your/data
./scripts/03-deepvariant.sh your_sample
./scripts/03a-gatk-haplotypecaller.sh your_sample
./scripts/03b-freebayes.sh your_sample
./scripts/03c-strelka2-germline.sh your_sample
```

GATK and FreeBayes also accept an `INTERVALS` environment variable to restrict calling to a specific region (see the chr22 example below). Strelka2 and TIDDIT always process the full genome.

An alternative aligner is also available:

| Script | Aligner | Output Directory | Description |
|---|---|---|---|
| `scripts/02-alignment.sh` | minimap2 2.28 | `aligned/` | Default. Faster, good for germline WGS. |
| `scripts/02a-alignment-bwamem2.sh` | BWA-MEM2 2.2.1 | `aligned_bwamem2/` | Produces XS tags needed by Strelka2. Slightly more accurate for somatic calling. |

---

## Benchmarking Modes

There are two ways to evaluate caller performance:

### Mode 1: Pairwise Concordance (No Truth Set Required)

Compare two or more caller VCFs against each other to measure agreement. This tells you how much the callers agree, but not which one is "right."

Use `bcftools isec` to compute the intersection and unique calls:

```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# Compare DeepVariant vs GATK (PASS variants only)
mkdir -p "${GENOME_DIR}/${SAMPLE}/benchmark"

docker run --rm -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
  bcftools isec -p /genome/${SAMPLE}/benchmark/dv_vs_gatk \
    -f PASS \
    /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    /genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz
```

This produces four files in the output directory:

| File | Contents |
|---|---|
| `0000.vcf` | Variants unique to caller A (DeepVariant) |
| `0001.vcf` | Variants unique to caller B (GATK) |
| `0002.vcf` | Shared variants (from caller A's perspective) |
| `0003.vcf` | Shared variants (from caller B's perspective) |
| `sites.txt` | Per-site presence/absence table |

Compute Jaccard similarity (shared / union):

```bash
SHARED=$(grep -c -v '^#' "${GENOME_DIR}/${SAMPLE}/benchmark/dv_vs_gatk/0002.vcf")
UNIQUE_A=$(grep -c -v '^#' "${GENOME_DIR}/${SAMPLE}/benchmark/dv_vs_gatk/0000.vcf")
UNIQUE_B=$(grep -c -v '^#' "${GENOME_DIR}/${SAMPLE}/benchmark/dv_vs_gatk/0001.vcf")
UNION=$((SHARED + UNIQUE_A + UNIQUE_B))
echo "Jaccard: $(echo "scale=4; $SHARED / $UNION" | bc)"
```

### Mode 2: Truth Set Benchmarking (GIAB)

Compare caller output against a known-correct truth set. This gives you absolute precision, recall, and F1 scores.

**Important:** Truth set benchmarking is only valid when the query VCF was generated from the **same biological sample** as the truth set. If the truth set is HG002, you must sequence HG002 DNA, align it, call variants, and then compare those calls against the HG002 truth VCF. Running hap.py with your personal sample against the HG002 truth set produces meaningless precision/recall numbers because the "true" variants are different for every individual.

The Genome in a Bottle (GIAB) consortium publishes validated truth sets for several reference samples. The most widely used is **HG002** (Ashkenazi Jewish male), which has the most comprehensive high-confidence calls.

---

## GIAB Truth Set Setup

### Download HG002 Truth Set

```bash
GENOME_DIR=/path/to/your/data
mkdir -p "${GENOME_DIR}/giab"

# HG002 truth VCF (GRCh38, v4.2.1)
wget -c -P "${GENOME_DIR}/giab" \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

wget -c -P "${GENOME_DIR}/giab" \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"

# HG002 high-confidence regions BED
wget -c -P "${GENOME_DIR}/giab" \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
```

The BED file defines the regions where the truth set is confident. Variants outside these regions are excluded from benchmarking because the truth status is unknown.

### Alternative: NA12878 (HG001)

The pipeline's [quick-test.md](quick-test.md) uses NA12878 (HG001), which is also a valid truth set:

```bash
wget -c -P "${GENOME_DIR}/giab" \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
```

HG002 is preferred for benchmarking because its truth set covers more difficult genomic regions.

### Directory Structure After Setup

```
${GENOME_DIR}/giab/
  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz           # Truth VCF
  HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi        # Truth VCF index
  HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed # High-confidence BED
```

---

## Running Truth Set Benchmarking with hap.py

[hap.py](https://github.com/Illumina/hap.py) (Illumina) is the standard benchmarking tool for SNP and indel callers. It decomposes complex variants, performs genotype matching, and reports precision/recall/F1 stratified by variant type.

```bash
# IMPORTANT: SAMPLE must be the GIAB sample that matches the truth set.
# If using HG002 truth, you must have sequenced and called variants on HG002.
SAMPLE=HG002
GENOME_DIR=/path/to/your/data
TRUTH_VCF="/genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONF_BED="/genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
REF="/genome/reference/Homo_sapiens_assembly38.fasta"

# Benchmark DeepVariant against truth set
docker run --rm \
  --cpus 4 --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    "$TRUTH_VCF" \
    "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    -r "$REF" \
    -f "$CONF_BED" \
    -o "/genome/${SAMPLE}/benchmark/deepvariant_vs_truth" \
    --threads 4

# Repeat for GATK
docker run --rm \
  --cpus 4 --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    "$TRUTH_VCF" \
    "/genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz" \
    -r "$REF" \
    -f "$CONF_BED" \
    -o "/genome/${SAMPLE}/benchmark/gatk_vs_truth" \
    --threads 4

# Repeat for FreeBayes
docker run --rm \
  --cpus 4 --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    "$TRUTH_VCF" \
    "/genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz" \
    -r "$REF" \
    -f "$CONF_BED" \
    -o "/genome/${SAMPLE}/benchmark/freebayes_vs_truth" \
    --threads 4
```

Each hap.py run produces several output files:

| File | Contents |
|---|---|
| `*.summary.csv` | Precision, recall, F1 by variant type (SNP, INDEL) |
| `*.extended.csv` | Detailed metrics by quality tier and subtype |
| `*.vcf.gz` | Annotated VCF with TP/FP/FN labels per variant |
| `*.roc.*` | ROC curve data for quality-score-based filtering |

---

## Interpreting Results

### Key Metrics

| Metric | What It Measures | Formula |
|---|---|---|
| **Precision** | How many called variants are real (low = too many false positives) | TP / (TP + FP) |
| **Recall** | How many real variants were found (low = missing real variants) | TP / (TP + FN) |
| **F1** | Harmonic mean of precision and recall (balanced accuracy) | 2 * (Precision * Recall) / (Precision + Recall) |
| **Jaccard** | Overlap between two call sets (pairwise mode) | Shared / (Shared + Unique_A + Unique_B) |

Where:
- **TP** (True Positive) = variant called correctly
- **FP** (False Positive) = variant called but does not exist in truth
- **FN** (False Negative) = variant exists in truth but was not called

### What the Numbers Mean

- **F1 > 0.99**: Excellent. State-of-the-art for SNPs on 30X WGS.
- **F1 0.95-0.99**: Good. Typical for indels, or SNPs with a less optimal caller.
- **F1 < 0.95**: Investigate. May indicate alignment issues, low coverage, or caller misconfiguration.
- **Jaccard > 0.95**: Two callers agree on the vast majority of calls.
- **Jaccard 0.85-0.95**: Moderate agreement. Differences are mostly in low-confidence regions or indels.
- **Jaccard < 0.85**: Significant disagreement. Worth investigating the unique calls from each caller.

### Typical Results (30X WGS, GRCh38)

These are representative numbers from published benchmarks and real pipeline runs. Your results will vary depending on coverage, sample quality, and alignment.

| Caller | SNP Precision | SNP Recall | SNP F1 | Indel Precision | Indel Recall | Indel F1 |
|---|---|---|---|---|---|---|
| DeepVariant 1.6 | 0.999+ | 0.998+ | 0.999 | 0.995+ | 0.993+ | 0.994 |
| GATK HC 4.6 | 0.998 | 0.997 | 0.998 | 0.985 | 0.980 | 0.983 |
| FreeBayes 1.3 | 0.990 | 0.998 | 0.994 | 0.950 | 0.970 | 0.960 |

Pairwise concordance (Jaccard similarity, PASS SNPs only):

| Pair | Typical Jaccard |
|---|---|
| DeepVariant vs GATK | 0.95 -- 0.98 |
| DeepVariant vs FreeBayes | 0.88 -- 0.93 |
| GATK vs FreeBayes | 0.87 -- 0.92 |

Key observations:
- DeepVariant and GATK agree on ~95-98% of SNPs. Most disagreements are in low-complexity regions and around indels.
- FreeBayes calls more variants than either (higher sensitivity), but many of the extra calls are false positives, pulling its precision down.
- Indel calling is harder than SNP calling for all tools. The gap between callers is wider for indels.
- FreeBayes is the most useful as a "second opinion" caller -- its unique calls include both real variants the other callers missed and a larger set of false positives.

---

## Practical Example: chr22-Only Benchmark

Running all four callers on the full genome takes 1-12 hours each (FreeBayes being the slowest at 8-12h single-threaded, Strelka2 the fastest at ~1h). For a quick script validation, you can restrict GATK and FreeBayes to chromosome 22 only (~1-2% of the genome, takes 3-8 minutes per caller), but **full-genome runs are strongly recommended for meaningful benchmarking**.

### Step 1: Run Callers on chr22 (GATK + FreeBayes)

```bash
export GENOME_DIR=/path/to/your/data
export SAMPLE=your_sample

# DeepVariant does not support INTERVALS — use the full VCF and filter afterward
# (or run on the full genome; it is already fast enough on modern hardware)

# GATK on chr22 only
INTERVALS=chr22 ./scripts/03a-gatk-haplotypecaller.sh "$SAMPLE"

# FreeBayes on chr22 only
INTERVALS=chr22 ./scripts/03b-freebayes.sh "$SAMPLE"
```

For DeepVariant, extract chr22 from an existing full-genome VCF:

```bash
docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools view -r chr22 \
    "/genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz" \
    -Oz -o "/genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz"

docker run --rm -v "${GENOME_DIR}:/genome" staphb/bcftools:1.21 \
  bcftools index -t "/genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz"
```

### Step 2: Pairwise Concordance

```bash
BCFTOOLS_IMAGE="staphb/bcftools:1.21"

# DeepVariant vs GATK on chr22
docker run --rm -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
  bcftools isec -p /genome/${SAMPLE}/benchmark/chr22_dv_vs_gatk \
    -f PASS \
    /genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz \
    /genome/${SAMPLE}/vcf_gatk/${SAMPLE}.vcf.gz

# DeepVariant vs FreeBayes on chr22
docker run --rm -v "${GENOME_DIR}:/genome" "$BCFTOOLS_IMAGE" \
  bcftools isec -p /genome/${SAMPLE}/benchmark/chr22_dv_vs_fb \
    -f PASS \
    /genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz \
    /genome/${SAMPLE}/vcf_freebayes/${SAMPLE}.vcf.gz
```

### Step 3: Truth Set Benchmark (if Using HG002)

```bash
# Benchmark each caller's chr22 output against GIAB truth set
for CALLER_DIR in vcf vcf_gatk vcf_freebayes; do
  LABEL=$(basename "$CALLER_DIR")
  docker run --rm \
    --cpus 4 --memory 16g \
    -v "${GENOME_DIR}:/genome" \
    jmcdani20/hap.py:v0.3.12 \
    /opt/hap.py/bin/hap.py \
      /genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
      "/genome/${SAMPLE}/${CALLER_DIR}/${SAMPLE}.vcf.gz" \
      -r /genome/reference/Homo_sapiens_assembly38.fasta \
      -f /genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
      -o "/genome/${SAMPLE}/benchmark/${LABEL}_vs_truth_chr22" \
      --threads 4 \
      -l chr22
done

# Strelka2 has a different output path structure
docker run --rm \
  --cpus 4 --memory 16g \
  -v "${GENOME_DIR}:/genome" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
    /genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    "/genome/${SAMPLE}/vcf_strelka2/results/variants/variants.vcf.gz" \
    -r /genome/reference/Homo_sapiens_assembly38.fasta \
    -f /genome/giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -o "/genome/${SAMPLE}/benchmark/vcf_strelka2_vs_truth_chr22" \
    --threads 4 \
    -l chr22
```

### Step 4: Compare Summary CSVs

```bash
echo "=== DeepVariant ==="
cat "${GENOME_DIR}/${SAMPLE}/benchmark/vcf_vs_truth_chr22.summary.csv"

echo "=== GATK ==="
cat "${GENOME_DIR}/${SAMPLE}/benchmark/vcf_gatk_vs_truth_chr22.summary.csv"

echo "=== FreeBayes ==="
cat "${GENOME_DIR}/${SAMPLE}/benchmark/vcf_freebayes_vs_truth_chr22.summary.csv"

echo "=== Strelka2 ==="
cat "${GENOME_DIR}/${SAMPLE}/benchmark/vcf_strelka2_vs_truth_chr22.summary.csv"
```

> **Note:** The file names above (`vcf_vs_truth_chr22.*`) are from the manual hap.py examples. The automated `benchmark-variants.sh --truth` script uses `happy_<Caller>.*` naming (e.g., `happy_DeepVariant.summary.csv`).

---

## Output Files Summary

After a full benchmarking run, expect these files in `${GENOME_DIR}/${SAMPLE}/benchmark/`:

```
benchmark/
  # Pairwise mode (benchmark-variants.sh without --truth):
  isec_DeepVariant_vs_GATK/        # Pairwise: DeepVariant vs GATK
    0000.vcf                       #   Unique to DeepVariant
    0001.vcf                       #   Unique to GATK
    0002.vcf                       #   Shared (DV perspective)
    0003.vcf                       #   Shared (GATK perspective)
    sites.txt                      #   Per-site presence table
  isec_DeepVariant_vs_FreeBayes/   # Pairwise: DeepVariant vs FreeBayes
    ...
  isec_GATK_vs_FreeBayes/          # Pairwise: GATK vs FreeBayes
    ...
  summary.txt                      # Human-readable results table
  comparison.tsv                   # Machine-readable TSV

  # Truth mode (benchmark-variants.sh --truth <vcf> --regions <bed>):
  happy_DeepVariant.summary.csv    # hap.py: DeepVariant precision/recall
  happy_DeepVariant.extended.csv
  happy_DeepVariant.vcf.gz         # TP/FP/FN annotated VCF
  happy_GATK.summary.csv           # hap.py: GATK precision/recall
  happy_GATK.extended.csv
  happy_GATK.vcf.gz
  happy_FreeBayes.summary.csv      # hap.py: FreeBayes precision/recall
  happy_FreeBayes.extended.csv
  happy_FreeBayes.vcf.gz
  summary.txt                      # Human-readable results table
  comparison.tsv                   # Machine-readable TSV
```

---

## Tips

- **PASS filter and normalization are applied automatically.** The `benchmark-variants.sh` script filters to PASS variants and normalizes VCFs with `bcftools norm -m-both` before running `bcftools isec`. This decomposes MNPs and left-aligns indels for fairer comparison. Without this, representation differences between callers (especially FreeBayes) depress Jaccard values.
- **hap.py handles normalization internally.** It decomposes complex variants and performs genotype-aware matching, so it is more accurate than raw `bcftools isec` for truth set evaluation.
- **chr22 is representative but not definitive.** It is a good proxy for genome-wide performance on autosomes, but does not test difficult regions (centromeres, segmental duplications, sex chromosomes).
- **Runtime scales with region size.** Full-genome hap.py takes 30-60 minutes per caller. chr22-only takes 1-3 minutes.

---

## Real-World Full-Genome Benchmark Results

These are actual results from a 30X WGS sample (Intel i5-14500, 64GB RAM).

### Variant Counts (Full Genome)

| Caller | Total | SNPs | Indels | MNPs/Other |
|---|---|---|---|---|
| DeepVariant 1.6.0 | 5,558,877 | — | — | — |
| GATK HC 4.6.1 | 4,662,941 | 3,829,614 | 836,747 | 0 |
| FreeBayes 1.3.6 | 20,059,293 | 18,825,881 | 867,375 | 625,201 |
| Strelka2 2.9.10 | 5,565,634 | — | — | — |
| TIDDIT 3.9.5 | 5,346 SVs | 2,299 DEL | 639 DUP | 204 INV |

FreeBayes calls **4x more variants** than DeepVariant. The ~16M FreeBayes-unique calls are overwhelmingly false positives. GATK is the most conservative. Strelka2 and DeepVariant are very close in total count.

### Pairwise Concordance (Full Genome, All Variants, Pre-Normalization)

**Note:** These numbers were computed on raw (unfiltered, unnormalized) VCFs. The `benchmark-variants.sh` script now applies PASS filtering and `bcftools norm -m-both` normalization before comparison, which will produce higher Jaccard values (especially for DV vs GATK, expected >0.95 for PASS SNPs).

| Caller A | Caller B | Shared | A Unique | B Unique | Jaccard |
|---|---|---|---|---|---|
| DeepVariant | GATK | 4,492,079 | 1,072,475 | 170,862 | **0.783** |
| DeepVariant | FreeBayes | 3,772,065 | 1,792,489 | 16,287,228 | 0.173 |
| GATK | FreeBayes | 3,495,334 | 1,167,607 | 16,563,959 | 0.165 |

Key observations:
- **DeepVariant vs GATK** share 4.5M variants (Jaccard 0.783). GATK has very few unique calls (171K) while DeepVariant has 1.07M unique — suggesting DeepVariant captures more real variants that GATK misses.
- **FreeBayes vs everything** has low concordance (~17%) because it calls 16.3M extra variants. These need aggressive quality filtering before use.
- **GATK is the most conservative** — nearly all GATK calls (96%) are also found by DeepVariant, making it a high-confidence subset.
- Filtering to **PASS SNPs only** would significantly increase concordance (expected Jaccard >0.95 for DV vs GATK).

### Runtime (Full Genome, 30X WGS)

| Caller | Runtime | Threads | Peak Memory |
|---|---|---|---|
| DeepVariant 1.6.0 | ~3-5 hours | 8 | 32 GB |
| GATK HC 4.6.1 | 8.6 hours (518 min) | 8 | 4.3 GB |
| FreeBayes 1.3.6 | 9.3 hours | 1 (single-threaded) | 12.8 GB |
| Strelka2 2.9.10 | 72 min | 8 | 1 GB |
| TIDDIT 3.9.5 | 7 min | 4 | <1 GB |

FreeBayes is the bottleneck — single-threaded with no parallelism flag. It also requires **32 GB memory** for full-genome runs (peaked at 12.8 GB but grows unpredictably through complex regions). Plan accordingly when running all callers.
