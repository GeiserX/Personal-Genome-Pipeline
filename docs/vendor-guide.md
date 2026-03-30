# Vendor Compatibility Guide

Your genomics data can come from many different providers. This guide explains what format each vendor delivers, how to get it into the pipeline, and what to watch out for.

## Quick Reference

| Vendor | Data Type | Typical Format | Genome Build | Pipeline Entry | Price (2025-2026) |
|---|---|---|---|---|---|
| Nebula / DNA Complete | 30X WGS | FASTQ + VCF | GRCh38 | Path A or C | $495 (30X) |
| Dante Labs | 30X WGS | FASTQ + BAM + VCF | GRCh38 | Any | $300-600 |
| Sequencing.com | 30X WGS | FASTQ + BAM + VCF | GRCh38 | Any | $399-799 |
| Novogene / BGI | 30X WGS | FASTQ | GRCh38 | Path A | $200-400 |
| Illumina DRAGEN (clinical) | 30X WGS | ORA / BAM + VCF | GRCh38 | Path D / B / C | $300-1000 |
| Full Genomes Corporation | 30X WGS | BAM + VCF | GRCh38 | Path B or C | ~$1000 |
| Oxford Nanopore | Long-read WGS | POD5 + BAM | GRCh38 | Not supported | $1000-3000 |
| PacBio HiFi | Long-read WGS | HiFi BAM | GRCh38 | Not supported | $1000-2000 |
| 23andMe | Genotyping array | TSV (~640K SNPs) | GRCh37 | Partial | $79-229 |
| AncestryDNA | Genotyping array | TSV (~700K SNPs) | GRCh37 | Partial | $99-199 |
| MyHeritage | Genotyping array | CSV (~643K SNPs) | GRCh37 | Partial | $79-199 |

---

## Illumina Short-Read WGS (Most Vendors)

Most consumer WGS vendors use Illumina sequencing platforms (NovaSeq 6000, NovaSeq X Plus). The data comes in standard formats that this pipeline handles directly.

### What You'll Receive

- **FASTQ files** (`.fastq.gz`): Raw sequencing reads, paired-end (R1 + R2). Typically 60-90 GB compressed per sample.
- **BAM file** (`.bam`): Aligned reads. 80-120 GB per sample. Already aligned to GRCh38.
- **VCF file** (`.vcf.gz`): Variant calls. 80-200 MB per sample. ~4.5-5.5 million variants.

### Getting Started

1. **If you have FASTQ:** Copy R1 and R2 files to `${GENOME_DIR}/${SAMPLE}/fastq/`. Start with step 2 (alignment).
2. **If you have BAM:** Copy to `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam`. Make sure the BAM index (`.bai`) is present. Start with step 3 (variant calling).
3. **If you have VCF:** Copy to `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz`. Make sure the index (`.tbi`) is present. Start with step 6 (ClinVar screen).

---

## Nebula Genomics / DNA Complete

Nebula was acquired by ProPhase Labs and rebranded as DNA Complete. They use **MGI/DNBSEQ** sequencing (BGI technology), not Illumina.

### What's Different

- **Read names** follow BGI format instead of Illumina format. This is purely cosmetic -- all alignment tools handle it correctly.
- **Quality scores** are the same encoding (Phred+33). No conversion needed.
- **Adapter sequences** differ from Illumina. If you're trimming adapters (not required for this pipeline), use the MGI adapter sequences.

### Data Access

Download your data from the DNA Complete portal. Both FASTQ and VCF are available. **Important:** Data access may require an active subscription. Download everything immediately after receiving results.

### Entry Point

Use Path A (FASTQ) for the most complete analysis, or Path C (VCF) if you only want annotation.

---

## Dante Labs

Italian company using Illumina NovaSeq. Standard Illumina output.

### Known Issues

- Delivery times have historically been unpredictable (weeks to months).
- Data download portal can be unreliable. Try different browsers or times of day.
- Some older samples were sequenced on BGI platforms. Check your delivery email for platform details.

### Entry Point

All three formats (FASTQ, BAM, VCF) are typically provided. Use whichever entry point matches your goals.

---

## Novogene / BGI Direct

Research-focused sequencing service. Cheapest option for 30X WGS (~$200-400).

### What You'll Receive

Typically FASTQ only (paired-end, gzipped). BAM and VCF may be available at extra cost or on request.

### BGI/MGI FASTQ Quirks

BGI read names look different from Illumina:
```
# Illumina:
@A00123:456:HXXXXXXX:1:1101:12345:67890 1:N:0:ATCGATCG

# BGI/DNBSEQ:
@V350012345L1C001R00100000001/1
```

This does **not** affect any pipeline step. BWA, minimap2, DeepVariant, and all other tools only use the sequence and quality lines.

### Entry Point

Path A (FASTQ). You'll need to run the full pipeline from alignment.

---

## Illumina DRAGEN (Clinical/Hospital)

If your WGS was done through a clinical lab or hospital, they likely used Illumina's DRAGEN pipeline for processing.

### ORA Format

Some labs deliver FASTQ files compressed in Illumina's proprietary **ORA format** (~5x smaller than gzipped FASTQ). You need the `orad` decompressor:

```bash
# Step 1 in this pipeline handles ORA decompression
./scripts/01-ora-to-fastq.sh $SAMPLE
```

See [docs/01-ora-to-fastq.md](01-ora-to-fastq.md) for details on obtaining the `orad` binary.

### DRAGEN VCF Notes

DRAGEN VCFs include non-standard annotations (e.g., `DRAGEN:` prefixed INFO fields). These are ignored by standard tools but may cause warnings. This is harmless.

If your lab provided a DRAGEN-called VCF, you can skip steps 2-3 and go directly to analysis (Path C). However, re-calling variants with DeepVariant (step 3) from the BAM may find additional variants that DRAGEN missed, especially in difficult regions.

### Entry Point

- **ORA files:** Path D (decompress first)
- **BAM:** Path B (variant calling + analysis)
- **VCF:** Path C (analysis only)

---

## CRAM Files

Some providers deliver CRAM instead of BAM (40-60% smaller). Convert to BAM first:

```bash
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.21 \
  samtools view -b \
    -T /genome/reference/Homo_sapiens_assembly38.fasta \
    -o /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
    /genome/${SAMPLE}/aligned/${SAMPLE}.cram

# Index the BAM
docker run --rm \
  -v ${GENOME_DIR}:/genome \
  staphb/samtools:1.21 \
  samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam
```

**Important:** CRAM decoding requires the same reference genome used for encoding. This pipeline uses `Homo_sapiens_assembly38.fasta` (GRCh38). If your CRAM was encoded against a different reference, you'll get errors.

---

## Long-Read Sequencing (Not Supported)

### Oxford Nanopore (MinION / PromethION)

Nanopore produces long reads (10-50 kb average) with different error profiles than Illumina. This pipeline's tools are optimized for short reads and will produce incorrect results with nanopore data.

**What you'd need instead:**
- Alignment: `minimap2 -ax map-ont` (not the default short-read preset)
- Variant calling: **Clair3** (not DeepVariant, though DeepVariant has an ONT model)
- SV calling: **Sniffles2** or **cuteSV** (not Manta)
- Basecalling: **Dorado** from raw POD5/FAST5 signal

### PacBio HiFi

PacBio HiFi reads are highly accurate (>Q20) and 10-20 kb long. Different tools required:
- Alignment: `pbmm2` or `minimap2 -ax map-hifi`
- Variant calling: **DeepVariant** (PacBio model) or **PEPPER-Margin-DeepVariant**
- SV calling: `pbsv` or Sniffles2

> A long-read pipeline branch may be added in the future. For now, these are the recommended tools.

---

## Genotyping Arrays (Partial Support)

### 23andMe, AncestryDNA, MyHeritage

These services use genotyping chips that test ~600,000-700,000 specific positions. This is **not whole genome sequencing** -- it covers ~0.02% of your genome.

### What You Can Do

1. **Convert to VCF** using tools like `plink` or custom scripts
2. **Run ClinVar screening** (step 6) on the converted VCF
3. **Run PharmCAT** (step 7) for pharmacogenomics -- though coverage of PGx positions varies by chip version

### What You Cannot Do

- Alignment (step 2) -- no raw reads exist
- Variant calling (step 3) -- genotypes are already called by the array
- Structural variants (steps 4, 18, 19) -- arrays can't detect SVs
- Telomere analysis (step 10) -- requires actual sequencing reads
- Functional annotation (step 13) -- limited value with only ~600K variants

### Imputation

You can dramatically increase coverage by imputing (predicting) missing genotypes using the [Michigan Imputation Server](https://imputationserver.sph.umich.edu/) or [TOPMed Imputation Server](https://imputation.biodatacatalyst.nhlbi.nih.gov/). This can fill in ~40 million positions from your 600K, but the accuracy depends on your ancestry and the reference panel used.

**Note:** Imputation servers typically require 20+ samples per job. Single-sample submissions may be rejected.

---

## Genome Build: GRCh37 (hg19) vs GRCh38 (hg38)

This pipeline uses **GRCh38 (hg38)** exclusively. If your data is on an older build:

### How to Check Your Build

```bash
# For BAM files -- look at the reference in the header:
samtools view -H your_file.bam | grep "^@SQ" | head -3

# GRCh38 uses "chr" prefix: SN:chr1, SN:chr2, etc.
# GRCh37 may lack "chr" prefix: SN:1, SN:2, etc.
# (Some GRCh37 builds do use chr prefix -- check chromosome lengths to be sure)
# chr1 length: GRCh38 = 248,956,422; GRCh37 = 249,250,621
```

### Converting from GRCh37 to GRCh38

**Best approach (recommended):** Extract FASTQ from BAM and re-align to GRCh38:
```bash
# Extract paired-end FASTQ from BAM
docker run --rm -v ${GENOME_DIR}:/genome staphb/samtools:1.21 \
  bash -c "samtools sort -n /genome/${SAMPLE}/old_hg19.bam | \
           samtools fastq -1 /genome/${SAMPLE}/fastq/${SAMPLE}_R1.fastq.gz \
                          -2 /genome/${SAMPLE}/fastq/${SAMPLE}_R2.fastq.gz -"

# Then run the pipeline from step 2 (alignment to GRCh38)
./scripts/02-alignment.sh $SAMPLE
```

**Alternative (quicker but less accurate):** Use Picard LiftoverVcf to convert VCF coordinates. This can introduce artifacts at complex regions and is not recommended for clinical use.

---

## File Size Reference

Know what to expect before downloading:

| File Type | Typical Size (30X WGS) | Notes |
|---|---|---|
| FASTQ (gzipped, paired) | 60-90 GB | Two files: R1 + R2 |
| ORA (Illumina compressed) | 15-20 GB | Same data as FASTQ, ~5x smaller |
| BAM (aligned) | 80-120 GB | Largest single file |
| CRAM (compressed aligned) | 40-60 GB | 40-60% smaller than BAM |
| VCF (variants) | 80-200 MB | Relatively small |
| gVCF (with reference blocks) | 3-10 GB | Much larger than VCF |
