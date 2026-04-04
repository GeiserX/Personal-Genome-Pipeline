# Quick Test: Verify Your Setup Before Running on Real Data

Don't commit 12+ hours and 500 GB to a full pipeline run before verifying everything works. This guide shows how to test with a small public dataset in under 30 minutes.

---

## Option A: Chromosome 22 Only (Recommended)

Chromosome 22 is the smallest autosome (~51 MB). Running the pipeline on chr22 data tests all tools with ~2% of the data, completing in minutes instead of hours.

### Step 1: Download Test Data

We'll use the Genome in a Bottle NA12878 sample (NIST reference standard):

```bash
export GENOME_DIR=/path/to/test/data
export SAMPLE=test_na12878
mkdir -p ${GENOME_DIR}/${SAMPLE}/vcf ${GENOME_DIR}/reference

# Download the GRCh38 reference (required, ~3.1 GB — skip if you already have it)
# See docs/00-reference-setup.md for full instructions

# Download a pre-called chr22 VCF from Genome in a Bottle
wget -O ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

wget -O ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
```

**Note:** This is the full-genome GIAB VCF (~250 MB). For a chr22-only test, extract just chr22:

```bash
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools view -r chr22 \
    /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    -Oz -o /genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  staphb/bcftools:1.21 \
  bcftools index -t /genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz
```

### Step 2: Run VCF-Only Steps

If you extracted chr22 above, back up the original and symlink the chr22 extract:
```bash
# Preserve the original full VCF (idempotent — skips if already backed up)
cd ${GENOME_DIR}/${SAMPLE}/vcf
if [ ! -f ${SAMPLE}_full.vcf.gz ]; then
  mv ${SAMPLE}.vcf.gz     ${SAMPLE}_full.vcf.gz
  mv ${SAMPLE}.vcf.gz.tbi ${SAMPLE}_full.vcf.gz.tbi
fi

# Point the pipeline at the chr22 extract
ln -sfn ${SAMPLE}_chr22.vcf.gz     ${SAMPLE}.vcf.gz
ln -sfn ${SAMPLE}_chr22.vcf.gz.tbi ${SAMPLE}.vcf.gz.tbi

# To restore the full VCF later:
#   rm ${SAMPLE}.vcf.gz ${SAMPLE}.vcf.gz.tbi
#   mv ${SAMPLE}_full.vcf.gz ${SAMPLE}.vcf.gz
#   mv ${SAMPLE}_full.vcf.gz.tbi ${SAMPLE}.vcf.gz.tbi
```

Then run the VCF-only steps (they expect `${SAMPLE}.vcf.gz`):
```bash
# ClinVar screen (~1 min)
./scripts/06-clinvar-screen.sh ${SAMPLE}

# PharmCAT (~2 min)
./scripts/07-pharmacogenomics.sh ${SAMPLE}

# ROH analysis (~1 min)
./scripts/11-roh-analysis.sh ${SAMPLE}
```

### Step 3: Verify Output

```bash
# Should see ClinVar hits
ls -la ${GENOME_DIR}/${SAMPLE}/clinvar/

# Should see PharmCAT HTML report (written alongside VCF)
ls -la ${GENOME_DIR}/${SAMPLE}/vcf/*.report.html

# Should see ROH output
ls -la ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt
```

If all three steps produce output, your Docker setup, reference data, and pipeline scripts are working correctly.

---

## Option B: Full Pipeline Test with Minimal BAM

If you want to test the BAM-dependent steps (Manta, CNVnator, Delly, TelomereHunter, indexcov), you need a BAM file. A chr22-only BAM is small enough (~3-4 GB) for a quick test:

```bash
# Download chr22 reads for NA12878 from 1000 Genomes
# (This is a ~3 GB download)
mkdir -p ${GENOME_DIR}/${SAMPLE}/aligned

wget -O ${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam \
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1000G_2504_highcov_chr22.bam"

# This URL may change. If it doesn't work, see:
# https://www.internationalgenome.org/data-portal/sample/NA12878
# Download any chr22 BAM for GRCh38
```

**Alternative: Extract chr22 from a full BAM** (if you already have one):
```bash
docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  staphb/samtools:1.20 \
  bash -c "samtools view -b /genome/${SAMPLE}/aligned/${SAMPLE}_sorted.bam chr22 \
    > /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam && \
    samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam"
```

Then test BAM-dependent steps:
```bash
./scripts/04-manta.sh ${SAMPLE}       # Manta SVs (~2 min on chr22)
./scripts/16-indexcov.sh ${SAMPLE}     # Coverage QC (~1 sec)
```

---

## What to Expect

### VCF-only steps (Option A)

| Step | Expected Runtime | Expected Output |
|---|---|---|
| ClinVar screen | < 1 min | 0-5 pathogenic hits for NA12878 |
| PharmCAT | 1-3 min | HTML report with gene calls |
| ROH analysis | < 1 min | ROH segments file |

### BAM-dependent steps (Option B)

| Step | Expected Runtime | Expected Output |
|---|---|---|
| Manta | 1-3 min | Small VCF with chr22 SVs |
| indexcov | < 5 sec | Coverage plots for chr22 |

---

## Common Test Failures

**"No such file" errors:** Check that `GENOME_DIR` is set and the downloaded files are in the expected paths.

**Docker image pull fails:** Some images are on quay.io, which occasionally has downtime. Wait and retry, or check the exact image tag in the step's documentation.

**0 ClinVar hits on test data:** The GIAB benchmark VCF may not overlap with the ClinVar pathogenic subset if you are using a chr22-only extract. This is expected — the test validates that the pipeline mechanics work, not that there are clinical findings.

**PharmCAT produces empty report:** The GIAB VCF uses a different sample name internally. PharmCAT should still run but may show "No data" for some genes. This is expected for test data.

---

## Once the Test Passes

You're ready to run on your real data:

1. Set `GENOME_DIR` to your actual data directory
2. Set `SAMPLE` to your actual sample name
3. Place your files according to the [directory structure](../README.md#directory-structure)
4. Run `./scripts/validate-setup.sh $SAMPLE` for a comprehensive pre-flight check
5. Start with the [Quick Start](../README.md#quick-start) path that matches your input data
