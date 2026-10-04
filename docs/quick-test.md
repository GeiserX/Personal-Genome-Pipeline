# Quick Test: Verify Your Setup Before Running on Real Data

Don't commit 12+ hours and 500 GB to a full pipeline run before verifying everything works. This guide shows how to test with a small public dataset in under 30 minutes.

The project's own scripted check is wider: CI runs most of the steps on a small slice of a public genome and checks what they write. [Testing](testing.md) explains that end-to-end run, and you can run the same scripts yourself.

---

## Option A: Chromosome 22 Only (Recommended)

Chromosome 22 is the smallest autosome (~51 MB), so a chr22 extract runs in minutes instead of hours.

Option A runs three VCF-only steps: the ClinVar screen (6), PharmCAT (7) and ROH analysis (11). It shows that Docker, the reference data and the scripts work together. It does not exercise alignment, variant calling or any BAM-dependent step; Option B adds two of those.

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
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

wget -O ${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz.tbi \
  "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
```

**Note:** This is the full-genome GIAB VCF (~250 MB). For a chr22-only test, extract just chr22:

```bash
source versions.env   # from the repository root
docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -r chr22 \
    /genome/${SAMPLE}/vcf/${SAMPLE}.vcf.gz \
    -Oz -o /genome/${SAMPLE}/vcf/${SAMPLE}_chr22.vcf.gz

docker run --rm --user root \
  -v "${GENOME_DIR}:/genome" \
  "${BCFTOOLS_IMAGE}" \
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

If you want to test BAM-dependent steps, you need an indexed BAM at `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` (plus `.bai`), which is where the scripts read it. A chr22-only BAM of a 30x genome is about 560 MB.

The command below reads only the chr22 reads of the 1000 Genomes 30x NA12878 alignment, the same person as the Option A VCF. The alignment is a CRAM file with an index, so samtools fetches just the chr22 part (a few hundred MB) instead of the 16 GB file. The first container uses the biocontainers samtools image because it ships CA certificates and can fetch over `https://`; the `staphb/samtools` image used elsewhere has none and fails with "Libcurl reported error 60". chr22 is the same sequence in every GRCh38 file, so the pipeline's reference decodes it.

The 1000 Genomes CRAM was aligned to a GRCh38 file with ALT, HLA and decoy contigs (3,366 sequences), and its header lists all of them, which Manta and `validate-setup.sh` refuse against the pipeline's 195-sequence reference. So the `awk` between the two containers replaces the header's sequence lines with the reference's own (from its `.dict`) and drops the few reads whose mate sits on a contig the reference does not have. That is enough for this mechanics test; a real sample from another reference is realigned instead ([realignment](realignment.md)).

```bash
source versions.env   # from the repository root
REF_FASTA=reference/GRCh38_no_alt_analysis_set.fasta   # see 00-reference-setup.md#the-reference-path-on-every-page
# Uses GENOME_DIR and SAMPLE from Option A; needs the reference FASTA, .fai and .dict
mkdir -p ${GENOME_DIR}/${SAMPLE}/aligned

docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  -w /tmp \
  quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
  samtools view -h -@ 4 -T "/genome/${REF_FASTA}" \
    https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram chr22 \
| awk -v dict="${GENOME_DIR}/${REF_FASTA%.fasta}.dict" '
    BEGIN { print "@HD\tVN:1.6\tSO:coordinate"
            while ((getline l < dict) > 0) if (l ~ /^@SQ/) { print l; split(l, f, "\t"); sub(/^SN:/, "", f[2]); keep[f[2]] = 1 } }
    /^@HD/ || /^@SQ/ { next }
    /^@/ { print; next }
    $7 == "=" || $7 == "*" || ($7 in keep) { print }' \
| docker run --rm -i --user "$(id -u):$(id -g)" \
  -v "${GENOME_DIR}:/genome" \
  "${SAMTOOLS_IMAGE}" \
  bash -c "samtools view -b -@ 4 -o /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam - &&
    samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam"
```

**Alternative: extract chr22 from a full BAM you already have** and that was aligned to the pipeline's reference (a BAM from another reference needs a [realignment](realignment.md) first). Put it at `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` first; the commands below move it aside before the link step replaces that name:
```bash
source versions.env   # from the repository root
cd ${GENOME_DIR}/${SAMPLE}/aligned
# Skipped on a rerun, when _sorted.bam already is the chr22 BAM (-ef follows links)
if [ ! -e ${SAMPLE}_full.bam ] && [ ! ${SAMPLE}_sorted.bam -ef ${SAMPLE}_chr22.bam ]; then
  mv ${SAMPLE}_sorted.bam     ${SAMPLE}_full.bam
  mv ${SAMPLE}_sorted.bam.bai ${SAMPLE}_full.bam.bai
fi

docker run --rm --user root \
  --cpus 4 --memory 4g \
  -v "${GENOME_DIR}:/genome" \
  "${SAMTOOLS_IMAGE}" \
  bash -c "set -euo pipefail
    samtools view -b -o /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam \
      /genome/${SAMPLE}/aligned/${SAMPLE}_full.bam chr22
    samtools index /genome/${SAMPLE}/aligned/${SAMPLE}_chr22.bam"
```

Either way, point the pipeline at the chr22 BAM, as Option A does for the VCF:
```bash
cd ${GENOME_DIR}/${SAMPLE}/aligned
ln -sfn ${SAMPLE}_chr22.bam     ${SAMPLE}_sorted.bam
ln -sfn ${SAMPLE}_chr22.bam.bai ${SAMPLE}_sorted.bam.bai
ls -lL ${SAMPLE}_sorted.bam ${SAMPLE}_sorted.bam.bai
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
3. Place your files according to the [directory structure](getting-started.md#directory-structure)
4. Run `./scripts/validate-setup.sh $SAMPLE` for a comprehensive pre-flight check
5. Start with the [Quick Start](getting-started.md#quick-start) path that matches your input data
