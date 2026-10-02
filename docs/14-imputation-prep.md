# Step 14: Imputation Preparation

## What This Does
Prepares a WGS VCF for submission to the Michigan Imputation Server (MIS) or TOPMed Imputation Server. Splits the VCF by chromosome, filters to PASS variants, and converts to the required format.

## Why
Imputation servers statistically infer missing genotypes using large reference panels. For WGS data, imputation is primarily useful for **phasing** (determining which alleles are on the same chromosome) rather than filling in missing variants. Phased data is required for haplotype-level analyses and accurate PRS calculation.

## Tool
- **bcftools** (samtools/bcftools) — for VCF filtering, splitting, and indexing

## Docker Image
- `BCFTOOLS_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Command
```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# Step 1: Filter to PASS variants only
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  "${BCFTOOLS_IMAGE}" \
  bcftools view -f PASS \
    /data/${SAMPLE}.vcf.gz \
    -Oz -o /data/${SAMPLE}_pass.vcf.gz

# Step 2: Index the filtered VCF
docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  "${BCFTOOLS_IMAGE}" \
  bcftools index -t /data/${SAMPLE}_pass.vcf.gz

# Step 3: Split by chromosome (chr1-22, autosomes only)
for CHR in $(seq 1 22); do
  docker run --rm \
    -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
    "${BCFTOOLS_IMAGE}" \
    bcftools view -r chr${CHR} \
      /data/${SAMPLE}_pass.vcf.gz \
      -Oz -o /data/imputation/chr${CHR}.vcf.gz

  docker run --rm \
    -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
    "${BCFTOOLS_IMAGE}" \
    bcftools index -t /data/imputation/chr${CHR}.vcf.gz
done

# Output: 22 per-chromosome VCF files in /data/imputation/
```

## Server Options
| Server | Panel | Samples | Build | URL |
|---|---|---|---|---|
| Michigan (MIS) | HRC r1.1 | 32,470 | GRCh37/38 | imputationserver.sph.umich.edu |
| TOPMed | TOPMed r2 | 132,070 | GRCh38 native | imputation.biodatacatalyst.nhlbi.nih.gov |

## Your data leaves the machine here

This is the only step whose output is meant to be sent somewhere else. The other steps keep your genome on your disk; an imputation server receives your genotypes for every chromosome you upload. Uploading is a separate decision you make, not something the step does.

Before you upload, read the server's data policy. The Michigan server's [security page](https://imputationserver.readthedocs.io/en/latest/data-sensitivity/) says it deletes the input once the job no longer needs it, keeps only the number of samples and markers, encrypts the results with a one-time password, and keeps them for 7 days for download. Check the TOPMed server's own terms; they are not the same service. Both need an account, so the upload is tied to your email address.

## Important Notes
- MIS requires a **minimum of 20 samples per job** — a single WGS sample is useful mainly for phasing, not imputation
- **TOPMed r2 panel is recommended for European ancestry** (132K samples, GRCh38 native — no liftover needed)
- Registration is required at the imputation server before submitting jobs
- Upload per-chromosome VCF files (not the whole-genome file)
- Servers accept `.vcf.gz` format — ensure files are bgzipped (bcftools output is bgzipped by default)
- Sex chromosomes (chrX) can be submitted separately with ploidy-aware settings
- Results include phased haplotypes and imputation quality scores (R-squared) — filter imputed variants with R2 < 0.3
- Create the output directory before running: `mkdir -p ${GENOME_DIR}/${SAMPLE}/vcf/imputation`
