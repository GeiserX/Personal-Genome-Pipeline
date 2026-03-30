# Step 14: Imputation Preparation

## What This Does
Prepares a WGS VCF for submission to the Michigan Imputation Server (MIS) or TOPMed Imputation Server. Splits the VCF by chromosome, filters to PASS variants, and converts to the required format.

## Why
Imputation servers statistically infer missing genotypes using large reference panels. For WGS data, imputation is primarily useful for **phasing** (determining which alleles are on the same chromosome) rather than filling in missing variants. Phased data is required for haplotype-level analyses and accurate PRS calculation.

## Tool
- **bcftools** (samtools/bcftools) — for VCF filtering, splitting, and indexing

## Docker Image
```
staphb/bcftools:1.21
```

## Command
```bash
SAMPLE=your_sample
GENOMA_DIR=/path/to/genome/data

# Step 1: Filter to PASS variants only
docker run --rm \
  -v ${GENOMA_DIR}/${SAMPLE}/vcf:/data \
  staphb/bcftools:1.21 \
  bcftools view -f PASS \
    /data/${SAMPLE}.vcf.gz \
    -Oz -o /data/${SAMPLE}_pass.vcf.gz

# Step 2: Index the filtered VCF
docker run --rm \
  -v ${GENOMA_DIR}/${SAMPLE}/vcf:/data \
  staphb/bcftools:1.21 \
  bcftools index -t /data/${SAMPLE}_pass.vcf.gz

# Step 3: Split by chromosome (chr1-22, autosomes only)
for CHR in $(seq 1 22); do
  docker run --rm \
    -v ${GENOMA_DIR}/${SAMPLE}/vcf:/data \
    staphb/bcftools:1.21 \
    bcftools view -r chr${CHR} \
      /data/${SAMPLE}_pass.vcf.gz \
      -Oz -o /data/imputation/chr${CHR}.vcf.gz

  docker run --rm \
    -v ${GENOMA_DIR}/${SAMPLE}/vcf:/data \
    staphb/bcftools:1.21 \
    bcftools index -t /data/imputation/chr${CHR}.vcf.gz
done

# Output: 22 per-chromosome VCF files in /data/imputation/
```

## Server Options
| Server | Panel | Samples | Build | URL |
|---|---|---|---|---|
| Michigan (MIS) | HRC r1.1 | 32,470 | GRCh37/38 | imputationserver.sph.umich.edu |
| TOPMed | TOPMed r2 | 132,070 | GRCh38 native | imputation.biodatacatalyst.nhlbi.nih.gov |

## Important Notes
- MIS requires a **minimum of 20 samples per job** — a single WGS sample is useful mainly for phasing, not imputation
- **TOPMed r2 panel is recommended for European ancestry** (132K samples, GRCh38 native — no liftover needed)
- Registration is required at the imputation server before submitting jobs
- Upload per-chromosome VCF files (not the whole-genome file)
- Servers accept `.vcf.gz` format — ensure files are bgzipped (bcftools output is bgzipped by default)
- Sex chromosomes (chrX) can be submitted separately with ploidy-aware settings
- Results include phased haplotypes and imputation quality scores (R-squared) — filter imputed variants with R2 < 0.3
- Create the output directory before running: `mkdir -p ${GENOMA_DIR}/${SAMPLE}/vcf/imputation`
