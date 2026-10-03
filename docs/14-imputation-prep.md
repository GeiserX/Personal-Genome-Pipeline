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
export GENOME_DIR=/path/to/your/data
./scripts/14-imputation-prep.sh your_sample
```

The script starts one container and makes one `bcftools view` pass per chromosome, chr1-22 and chrX. Each file keeps the PASS records (and records with no filter), is written with its index under a temporary name, and is renamed when both are complete: the old index is removed first and the new one is moved in last, so an index never sits beside a VCF it was not built from. A chromosome the VCF has no record on gets no file (a file an earlier run left for it is removed), and the log names it. The commands it runs:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
mkdir -p ${GENOME_DIR}/${SAMPLE}/imputation/mis_ready

docker run --rm \
  -v ${GENOME_DIR}/${SAMPLE}:/data \
  "${BCFTOOLS_IMAGE}" \
  bash -c 'for chr in $(seq -f "chr%g" 1 22) chrX; do
    bcftools view -f PASS,. -r "$chr" -Oz --write-index=tbi \
      -o /data/imputation/mis_ready/'"${SAMPLE}"'_${chr}.vcf.gz /data/vcf/'"${SAMPLE}"'.vcf.gz
  done'

# Output: 23 per-chromosome VCFs, each with its .tbi, in ${SAMPLE}/imputation/mis_ready/
```

**By default the input is variant-only.** The step reads the VCF from step 3, which lists only the sites where the sample differs from the reference. A site that is missing from it is not a confirmed homozygous-reference genotype, and an imputation server treats it as missing.

**With panel sites, hom-ref genotypes come from the gVCF.** Step 3 also writes a gVCF (`vcf/${SAMPLE}.g.vcf.gz`), which records where the sample matches the reference. Give the step the reference panel's sites and it genotypes each of them from the gVCF: a variant call, a 0/0 call where the sample matches the reference, and nothing where the site was not covered. The per-chromosome files then hold those sites only, which are the sites the server imputes from.

```bash
# CHROM<TAB>POS per line, or a VCF of the panel's sites, inside GENOME_DIR
IMPUTATION_SITES=${GENOME_DIR}/reference/panel_sites.tsv ./scripts/14-imputation-prep.sh your_sample
```

A panel of tens of millions of sites needs a few GB of memory for this. Without `IMPUTATION_SITES` and with a gVCF present, the step prints how to use it.

## Server Options
| Server | Panel | Samples | Build | URL |
|---|---|---|---|---|
| Michigan (MIS) | HRC r1.1 | 32,470 | GRCh37/38 | imputationserver.sph.umich.edu |
| TOPMed | TOPMed r2 | 132,070 | GRCh38 native | imputation.biodatacatalyst.nhlbi.nih.gov |

## Your data leaves the machine here

This is the only step whose output is meant to be sent somewhere else. The other steps keep your genome on your disk; an imputation server receives your genotypes for every chromosome you upload. Uploading is a separate decision you make, not something the step does.

Before you upload, read the server's data policy. The Michigan server's [security page](https://genepi.github.io/michigan-imputationserver/data-sensitivity/) says it deletes the input once it is no longer needed, keeps only the number of samples and markers, and encrypts the results with a one-time password. Its [getting started guide](https://genepi.github.io/michigan-imputationserver/getting-started/) says the results are deleted 7 days after the job ends. Neither page gives an exact time for deleting inputs or backups, or says whether data is encrypted at rest. Check the TOPMed server's own terms; they are not the same service. Both need an account, so the upload is tied to your email address.

## Important Notes
- MIS requires a **minimum of 20 samples per job** — a single WGS sample is useful mainly for phasing, not imputation
- **TOPMed r2 panel is recommended for European ancestry** (132K samples, GRCh38 native — no liftover needed)
- Registration is required at the imputation server before submitting jobs
- Upload per-chromosome VCF files (not the whole-genome file)
- Servers accept `.vcf.gz` format — ensure files are bgzipped (bcftools output is bgzipped by default)
- chrX is prepared too; servers usually take it as a separate job with ploidy-aware settings
- Results include phased haplotypes and imputation quality scores (R-squared) — filter imputed variants with R2 < 0.3
