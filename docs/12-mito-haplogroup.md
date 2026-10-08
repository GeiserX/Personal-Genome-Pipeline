# Step 12: Mitochondrial Haplogroup

## What This Does
Determines the maternal lineage from mitochondrial DNA (mtDNA) variants, and checks the mtDNA for a second person's reads (contamination).

## Why
The mitochondrial haplogroup reveals deep maternal ancestry. Some haplogroups have known population-level health associations (e.g., longevity, metabolic traits).

A sample with DNA of two people shows two mtDNA haplogroups at once: the variants of the second one appear as heteroplasmies at a similar level. haplocheck looks for that pattern. It is a check of the sample, not of the person: `YES` means the calls of this sample, above all its heterozygous ones, may mix two people.

## Tool
- **haplogrep3** (Medical University of Innsbruck): the haplogroup
- **haplocheck** 1.3.3 (Weissensteiner et al., Genome Res 2021): the contamination check

## Docker Image
- `HAPLOGREP3_IMAGE`, `HAPLOCHECK_IMAGE`, and `BCFTOOLS_IMAGE` to prepare the calls

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

> The image is the Bioconda build of haplogrep3 3.2.2. It classifies without a network: the image test runs it with `--network none`. It replaced a digest-pinned build of 3.2.1 from a personal Docker Hub account; on the image test's chrM calls of the HG002 fixture both give H5a7 with quality 1.0000 and the same found and remaining polymorphisms.

## Input

The step reads, first match:

1. `${SAMPLE}/mito/${SAMPLE}_chrM_filtered.vcf.gz`, step 20's Mutect2 calls in mitochondrial mode. Mutect2 is the caller made for chrM: it reports each allele's fraction, so haplogrep3 sees heteroplasmies and haplocheck can run. The step keeps the PASS records (possible NuMTs are filtered by step 20) and splits multi-allelic records.
2. Otherwise the chrM records of `${SAMPLE}/vcf/${SAMPLE}.vcf.gz` (step 03). DeepVariant calls chrM as a diploid nuclear contig, so this is the fallback, and haplocheck does not run on it.

Run step 20 before step 12 for the Mutect2 input. The Nextflow pipeline does the same: with `mito_variants` in `--tools` (it is in the default set), `mito_haplogroup` reads `MITO_VARIANTS`' calls and runs `HAPLOCHECK`.

## Command
```bash
./scripts/12-mito-haplogroup.sh your_sample
```

By hand, on step 20's calls:

```bash
source versions.env   # from the repository root
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data
M=/genome/${SAMPLE}/mito

# PASS records, one allele per record
docker run --rm --network none -v "${GENOME_DIR}:/genome" "${BCFTOOLS_IMAGE}" bash -c \
  "bcftools view -f PASS ${M}/${SAMPLE}_chrM_filtered.vcf.gz | bcftools norm -m-any -Oz -o ${M}/${SAMPLE}_chrM_for_haplogroup.vcf.gz"

docker run --rm --network none -v "${GENOME_DIR}:/genome" "${HAPLOGREP3_IMAGE}" \
  haplogrep3 classify --tree phylotree-fu-rcrs@1.2 \
    --input ${M}/${SAMPLE}_chrM_for_haplogroup.vcf.gz \
    --output ${M}/${SAMPLE}_haplogroup.txt --extend-report

docker run --rm --network none -v "${GENOME_DIR}:/genome" "${HAPLOCHECK_IMAGE}" \
  haplocheck --out ${M}/${SAMPLE}_haplocheck.txt ${M}/${SAMPLE}_chrM_for_haplogroup.vcf.gz
```

## Output

| File | Contents |
|---|---|
| `mito/${SAMPLE}_haplogroup.txt` | haplogrep3's haplogroup and quality |
| `mito/${SAMPLE}_haplocheck.txt` | haplocheck's report (Mutect2 input only): `Contamination Status` (`YES`, `NO` or `ND`, not determined), `Contamination Level`, the major and minor haplogroups it found |
| `mito/${SAMPLE}_chrM_for_haplogroup.vcf.gz` | the chrM calls both tools read |

Both reports (step 24 and the text report) show the haplogroup and the contamination status, or "not checked" when the input was the step 03 VCF.

## Runtime
About a minute.

## Interpretation
- Common European haplogroups: H, U, J, T, K, V, W, X
- Output includes quality score (0-1): >0.9 = high confidence
- Discordant variants may indicate heteroplasmy (mixture of mtDNA types)
- Contamination `YES` with a level of a few percent: compare with step 33's VerifyBamID2 estimate, which reads the nuclear genome. The two can differ when the mtDNA copy number differs between the two people's cells.
