# Step 5: Structural Variant Annotation (AnnotSV)

## What This Does
Classifies every structural variant from Manta using ACMG guidelines (class 1-5), adding gene overlap, population frequency, and clinical significance.

## Why
Raw Manta output contains thousands of SVs with no clinical interpretation. AnnotSV tells you which ones matter by cross-referencing known pathogenic SVs, gene databases, and population data.

## Tool
- **AnnotSV** — ACMG-compliant structural variant annotation and classification

## Docker Image
```
quay.io/biocontainers/annotsv:3.5.10--hdfd78af_0
```

## Annotation Data
The image holds AnnotSV's code only. Its annotation data (genes, known pathogenic SVs, population frequencies) is a separate 5.3 GB download, and AnnotSV exits with an error without it. `./scripts/setup.sh` downloads and unpacks it into `${GENOME_DIR}/annotsv_annotations/`. To do it by hand:

```bash
curl -fL -C - -o ${GENOME_DIR}/Annotations_Human_3.5.tar.gz \
  https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
mkdir -p ${GENOME_DIR}/annotsv_annotations
tar -xzf ${GENOME_DIR}/Annotations_Human_3.5.tar.gz -C ${GENOME_DIR}/annotsv_annotations
```

The script checks for `${GENOME_DIR}/annotsv_annotations/Annotations_Human/Genes/GRCh38` and stops with a pointer to `setup.sh` when it is missing. `run-all.sh` reports step 5 as skipped in that case.

## Command
```bash
./scripts/05-annotsv.sh your_sample
```

What the script runs:

```bash
docker run --rm --user root \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  quay.io/biocontainers/annotsv:3.5.10--hdfd78af_0 \
  AnnotSV \
    -SVinputFile /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
    -outputFile /genome/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv \
    -genomeBuild GRCh38 \
    -annotationMode both \
    -annotationsDir /genome/annotsv_annotations
```

To annotate another SV VCF (for example Sniffles2 output), set `SV_VCF` to its host path; the file must be inside `${GENOME_DIR}`.

## Output
- `${GENOME_DIR}/${SAMPLE}/annotsv/${SAMPLE}_sv_annotated.tsv` — main annotated output (one row per SV, with ACMG class)
- Columns include: SV type, coordinates, overlapping genes, DGV frequency, ACMG classification, ClinVar hits

## ACMG Classification
| Class | Meaning | Action |
|---|---|---|
| 1 | Benign | Ignore |
| 2 | Likely benign | Ignore |
| 3 | Variant of uncertain significance (VUS) | Review if in known disease gene |
| 4 | Likely pathogenic | Investigate — check gene, inheritance, phenotype |
| 5 | Pathogenic | Investigate — known disease-causing SV |

## Important Notes
- **Class 4-5 = pathogenic/likely pathogenic** — these require manual review
- SVs >5MB in short-read WGS are usually artifacts from segmental duplications — do not trust large calls blindly
- Most SVs will be class 2-3 (benign/VUS) — this is normal for a healthy genome
- Input must be from Manta step 4 (`diploidSV.vcf.gz`), not the unfiltered candidates
- The annotation databases are not in the Docker image: they live in `${GENOME_DIR}/annotsv_annotations/` (see Annotation Data above)
