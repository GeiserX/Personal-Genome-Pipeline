# Step 5: Structural Variant Annotation (AnnotSV)

## What This Does
Classifies every structural variant from Manta using ACMG guidelines (class 1-5), adding gene overlap, population frequency, and clinical significance.

## Why
Raw Manta output contains thousands of SVs with no clinical interpretation. AnnotSV tells you which ones matter by cross-referencing known pathogenic SVs, gene databases, and population data.

## Tool
- **AnnotSV** — ACMG-compliant structural variant annotation and classification

## Docker Image
```
getwilds/annotsv:3.4.4
```

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

docker run --rm \
  --cpus 4 --memory 8g \
  -v ${GENOME_DIR}:/genome \
  getwilds/annotsv:3.4.4 \
  AnnotSV \
    -SVinputFile /genome/${SAMPLE}/manta/results/variants/diploidSV.vcf.gz \
    -genomeBuild GRCh38 \
    -outputFile /genome/${SAMPLE}/annotsv/${SAMPLE}_annotsv \
    -outputDir /genome/${SAMPLE}/annotsv
```

## Output
- `${SAMPLE}_annotsv.tsv` — main annotated output (one row per SV, with ACMG class)
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
- AnnotSV bundles its own annotation databases inside the Docker image — no separate download needed
