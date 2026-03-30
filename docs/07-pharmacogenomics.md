# Step 7: Pharmacogenomics (PharmCAT)

## What This Does
Clinical-grade pharmacogenomic analysis — determines how you metabolize drugs based on your DNA. The same tool hospitals use for precision medicine.

## Why
Identifies which drugs work well, which need dose adjustments, and which to avoid entirely. Covers 23 pharmacogenes affecting hundreds of medications.

## Tool
- **PharmCAT** v2.15.5 (Pharmacogenomics Clinical Annotation Tool, CPIC/PharmGKB)

## Docker Image
```
pgkb/pharmcat:2.15.5
```

## Command
```bash
SAMPLE=your_sample
GENOME_DIR=/path/to/your/data

# PharmCAT needs the reference genome for VCF preprocessing
docker run --rm \
  --cpus 2 --memory 4g \
  -v ${GENOME_DIR}/${SAMPLE}/vcf:/data \
  -v ${GENOME_DIR}/reference:/ref \
  pgkb/pharmcat:2.15.5 \
  java -jar /pharmcat/pharmcat.jar \
    -vcf /data/${SAMPLE}.vcf.gz \
    -refFasta /ref/Homo_sapiens_assembly38.fasta \
    -o /data/ \
    -bf ${SAMPLE}

# Output: ${SAMPLE}.report.html (interactive HTML report)
```

## Output
- HTML report with drug recommendations per gene
- Covers CYP2C19, CYP2D6, CYP2B6, CYP3A4/5, UGT1A1, DPYD, NAT2, TPMT, etc.
- Star allele calls with metabolizer status (Poor/Intermediate/Normal/Rapid/Ultra-rapid)

## Key Genes
| Gene | Drugs Affected | Example |
|---|---|---|
| CYP2C19 | SSRIs, PPIs, clopidogrel | *1/*17 = rapid → SSRIs fail faster |
| CYP2D6 | 25% of all drugs, opioids, tamoxifen | Complex — may need BAM-based calling |
| UGT1A1 | Irinotecan, bilirubin clearance | *28/*28 = Gilbert's syndrome |
| DPYD | 5-FU, capecitabine (chemo) | Poor = lethal toxicity |
| NAT2 | Isoniazid, hydralazine | Slow acetylator = increased toxicity |

## Limitations
- **CYP2D6** often returns `Not called` — gene has pseudogene homology that confounds VCF-based calling. Use Cyrius or StellarPGx (BAM-based) if CYP2D6 is critical.
- PharmCAT may disagree with lab reports on complex haplotypes (e.g., NAT2). When in doubt, trust PharmCAT + raw VCF over lab transcription.
