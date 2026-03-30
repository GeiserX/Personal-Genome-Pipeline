# Recommended Resources

Free resources for understanding your genomic data, organized from most accessible to most technical.

---

## Start Here (No Background Needed)

### Understanding Genetics Basics

- **[Your Genome](https://www.yourgenome.org/)** — Wellcome Connecting Science. Excellent visual explainers on DNA, genes, inheritance, and genetic testing. Start here if you are new to genetics.
- **[MedlinePlus Genetics](https://medlineplus.gov/genetics/)** — NIH's plain-language genetics reference. Look up any gene or condition mentioned in your results.
- **[OMIM (Online Mendelian Inheritance in Man)](https://www.omim.org/)** — The definitive catalog of human genes and genetic disorders. Search for any gene from your ClinVar report to understand the associated conditions.

### Understanding Your Specific Results

- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** — Look up any variant ID (rsXXXXX) from your ClinVar screen to see its clinical significance, supporting evidence, and review status. Pay attention to the star rating and the number of submitters.
- **[gnomAD Browser](https://gnomad.broadinstitute.org/)** — Look up any variant to see how common it is in the general population. If a "pathogenic" variant has a gnomAD AF > 1%, the classification may be wrong.
- **[PharmGKB](https://www.pharmgkb.org/)** — The pharmacogenomics knowledge base. Look up any gene from your PharmCAT report to understand which drugs are affected and what the clinical recommendations are.
- **[CPIC Guidelines](https://cpicpgx.org/guidelines/)** — Clinical Pharmacogenetics Implementation Consortium. The clinical guidelines that PharmCAT's recommendations are based on. Organized by gene-drug pair.

---

## Intermediate (Some Biology Background Helpful)

### Variant Interpretation

- **[ClinGen](https://www.clinicalgenome.org/)** — The Clinical Genome Resource. Maintains gene-disease validity assessments and variant interpretation standards. Useful for understanding how variants are classified.
- **[ACMG SF v3.2 Gene List](https://www.nature.com/articles/s41436-023-02171-w)** — The 81 genes that ACMG recommends reporting incidental findings for. These are the genes CPSR screens. Understanding this list helps you prioritize findings.
- **[MitoMap](https://www.mitomap.org/)** — The human mitochondrial genome database. Look up any mitochondrial variant from step 20 to check its disease association and heteroplasmy significance.
- **[GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1116/)** — Expert-authored, peer-reviewed descriptions of genetic conditions. If you find a pathogenic variant, the GeneReviews entry for the associated condition is the best single resource for understanding it.

### Structural Variants and CNVs

- **[DECIPHER](https://www.deciphergenomics.org/)** — Database of genomic variation and associated phenotype. Useful for checking if a CNV or SV overlaps with known pathogenic regions.
- **[DGV (Database of Genomic Variants)](http://dgv.tcag.ca/)** — Catalog of structural variants found in healthy populations. If your SV is in DGV, it is likely benign.
- **[ClinGen Dosage Sensitivity](https://dosage.clinicalgenome.org/)** — Which genes are sensitive to having too few copies (haploinsufficiency) or too many copies (triplosensitivity). Essential for interpreting CNV results.

---

## Advanced (For Researchers and Power Users)

### Tools and Methods

- **[GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)** — The gold standard for variant calling pipelines. Useful for understanding the theory behind this pipeline's approach.
- **[DeepVariant paper (2018)](https://www.nature.com/articles/nbt.4235)** — The original publication describing the deep learning variant caller used in step 3.
- **[Manta paper (2016)](https://academic.oup.com/bioinformatics/article/32/8/1220/1743909)** — How Manta detects structural variants from paired-end reads.
- **[VEP documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)** — Complete reference for understanding VEP's annotation fields and consequence types.

### Population Genetics

- **[1000 Genomes Project](https://www.internationalgenome.org/)** — Catalog of human genetic variation from 2,504 individuals across 26 populations.
- **[gnomAD paper (2020)](https://www.nature.com/articles/s41586-020-2308-7)** — Understanding allele frequencies and how gnomAD's constraint metrics (pLI, LOEUF) predict gene essentiality.

### Online Courses (Free)

- **[Coursera: Genomic Data Science](https://www.coursera.org/specializations/genomic-data-science)** — Johns Hopkins. 8-course specialization covering the full analysis pipeline from FASTQ to interpretation.
- **[MIT OpenCourseWare: Computational Biology](https://ocw.mit.edu/courses/6-047-computational-biology-fall-2015/)** — MIT 6.047. More theoretical but excellent for understanding algorithms behind variant calling and annotation.

---

## Useful Command-Line References

### bcftools Cheat Sheet

```bash
# Count variants by type
bcftools stats file.vcf.gz | grep "^SN"

# Extract specific fields
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' file.vcf.gz

# Filter to PASS variants only
bcftools view -f PASS file.vcf.gz

# Filter to a specific region
bcftools view -r chr22:16000000-17000000 file.vcf.gz

# Filter to heterozygous variants only
bcftools view -g het file.vcf.gz

# Filter to homozygous ALT variants only
bcftools view -g hom file.vcf.gz

# Count variants per chromosome
bcftools view -f PASS file.vcf.gz | grep -v '^#' | cut -f1 | sort | uniq -c | sort -rn
```

### samtools Cheat Sheet

```bash
# Quick alignment summary
samtools flagstat file.bam

# Per-chromosome read counts
samtools idxstats file.bam

# Extract reads from a specific region
samtools view -b file.bam chr22:16000000-17000000 > region.bam

# View BAM header (reference info, read groups)
samtools view -H file.bam

# Calculate average depth
samtools depth -a file.bam | awk '{sum+=$3; n++} END {print sum/n}'
```

---

## Community

- **[r/bioinformatics](https://www.reddit.com/r/bioinformatics/)** — General bioinformatics discussion
- **[r/genomics](https://www.reddit.com/r/genomics/)** — Genomics-focused community
- **[BioStars](https://www.biostars.org/)** — Q&A forum for bioinformatics (like StackOverflow for genomics)
- **[SEQanswers](http://seqanswers.com/)** — Sequencing technology discussion forum
