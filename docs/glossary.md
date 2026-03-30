# Glossary

Terms used throughout this pipeline's documentation, explained for people who are not bioinformaticians.

---

**Allele** — One of the possible versions of a gene or variant at a specific position. You have two alleles at most positions (one from each parent). Example: at a SNP position, the reference allele might be `A` and the alternate allele `G`.

**Allele Frequency (AF)** — How common a variant is in a population. AF=0.01 means 1% of people carry it. Common in gnomAD and ClinVar annotations.

**ALT** — The alternate allele; the version that differs from the reference genome. In a VCF file, `REF=A, ALT=G` means the reference has `A` but this sample has `G`.

**Autosomal Dominant** — A variant in one copy of a gene is enough to cause the condition. Contrast with autosomal recessive.

**Autosomal Recessive** — Both copies of a gene must be affected (homozygous or compound heterozygous) for the condition to manifest. Having one copy makes you a carrier.

**BAM** — Binary Alignment Map. A binary file containing sequencing reads aligned to a reference genome. Typically 80-120 GB for 30X WGS. Requires a `.bai` index file.

**Benign** — A variant classification meaning the variant is not disease-causing. The opposite of pathogenic.

**BND** — Breakend. A structural variant type where read pairs map to unexpected locations (different chromosomes or distant positions). Most are artifacts.

**Carrier** — A person who has one copy of a recessive disease variant (heterozygous). Carriers are typically unaffected but can pass the variant to children.

**chrM** — Chromosome M (mitochondrial). The small circular genome inside mitochondria (~16,569 bp). Inherited maternally.

**ClinVar** — NCBI's public database of clinically significant genetic variants. Each entry includes a significance classification (pathogenic, benign, VUS, etc.) and supporting evidence.

**CNV** — Copy Number Variant. A region of the genome that is deleted or duplicated relative to the reference. Detected by CNVnator (depth-based) and Manta/Delly (read-pair based).

**Compound Heterozygous** — Having two different pathogenic variants in the same gene, one from each parent. Can cause autosomal recessive disease even though neither variant is homozygous.

**Coverage / Depth** — The average number of sequencing reads covering each position. "30X" means each position is read ~30 times on average. Higher coverage = more accurate variant calling.

**CRAM** — Compressed Reference Alignment Map. A compressed version of BAM (40-60% smaller) that stores reads as differences from the reference genome.

**De Novo** — A variant that appeared for the first time in an individual (not inherited from either parent). Typically 50-100 de novo SNVs per genome.

**DEL** — Deletion. A structural variant where a segment of DNA is missing compared to the reference.

**DeepVariant** — Google's deep learning-based variant caller. Uses a neural network trained on labeled sequencing data to call SNPs and indels.

**DUP** — Duplication. A structural variant where a segment of DNA is copied one or more extra times.

**FASTQ** — A text-based format for storing raw sequencing reads and their quality scores. Paired-end sequencing produces two files (R1 and R2).

**gnomAD** — Genome Aggregation Database. Contains allele frequencies from ~140,000 exomes and ~76,000 genomes. The standard reference for how common a variant is in healthy populations.

**GRCh37 / hg19** — The older human genome reference assembly (2009). Some vendors still deliver data aligned to this build.

**GRCh38 / hg38** — The current human genome reference assembly (2013, with updates). This pipeline uses GRCh38 exclusively.

**Haplogroup** — A group of similar genetic variants inherited together, used to trace ancestry through maternal (mitochondrial) or paternal (Y-chromosome) lineages.

**Heteroplasmy** — Having a mixture of different mitochondrial DNA sequences within the same cell. Low-level heteroplasmy (<5%) is normal and increases with age.

**Heterozygous (het, 0/1)** — Having two different alleles at a position (one reference, one alternate). Written as `0/1` in VCF genotype fields.

**HLA** — Human Leukocyte Antigen. A highly polymorphic gene complex on chromosome 6 involved in immune function. Difficult to type accurately from short-read WGS.

**Homozygous (hom, 1/1)** — Having two copies of the same alternate allele at a position. Written as `1/1` in VCF genotype fields.

**Imputation** — Statistical prediction of missing genotypes based on known patterns of genetic variation (linkage disequilibrium). Used to fill gaps in genotyping array data.

**Indel** — An insertion or deletion of one or more bases. Distinct from structural variants (which are larger, typically >50 bp).

**INV** — Inversion. A structural variant where a segment of DNA is flipped in orientation relative to the reference.

**Likely Pathogenic** — A variant classification one step below "Pathogenic." Strong evidence of disease association but not yet definitive.

**Loss of Function (LoF)** — A variant that is predicted to destroy the function of a gene (stop-gain, frameshift, splice donor/acceptor disruption). Also called "HIGH impact" by VEP.

**MAF** — Minor Allele Frequency. The frequency of the less common allele at a variant position. Often used interchangeably with AF for rare variants.

**Manta** — An SV caller that uses paired-end read distance and split reads to detect structural variants. Fast and accurate for deletions and duplications.

**Missense** — A single nucleotide change that results in a different amino acid in the protein. May or may not affect protein function.

**PASS** — A VCF filter status indicating the variant passed all quality filters. Non-PASS variants have quality concerns and should generally be excluded from analysis.

**Pathogenic** — A variant classification meaning the variant is known to cause disease. The strongest clinical significance level in ClinVar.

**Penetrance** — The proportion of people with a pathogenic variant who actually develop the associated disease. "Incomplete penetrance" means not everyone with the variant gets sick.

**PGx / Pharmacogenomics** — The study of how genetic variants affect drug metabolism and response. Used to guide medication selection and dosing.

**PolyPhen** — A tool that predicts whether an amino acid change is damaging to protein function. Scores range from 0 (benign) to 1 (probably damaging).

**REF** — The reference allele; the version found in the reference genome at that position.

**ROH** — Run of Homozygosity. A long stretch of DNA where both copies are identical. Large ROH segments can indicate parental relatedness.

**SIFT** — A tool that predicts whether an amino acid change affects protein function. Scores < 0.05 are predicted "deleterious"; > 0.05 are "tolerated."

**SNP / SNV** — Single Nucleotide Polymorphism / Variant. A single base-pair change (e.g., A>G). The most common type of genetic variant (~4-5 million per genome).

**Star Allele** — A standardized naming system for pharmacogene variants (e.g., CYP2C19\*2, CYP2D6\*4). Used by PharmCAT and pharmacogenomics databases.

**STR** — Short Tandem Repeat. A DNA sequence of 2-6 bases repeated many times in a row (e.g., CAGCAGCAG...). Pathogenic expansions cause diseases like Huntington's.

**SV** — Structural Variant. A large-scale genomic rearrangement (>50 bp): deletion, duplication, inversion, insertion, or translocation.

**Tabix / TBI** — An index format for coordinate-sorted, block-compressed files (VCF, BED). Allows rapid random access to specific genomic regions. The `.tbi` file must accompany `.vcf.gz` files.

**VCF** — Variant Call Format. A tab-delimited text file listing genetic variants with their position, alleles, quality scores, and annotations. The standard format for variant data.

**VEP** — Variant Effect Predictor. Ensembl's tool for annotating variants with functional consequences, population frequencies, and pathogenicity predictions.

**VUS** — Variant of Uncertain Significance. Not enough evidence to classify as pathogenic or benign. Not clinically actionable. May be reclassified in the future as more data accumulates.

**WGS** — Whole Genome Sequencing. Sequencing the entire ~3.1 billion base pairs of a genome. Produces FASTQ files that are processed through this pipeline.

**X-linked** — A gene located on the X chromosome. Males (XY) have only one copy, so a single pathogenic variant can cause disease. Females (XX) have two copies and may be carriers.
