# Look Up a Variant, an rsID or a Gene

## What This Does

Answers "what is my genotype at rs4988235?", "what do I have at chr2:135851076?" or "what do I carry in CFTR?" from the files the pipeline already wrote. You do not open a VCF by hand.

- **A site** (an rsID or a GRCh38 position): your genotype from the gVCF of [step 3](03-variant-calling.md), its read depth and genotype quality, and, when you carry a variant there, its gene, consequence, ClinVar significance and gnomAD frequency from the annotated VCF of [step 13](13-vep-annotation.md).
- **A gene**: every variant called in it with the same details, and how many of the gene's bases have a call at all.
- Both end with what the reports already say about the gene or the rsID: PharmCAT and CPIC, the ClinVar screen, the couple carrier screen, the traits panel and every other section of `summary.json` that names it.

## Why

Most questions about a genome are about one variant or one gene, often one a doctor, a paper or a relative named. The answer sits in two large compressed VCFs, and reading them right is easy to get wrong: a variant-only VCF has no line where you match the reference, so "no line" does not mean "reference". The gVCF does have those lines (reference blocks), and the lookup reads it that way.

## Tool

- **bcftools** (`bcftools view -r`) cuts the region out of the indexed gVCF and annotated VCF.
- **`bin/genome_lookup.py`** (Python standard library) finds the region and writes the answer.

## Docker Image

- `BCFTOOLS_IMAGE` and `PYTHON_IMAGE`, pinned in `versions.env`; [Image versions](versions.md) lists the tags. The command adds no image.

## Input

| File | Needed | What it gives |
|---|---|---|
| `${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz` (step 3) | yes (else the VCF) | the genotype at any position |
| `${SAMPLE}/vep/${SAMPLE}_annotated.vcf.gz` (step 30) or `_vep.vcf.gz` (step 13) | no | gene, consequence, HGVS, ClinVar, gnomAD of your variants |
| `${SAMPLE}/summary.json` (`generate-report.sh`) | no | what the reports say |
| `reference/consumer/23andme_v5_sites_dbsnp<build>.tsv.gz` (`setup.sh --consumer-sites`) | no | the position of 614,138 rsIDs, also where you match the reference |
| `reference/gencode.v<release>.basic.genes.gtf` (`setup.sh`) | for gene queries | the gene's GRCh38 span |

Without a gVCF the command reads the VCF and says so: a site where you match the reference then reads as a no call.

How an rsID is found, in this order: the consumer site table, then `assets/traits.tsv` and `assets/bloodtype_markers.tsv`, then the ID column and VEP's `Existing_variation` of the annotated VCF. The annotated VCF only holds variants you carry, so an rsID where you match the reference is found only through the site table. If none of them has it, look the rsID up on [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) and ask for its GRCh38 position instead.

## Command

```bash
GENOME_DIR=/path/to/data ./scripts/pgp-lookup.sh your_name rs4988235
GENOME_DIR=/path/to/data ./scripts/pgp-lookup.sh your_name chr2:135851076 CFTR   # several queries at once
```

On Nextflow output, which keeps the same per-sample folders, use the output directory as `GENOME_DIR`. When the reference folder is not inside it, point `REF_DIR` at it:

```bash
GENOME_DIR=./results REF_DIR=/path/to/data/reference ./scripts/pgp-lookup.sh your_name CFTR
```

A position is `chr7:117559590`, `7:117559590` or `chr7:117,559,590`; `chrM` and `chrMT` are the same. A gene is its HGNC symbol as GENCODE names it (`CFTR`).

## Output

The answer is printed; nothing is written beside your results. A made-up example, shortened (the impact, HGVSc and rsID lines are left out):

```text
== rs1001: chr7:250 (GRCh38; from the annotated VCF (a variant called in this sample))
  genotype: C/T (heterozygous; REF C, ALT T)
  evidence: depth 33, genotype quality 40, QUAL 42.5
  variant: chr7:250 C>T C/T (heterozygous)
    gene: GENEA
    consequence: missense_variant
    HGVSp: ENSP1:p.Arg4Cys
    ClinVar (pipeline's ClinVar file): Pathogenic
    ClinVar review status: criteria provided, multiple submitters
    gnomAD genomes AF: 0.0003
  in the reports:
  PharmCAT (step 07):
    - gene=GENEA, diplotype=*1/*2, phenotype=Intermediate
```

How to read the genotype line:

| Genotype line | Meaning |
|---|---|
| `A/A (homozygous reference, inside the reference block ...)` | you match the reference here; the base is shown when the site table or a pipeline table knows it, else `ref/ref` |
| `C/T (heterozygous ...)`, `T/T (homozygous ...)` | a variant call, the alleles in the order of the gVCF's GT field |
| `no call (no gVCF record covers this position ...)` | no reads here, or the caller made no call: unknown, not the reference |
| `no call (reference block ... has no genotype)` | the caller saw the position but could not call it |
| `this position lies inside the variant called at ...` | the position is part of a longer variant (a deletion) that starts before it |
| `filter: RefCall` (or another filter) | the caller did not pass the call; read it as uncertain |

For a gene, the line `N of M bases have a genotype call` tells you how much of the gene was readable. A variant can only be missing from the list where there is no call.

Two ClinVar lines can appear. `ClinVar (pipeline's ClinVar file)` comes from the ClinVar release step 13 annotated with, the one [step 6](06-clinvar-screen.md) screens against. `clinical significance (VEP cache)` is Ensembl's copy inside the VEP cache, fixed with the cache release, so it can be older than the first.

## Runtime

A site or a gene takes a few seconds with indexed files: bcftools reads only the region, and a container starts for each pass. An rsID that is in neither the site table nor the pipeline tables makes the command read the whole annotated VCF once, which takes minutes on a whole genome. A file without a `.tbi` index is streamed whole, with a note: the answer is the same and only the records in the region are kept in memory, but on a whole genome it takes minutes. `bcftools index -t` makes it fast.

## Notes

- Positions are GRCh38, the pipeline's reference. A position from a GRCh37 (hg19) source points elsewhere: convert it first, or ask by rsID.
- The command reads the gVCF's first sample. It changes nothing in your results.
- The reports part reads `summary.json` as it is. After rerunning a step, rerun `generate-report.sh` so the lookup sees the new result. Sections that did not run are not read; a section from an earlier run is marked stale.
- An answer here is a reading of your own files, not a clinical result. A finding that matters goes to a genetic counsellor or physician and is confirmed by a clinical test.
