# Step 49: Clonal Haematopoiesis (CHIP)

## What This Does

Looks for mutations that blood stem cells acquired during life and that now sit in a clone of the blood: clonal haematopoiesis of indeterminate potential (CHIP). It calls low-VAF changes in 13 CHIP genes (DNMT3A, TET2, ASXL1, JAK2, TP53, PPM1D, SF3B1, SRSF2, U2AF1, CBL, GNB1, IDH1, IDH2), keeps the ones in each gene's published classes, and says how small a clone it could have found at the sample's depth. Opt-in.

This is not a diagnosis. CHIP becomes common with age, and most people who have it never get a blood cancer.

## Why

A whole genome from blood holds the blood's acquired mutations as well as the inherited ones. Jaiswal et al. found such mutations in 9.5% of people aged 70 to 79 in exome data of 17,182 people, most of them in DNMT3A, TET2 and ASXL1, and they were linked to a higher risk of blood cancer, of coronary heart disease and of death (PMID 25426837; Genovese et al., PMID 25426838, found the same). Bick et al. found CHIP in 4,229 of 97,691 whole genomes of the TOPMed programme (PMID 33057201). The term CHIP and its threshold, a VAF of at least 2% in a gene recurrently mutated in blood cancers, are from Steensma et al. (PMID 25931582).

The germline caller of step 3 is built for 50% and 100% alleles and does not report a change carried by a few reads. Step 29 runs Mutect2 on these genes but leaves its calls unread. This step reads them.

## Tool

- GATK Mutect2, tumour-only, then LearnReadOrientationModel and FilterMutectCalls
- bcftools csq, which names each change on the gene's MANE Select transcript, and bcftools mpileup for the depth
- `bin/chip_report.py` (standard library Python), applying the gene table `assets/chip_genes.tsv`

## Docker Image

- `GATK_IMAGE`, `BCFTOOLS_IMAGE` and `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps)
- the reference with its `.fai` and `.dict`
- Ensembl's GRCh38 gene annotation, `reference/ensembl/Homo_sapiens.GRCh38.<release>.chr.gff3.gz` (`ENSEMBL_GFF_RELEASE` in `versions.env`, about 110 MB). Install it once with `./scripts/setup.sh --chip-data "${GENOME_DIR}"`; it is checked against `ENSEMBL_GFF_SHA256`.
- Optional, and the same files as [step 29](29-mutect2-somatic.md): gnomAD's AF-only germline resource `somatic/af-only-gnomad.hg38.vcf.gz` and a panel of normals `somatic/1000g_pon.hg38.vcf.gz`, each with its `.tbi`. Without gnomAD, inherited variants are left out by their VAF only.

## Command

```bash
./scripts/setup.sh --chip-data "${GENOME_DIR}"   # once
./scripts/49-chip.sh your_name
```

With `run-all.sh`: `TOOLS=...,chip`. It is opt-in, so a run lists step 49 only when `TOOLS` names it, and it is skipped without a BAM or without the annotation. In the Nextflow pipeline: `chip` in `--tools` with `--chip_gff` (and optionally `--chip_germline_resource` and `--chip_pon`), processes `CHIP_PREPARE`, `CHIP_MUTECT2`, `CHIP_ANNOTATE` and `CHIP_REPORT`. `CHIP_GENES=FILE` runs the script with another table in the same layout.

## How it decides

1. For each gene of `assets/chip_genes.tsv` the step takes the transcript Ensembl tags `MANE_Select`, and calls inside its coding exons, 10 bases wider on each side so the splice donor and acceptor sites are in. That is 162 exons and about 35 kb with Ensembl 116.
2. Mutect2 calls in tumour-only mode, with the read orientation model, and with gnomAD and the panel of normals when they are installed.
3. Each call that changes the protein is put in the gene's class from the table:

| Gene | Counted | Source (PMID) |
|---|---|---|
| DNMT3A | truncating anywhere; missense at R882 | 25426837, 21067377 |
| TET2 | truncating anywhere | 25426837, 19474426 |
| ASXL1 | truncating anywhere | 25426837, 19388938 |
| JAK2 | V617F | 25426837, 15858187 |
| TP53 | truncating anywhere; missense in the DNA-binding core domain, residues 102 to 292 | 25426837, 8023157 |
| PPM1D | truncating in the last exon (6) | 23242139, 29954749 |
| SF3B1 | missense at K700 | 21995386, 23335386 |
| SRSF2 | missense at P95 | 21909114, 23335386 |
| U2AF1 | missense at S34 and Q157 | 22158538, 23335386 |
| CBL | missense in exons 8 and 9 (codons 366 to 477, the linker and RING finger) | 19620960, 19276253 |
| GNB1 | missense at K57 | 25485910 |
| IDH1 | missense at R132 | 19657110, 26016821 |
| IDH2 | missense at R140 and R172 | 20171147, 26016821 |

"Truncating" is a frameshift, a stop gained, or a change of a splice donor or acceptor site. Every codon and exon is numbered on the MANE Select transcript; the residue at each hotspot was checked against Ensembl's protein sequence of that transcript. These are the best-known classes only, narrower than the variant lists of the large CHIP studies: a change outside them is listed as `other_change` and not counted.

4. Each change then gets a status:

| Status | Rule |
|---|---|
| `chip` | passed FilterMutectCalls, in the gene's class, VAF under 35%, at least 3 reads with the change |
| `likely_germline` | VAF 35% or more. An inherited variant sits near 50% or 100%, and without a matched normal sample the two cannot be told apart |
| `filtered` | FilterMutectCalls did not pass it; its filters are in the `note` column |
| `low_support` | fewer than 3 reads with the change |
| `other_change` | a protein change outside the gene's class |

Synonymous, intronic and UTR changes are left out.

## The detection limit

The step counts the depth over the same exons with bcftools mpileup at Mutect2's default read filters (mapping quality 20 and up, base quality 10 and up, duplicates left out), and gives each gene a detection limit: the smallest VAF at which a clone shows at least 3 reads with the change 8 times in 10 at that gene's mean depth (a binomial count). At 20x that is about 20%, at 30x about 14%, at 40x about 10%, at 60x about 7% and at 100x about 4%.

So at the depth of a usual 30x genome, a clone has to hold about 14% of the alleles to be found reliably. CHIP starts at 2%, so a clone between 2% and the limit is likely missed: a result with no call does not rule out CHIP. A gene whose depth is far below the genome's (reads that map with low quality there) gets a higher limit, or `too_shallow` when it has fewer than 3 reads.

## Output

| File | Contents |
|---|---|
| `chip/${SAMPLE}_chip.tsv` | one row per protein change in a CHIP gene: `gene`, `chrom`, `pos`, `ref`, `alt`, `consequence`, `protein_change` (for example `R882H`, `Q1200fs`), `exon`, `vaf`, `alt_reads`, `depth`, `filter`, `status`, `class` (truncating or missense, for a change in the gene's class), `note` |
| `chip/${SAMPLE}_chip_depth.tsv` | one row per gene: `transcript`, `chrom`, `coding_bases`, `mean_depth`, `detection_limit_vaf`, `status` (covered, too_shallow, not_on_reference) |
| `chip/${SAMPLE}_chip_filtered.vcf.gz` | Mutect2's calls after FilterMutectCalls |
| `chip/${SAMPLE}_chip.gff3`, `chip/${SAMPLE}_chip.bed` | the transcripts and exons it called on (the bash step) |

Both reports show a "Clonal Haematopoiesis (CHIP)" section: the calls, the detection limit at the sample's depth, the genes too shallow to read, and how many other changes were listed.

## Runtime

Not measured on a whole genome. Mutect2 reads about 35 kb of the BAM, which takes minutes, not hours.

## Notes

- A `chip` row is a finding worth showing a haematologist, together with a blood count. It does not mean a blood cancer, and the risks the studies report are averages over many people.
- Tumour-only calling cannot prove a change is acquired. The VAF rule removes most inherited variants, but a change at 35% or more is never counted, and a real large clone there is missed too.
- Steps 29 and 49 call the same genes in two ways: step 29 calls gene bodies and leaves the calls unread, step 49 calls the coding exons and reads them.
