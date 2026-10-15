# Step 41: Blood Groups (predicted)

## What This Does

Predicts the red cell and platelet antigens the genome encodes: ABO (with A1 and A2), RhD, RhCE (C/c, E/e), Kell (K/k, Kp, Js), Duffy (with the GATA-box variant), Kidd, MNS (M/N, S/s), Lutheran, Diego, Dombrock, Colton, and the platelet antigens HPA-1, -2, -3, -4, -5 and -15. Opt-in, and the last optional step.

The result is a prediction from the DNA, not a serological typing. It is information about the sample, and the table states only what the variants predict.

## Why

Most blood group antigens differ by one known variant, so the genome predicts them. Nothing else in the pipeline reports them, and working out even ABO by hand from a VCF means knowing which reference allele GRCh38 carries (it carries the O allele's c.261delG, so an A or a B shows up as an insertion).

## Tool

- `bin/bloodtype_report.py` (standard library Python), applying the marker table `assets/bloodtype_markers.tsv`
- samtools `bedcov` for the read depth over the marker sites and RHD

The variants and the antigens they define are the ISBT allele definitions, as RBCeq2 2.4.4 packages them (database v2.5.1, MIT licence, Australian Red Cross Lifeblood). The table's header names the release and the checksum of the file they were read from. Each position and REF base was checked against Ensembl GRCh38 and its dbSNP id. The HPA positions come from dbSNP through Ensembl, because that RBCeq2 database lists them at their GRCh37 positions.

Why not run RBCeq2 or bloodTyper itself: RBCeq2 is on PyPI only, with no Bioconda package and no container image, and its README says it does not type RH from ordinary short-read calls (it needs long reads or a caller that writes copy numbers). bloodTyper is on neither PyPI nor Bioconda; the one Docker Hub image is a personal upload. A small table applied in the pinned Python image needs no new image and reads RhD from depth.

## Docker Image

- `PYTHON_IMAGE` for the rules, `SAMTOOLS_IMAGE` for the depth

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz`: one sample, left-aligned records (DeepVariant writes them so; split multi-allelic records are read as they are)
- `${SAMPLE}/aligned/${SAMPLE}_sorted.bam` with its `.bai` (`ALIGN_DIR` as for the other BAM steps). Optional: without it RhD is not determined, and a site the VCF has no record for is read as the reference allele at confidence medium.

## Command

```bash
./scripts/41-bloodtype.sh your_name
```

With `run-all.sh`: `TOOLS=...,bloodtype`. It is opt-in and needs no data, so a run lists step 41 only when `TOOLS` names it. In the Nextflow pipeline: `bloodtype` in `--tools` (processes `BLOODTYPE_DEPTH` and `BLOODTYPE`).

## How each system is read

| System | Variants | Rule |
|---|---|---|
| ABO | c.261delG, c.802G>A, c.803G>C, c.796C>A, c.1061delC | GRCh38 carries c.261delG (O.01): each insertion there is a functional allele. c.802G>A makes one of them O.02, c.803G>C makes one B, the rest are A; c.1061delC makes an A an A2. A1 with A2 reads A1 |
| RhD | read depth over RHD; c.487-19_504dup with c.609G>A | Copies = 2 x RHD depth / depth of the two flanks outside the RHD deletion (0 under 0.5, 1 under 1.5, else 2). The RHD pseudogene (RHD*08N.01) is an inactive copy. Exons 3 to 7 without RHD reads in a gene that is there is the pattern of an RHD-RHCE hybrid gene: RhD is then not determined |
| RhCE | c.307C>T (C/c), c.676G>C (E/e) | one allele per copy of the variant |
| Kell | c.578C>T (K/k), c.841C>T (Kp), c.1790T>C (Js) | the same |
| Duffy | c.125G>A (Fya/Fyb), c.-67T>C (GATA box), c.265C>T (FY*X) | The GATA box stops expression on red cells; with no phase it is placed on a Fyb chromosome first, where it usually lies (FY*02N.01). FY*X gives a weak Fyb (`b+w`) |
| Kidd, MNS, Lutheran, Diego, Dombrock, Colton | one variant each pair (c.838G>A; c.59T>C and c.143C>T; c.230G>A; c.2561C>T; c.793G>A; c.134C>T) | one allele per copy |
| HPA-1, -2, -3, -4, -5, -15 | ITGB3 c.176T>C, GP1BA c.482C>T, ITGA2B c.2621T>G, ITGB3 c.506G>A, ITGA2 c.1600G>A, CD109 c.2108C>A | one allele per copy; GRCh38 carries HPA-15b |

A site the VCF has no record for is the reference allele when the BAM has at least 8 reads there at MAPQ 20, and not determined below that. Without the BAM it is read as reference, and the confidence says so. A contig the BAM has no `@SQ` line for (a BAM cut to a few chromosomes) gets no depth rows, because `samtools bedcov` stops on one: the sites there with no VCF record are not determined, and RhD is not determined without chr1.

## Output

| File | Contents |
|---|---|
| `bloodtype/${SAMPLE}_bloodtype.tsv` | one row per system: `system`, `phenotype` (for example `O`, `D+`, `Fy(a+b-)`, `K-k+ Kp(a-b+) Js(a-b+)`), `genotype` (the alleles, for example `O.01/O.01`, `RHD/RHD deleted (RHD copies by read depth 0.91)`), `confidence`, `note` |
| `bloodtype/${SAMPLE}_bloodtype_evidence.tsv` | one row per site (status: called, covered, assumed, low depth or no call; alt copies; read depth) and per RHD region (mean depth and copies) |
| `bloodtype/${SAMPLE}_bloodtype_depth.tsv` | samtools `bedcov -Q 20 -j` over the marker table (with a BAM) |

Confidence:

- **high**: every site was called in the VCF or covered by reads
- **medium**: a site was read as reference with no BAM depth to confirm it
- **low**: a site in a region where short reads are unreliable, named in the note (RHCE exon 2 for C/c, which is RHD-like on a C allele; GYPA, GYPB and GYPE, which are near copies), markers that disagree, or another allele at a site
- **not determined**: a needed site had no call or too few reads; RhD without the BAM, with flanks under 10x, or with the hybrid pattern

Both reports show the table as their "Blood Groups (predicted)" section.

## Runtime

Not measured on a whole genome. The step reads the VCF once, line by line, and samtools reads the BAM at 41 small regions.

## Notes

- RhD here is the number of RHD copies, not their kind: weak D and partial D (RHD missense alleles) are not assessed, and neither are the RHD-CE-D hybrids that keep exons 3 to 7. The common RhD-negative allele, the deletion of the whole gene (RHD*01N.01), is what the depth sees.
- Phase is not read. ABO assumes A and B sit on different chromosomes, Duffy places the GATA box on Fyb first.
- Rare and null alleles outside the table (for example Kidd's Jk(a-b-) alleles, ABO's weak subgroups, cis-AB) are not assessed.
- On the e2e fixture the case reads HG002's reads at the marker sites from GIAB's 60x BAM and HG002's genotypes from GIAB's v4.2.1 benchmark, because the fixture's slices hold none of these genes. Those published calls give O (no c.261 insertion), and the reads give two RHD copies (D+). A copy of the same reads without any RHD read gives D-, and with half of them D+ with one copy deleted.
