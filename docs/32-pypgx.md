# Step 32: Comprehensive Pharmacogenomics (pypgx)

## What This Does

Comprehensive pharmacogenomic star allele calling for 23 clinically actionable genes, including structural variation (SV) detection from BAM read depth. Complements PharmCAT (step 7) by covering genes that VCF-only callers miss entirely, including gene deletions and duplications. See Limitations before trusting a CYP2D6 copy-number call.

## Why

PharmCAT (step 7) calls ~23 genes from VCF data alone. This works well for simple SNP-based star alleles but fails for genes with structural variation — CYP2D6, CYP2A6, GSTM1, and GSTT1 all have common whole-gene deletions and duplications that VCF callers cannot represent. CYP2D6 alone affects 25% of all prescribed drugs, and PharmCAT frequently returns "Not called" for it.

Cyrius (step 21) was designed specifically for CYP2D6 but can fail on some WGS samples due to CYP2D7 pseudogene homology. pypgx calls copy number from read depth, a different method; neither tool settles CYP2D6 alone (see Limitations).

pypgx also calls genes absent from PharmCAT entirely: COMT, MTHFR, ABCB1, GSTM1, GSTT1, and IFNL3.

## Tool

- **pypgx** v0.26.0 (Sboner Lab, Weill Cornell Medicine)
- License: Apache-2.0 (GPL-3.0 compatible)
- Publication: [Lee et al., 2019](https://doi.org/10.1002/cpt.1552)

## Docker Image

- `PYPGX_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Prerequisites

### pypgx-bundle (required, one-time download)

pypgx requires a companion data bundle containing 1000 Genomes phasing panels (for Beagle haplotype estimation) and CNV classifier models. The bundle is **not included** in the Docker image and must be downloaded separately (~370 MB):

```bash
cd ${GENOME_DIR}/reference
git clone --branch 0.26.0 --depth 1 https://github.com/sbslee/pypgx-bundle.git
```

The bundle version must match the pypgx version (`0.26.0`). The script validates the bundle exists at `${GENOME_DIR}/reference/pypgx-bundle/` and exits with download instructions if missing.

### Input files

- BAM from alignment (step 2): `${GENOME_DIR}/${SAMPLE}/aligned/${SAMPLE}_sorted.bam` (+ `.bai` index)
- VCF from variant calling (step 3): `${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}.vcf.gz` (+ `.tbi` index)

## Command

```bash
./scripts/32-pypgx.sh your_name
```

## Gene List

The script calls 23 curated genes spanning CPIC Level A/B recommendations and key genes PharmCAT misses.

### BAM-based genes (structural variation detection)

These genes have common whole-gene deletions, duplications, or hybrid alleles that cannot be detected from VCF data alone. pypgx uses BAM read depth to identify copy number changes.

| Gene | Why BAM-based | Clinical Impact |
|---|---|---|
| CYP2D6 | Deletions (*5), duplications (*1x2, *2x2), CYP2D7 hybrids | 25% of all drugs; codeine, tamoxifen, SSRIs |
| CYP2A6 | Whole-gene deletion (*4), duplications | Nicotine metabolism, tegafur |
| GSTM1 | Homozygous deletion (null genotype) in ~50% of population | Detoxification, carcinogen metabolism |
| GSTT1 | Homozygous deletion (null genotype) in ~20% of population | Detoxification, drug conjugation |

### VCF-based genes (SNP/indel star alleles)

These are called using both BAM and VCF data. Star alleles are determined from variant calls.

| Gene | Key Drugs | CPIC Level | Notes |
|---|---|---|---|
| CYP1A2 | Caffeine, clozapine, theophylline | B | |
| CYP2B6 | Efavirenz, methadone | A | |
| CYP2C9 | Warfarin, phenytoin, NSAIDs | A | |
| CYP2C19 | Clopidogrel, SSRIs, PPIs | A | |
| CYP3A4 | Tacrolimus (with CYP3A5) | B | |
| CYP3A5 | Tacrolimus | A | |
| CYP4F2 | Warfarin (vitamin K cycle) | A | |
| DPYD | Fluorouracil, capecitabine | A | |
| TPMT | Azathioprine, mercaptopurine | A | |
| NUDT15 | Azathioprine, mercaptopurine | A | |
| UGT1A1 | Irinotecan, atazanavir, bilirubin | A | |
| SLCO1B1 | Simvastatin, statins | A | |
| VKORC1 | Warfarin | A | |
| NAT2 | Hydralazine; isoniazid | A (hydralazine), C (isoniazid) | |
| COMT | Catecholamine metabolism | -- | Not in PharmCAT |
| MTHFR | Folate metabolism, methotrexate | -- | Not in PharmCAT |
| ABCB1 | Drug efflux (broad substrate range) | -- | Not in PharmCAT |
| G6PD | Rasburicase, primaquine, dapsone | A | |
| IFNL3 | Peginterferon (historical, DAAs replaced) | A | |

## CYP2D6 Structural Variation Detection

CYP2D6 is the most complex pharmacogene. It has a tandemly duplicated pseudogene (CYP2D7) with >90% sequence identity, and common structural variants in the general population:

- **Gene deletion (*5)**: Entire CYP2D6 removed. Homozygous = poor metabolizer.
- **Gene duplication (*1x2, *2x2, etc.)**: Extra functional copies. Can produce ultra-rapid metabolizer status.
- **CYP2D6/CYP2D7 hybrids (*36, *13, etc.)**: Recombination between gene and pseudogene.

pypgx detects these by analyzing read depth across the CYP2D6/CYP2D7 locus. A drop in coverage indicates deletion; elevated coverage indicates duplication. VCF data alone cannot represent these copy number changes, which is why PharmCAT returns "Not called" for CYP2D6 in most WGS samples.

The BAM-based calling uses both `--variants` and `--depth-of-coverage` per the upstream pypgx WGS workflow, combining SNV/haplotype information with read-depth SV detection for the most complete genotype call.

## CYP2D6 depth check

A copy-number call from depth is only as good as the depth. Before pypgx runs, the step measures, with mosdepth (`${MOSDEPTH_IMAGE}`), the mean depth over CYP2D6 (chr22:42,123,193-42,132,032, the region Cyrius uses) and over two 50 kb flanks outside the CYP2D6-CYP2D8 cluster (chr22:42.05-42.10 Mb and 42.20-42.25 Mb), once for all reads and once for reads with MAPQ >= 1 (callers ignore MAPQ 0 reads). `bin/cyp2d6_depth_check.py` then judges it:

- **unreliable** when more than 40% of the reads at CYP2D6 (if it has at least 15% of the flank depth) or in the flanks have MAPQ 0: the aligner placed them on more than one sequence. The step writes `CYP2D6 copy number unreliable: reads are multi-mapped (was this BAM aligned to a reference with ALT contigs?)` to its log and to `${SAMPLE}_cyp2d6_depth_check.tsv`, and the summary's CYP2D6 row says `Indeterminate` (pypgx's own call stays in `CYP2D6/results.zip`).
- **ok** otherwise. A real deletion leaves few reads at CYP2D6, but those map uniquely and the flanks keep their depth, so it passes.

Beside the verdict, the check gives three more figures, which do not change it. Cyrius's region includes REP6 and the 3' flank of the gene, where a sizeable share of the depth is MAPQ 0 even on a clean BAM, so its MAPQ >= 1 share reads worse than the gene is. `cyp2d6_depth_check.py measure` counts, from the BAM itself (pysam, in the pypgx image), the depth over that region and over the gene body alone, CYP2D6's MANE Select transcript (chr22:42,126,499-42,130,810, without REP6), for all reads, MAPQ >= 1 and MAPQ >= 20, leaving out the reads mosdepth leaves out (unmapped, secondary, QC-failed, duplicate). The check file then has `gene_mapq1_share` (the share the verdict reads), `gene_mapq20_share`, `body_depth_all`, `body_mapq1_share` and `body_mapq20_share`. A share reads `not measured` when the measurement did not run, and `no reads` when there was no depth to divide by. The two shares on real 30x genomes are still to be measured; the verdict rule stays the one above until they are.

Why the share of uniquely mapped reads, not the ratio of CYP2D6 to flank depth: on the e2e fixture's CYP2D slice mapped to the default no-ALT reference, CYP2D6 has 0.62 of the flank depth with all reads and 0.57 at MAPQ >= 1, its normal value; mapped to the Broad hg38 FASTA with ALT contigs, a third of the reads go to `chr22_KI270879v1_alt` and CYP2D6 keeps 0.28 of the flank depth with all reads. A ratio cut at 0.6 flags the clean BAM and passes the broken one. The share of MAPQ >= 1 reads at CYP2D6 is 0.91 without ALT contigs and 0.005 with them. The ALT depth A/B workflow (`scripts/ci/alt-depth-ab.sh`) runs the check on both mappings and fails unless it passes the first and flags the second.

[Step 36](36-pgx-consensus.md) reads the check: with `unreliable`, no CYP2D6 call reaches PharmCAT, whatever the callers say. The summary is written as `${SAMPLE}_pypgx_summary.tsv.partial` and renamed only after the check's verdict is applied, so a run that stops in between leaves no unchecked CYP2D6 call.

## Output

All output is written to `${GENOME_DIR}/${SAMPLE}/pypgx/`.

| File | Contents |
|---|---|
| `<gene>/results.zip` | Per-gene pypgx archive with genotype data |
| `${SAMPLE}_pypgx_summary.tsv` | Consolidated: gene, diplotype, phenotype, pypgx's copy-number call, source |
| `${SAMPLE}_cyp2d6_depth_check.tsv` | The CYP2D6 depth check: status (`ok` or `unreliable`), its message, the four depths and the MAPQ shares (see the depth check above) |
| `cyp2d6_depth/` | mosdepth's region files the check read, and `mapq_shares.tsv`, the depths `measure` counted |
| `${SAMPLE}_pharmcat_comparison.tsv` | Side-by-side comparison with PharmCAT. Step 27 (CPIC lookup) writes it here when steps 7 and 32 have both run; rerunning step 32 alone keeps it |

### Summary TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| Diplotype | Star allele call (e.g., \*1/\*4); `Indeterminate` for CYP2D6 when the depth check found multi-mapped reads; `FAILED` for a gene pypgx wrote no result for; `N/A` for a gene not callable on this reference (see Status) |
| Phenotype | Metabolizer status (e.g., Intermediate Metabolizer). `no pypgx phenotype rule` for a gene pypgx has no phenotype rule for (no `PhenotypeMethod` in its own gene table; in 0.26.0 among the genes run here: CYP1A2, CYP2A6, CYP3A4, CYP4F2, VKORC1, NAT2, COMT, MTHFR, ABCB1, GSTM1, GSTT1), where pypgx itself writes `Indeterminate` |
| CNV_call | The copy-number call pypgx itself made from read depth (its `CNV` field, for example `Normal` or `WholeDel1`). BAM-based genes only; `N/A` for VCF-based genes, and `N/A (CYP2D6 depth check)` for CYP2D6 when the depth check failed. Earlier versions guessed a Yes/No flag from the allele names, which flagged any name containing `*5` and missed the GSTM1/GSTT1/CYP2A6 whole-gene deletions |
| Source | BAM (SV genes) or VCF (variant-based genes) |
| Status | `called`, `no call` (pypgx ran and gave no diplotype), `no call (CYP2D6 depth check)`, `failed`, or `not callable on this reference` (GSTT1 on a reference without ALT contigs, see Limitations) |
| Defining_variants | For a called gene: how many of the variants that define its star alleles (pypgx's core variants for each allele, `pypgx.list_variants`) the VCF has, as `2/2 in the VCF`; `none` when the alleles are defined by no variant (\*1, a whole-gene deletion) |
| Min_GQ | The lowest GQ among those variants in the VCF |
| Weakest_variant | That variant with its GQ and AD, for example `chr10:94781859 G>A GQ 12 AD 26,3` |
| Defining_filter | `PASS` when every defining variant is in the VCF with FILTER `PASS`; otherwise the ones that are not PASS (`RefCall`, `LowQual`...) or not in the VCF |
| Quality | `low-confidence: check the reads` when a defining variant has a GQ below the floor, `ok` otherwise, `N/A` without defining variants |

### Low-GQ star alleles

A star allele can rest on a single heterozygous variant called with a low genotype quality. Such a call prints like any other diplotype, so the summary looks each defining variant up in the VCF, and a call with one below GQ 20 gets `low-confidence: check the reads` in its Quality column and in both reports, with the weakest variant and its AD. Set the floor with `PYPGX_GQ_FLOOR=15 ./scripts/32-pypgx.sh your_name`, or for the Nextflow module with `withName: 'PYPGX' { ext.gq_floor = 15 }` in a config passed with `-c`. Mapping the defining variants back to the VCF uses pypgx's allele table; a variant pypgx names differently from the VCF (a multi-base change written another way) reads as `not in the VCF`, so check such a call by hand.

### In the reports

The pypgx section of both reports also lists the genes not assessable on this reference, the low-confidence calls and the depth check with its MAPQ shares. The Nextflow HTML report does not get the depth check file yet, so there the depth-check line appears only in the bash report. A **PGx copy number** block lists pypgx's copy-number call for each BAM-based gene and UGT2B17's copy number. UGT2B17 is a common whole-gene deletion with drug relevance that pypgx is not run on here: the block reads it from CNVpytor's calls (step 18) over the gene's span in the GENCODE gene lines `setup.sh` installs (`reference/gencode.v<N>.basic.genes.gtf`), as `no copy` (depth level below 0.25), `one copy` (below 0.75), `two copies` (no call over the gene) or `more than two copies` (a duplication), with the call's region, depth level, e-value and how much of the gene it covers. Without CNVpytor, or without the gene lines next to the results (a Nextflow report), it says `not assessed` and why. The CYP2D6 block shows pypgx's call labelled `probable (one caller: pypgx)` when step 36 allows it ([step 36](36-pgx-consensus.md)); it never reaches PharmCAT.

### PharmCAT comparison TSV columns

| Column | Description |
|---|---|
| Gene | Gene symbol |
| PharmCAT_diplotype | Diplotype from step 7 |
| pypgx_diplotype | Diplotype from this step |
| Match | Yes, No, pypgx only, or PharmCAT only. Allele names are compared after one rename: PharmCAT's `rs9923231 reference (C)` and pypgx's `Reference` are the same allele, and so are `rs9923231 variant (T)` and `rs9923231` |
| Called_by | `both` when both tools called the gene (Match is Yes or No), otherwise which tool did. PharmCAT's `Unknown/Unknown` and pypgx's `FAILED` and `N/A` count as no call; a gene neither tool called is left out |
| PharmCAT_phenotype | PharmCAT's phenotype for the gene, as step 27's phenotypes table shows it |
| pypgx_phenotype | The Phenotype column of the summary above |

## Runtime

5.5 to 10 minutes for 23 genes, peak memory about 4 GB, measured in three 30-40x WGS runs ([runtime](hardware-requirements.md#runtime-per-step)). All genes run sequentially inside a single Docker container to avoid repeated container startup overhead.

## Resource Requirements

```
--cpus 4 --memory 8g
```

The BAM-based SV detection (CYP2D6, CYP2A6, GSTM1, GSTT1) is the most memory-intensive phase due to read-depth calculation across the locus.

## Comparison with PharmCAT

| Aspect | PharmCAT (step 7) | pypgx (step 32) |
|---|---|---|
| Input | VCF only | BAM + VCF |
| CYP2D6 SVs | Cannot detect | Read-depth detection |
| GSTM1/GSTT1 | Not called | Deletion detection |
| COMT, MTHFR | Not covered | Covered |
| Drug recommendations | Yes (HTML/JSON report) | No (star alleles only) |
| CPIC integration | Built-in | Manual lookup via [CPIC guidelines](https://cpicpgx.org/guidelines/) |
| Validation | Widely used (research tool) | Less extensively validated |

The comparison reads PharmCAT 3.x reports in both shapes (a flat `genes` map or one nested by source), like step 27. If a PharmCAT report is found but cannot be parsed, or yields no gene, the step exits non-zero instead of writing a comparison where every gene is "pypgx only" and the report shows 0 conflicts.

The two tools are complementary. PharmCAT provides drug recommendations for the genes it covers; pypgx extends coverage to genes PharmCAT cannot call. When both tools call the same gene, concordance is expected for simple genotypes but discrepancies can occur for complex haplotypes. Neither tool is definitively "correct" in all cases — discrepancies should be investigated by examining the underlying variant calls. Neither tool constitutes a clinical test; results should be confirmed by a certified pharmacogenomics laboratory before making prescribing decisions.

## Note on Aldy

[Aldy](https://github.com/0xTCG/aldy) is widely considered the best CYP2D6 caller, with superior handling of complex structural rearrangements and hybrid alleles. However, Aldy is released under a custom academic-only license that is incompatible with GPL-3.0 redistribution. pypgx (Apache-2.0) is the GPL-compatible alternative used in this pipeline. If you are using this pipeline for personal/academic analysis and Aldy's license terms are acceptable, it can be run separately.

## Limitations

- pypgx gene coverage (88 total) is broader than the 23 curated here. The curated list focuses on CPIC Level A/B genes and key PharmCAT gaps. Edit the `BAM_GENES` and `VCF_GENES` variables in the script to add more.
- Star allele definitions evolve. pypgx 0.26.0 uses a specific PharmVar database snapshot that may not include the latest allele definitions.
- SV detection accuracy depends on sequencing depth. 30X WGS is adequate; lower depths produce less reliable copy number calls.
- pypgx does not produce drug recommendations directly. Its CYP2D6 call reaches PharmCAT, and so the CPIC recommendations of step 27, only through [step 36](36-pgx-consensus.md), when Cyrius (step 21, after `setup.sh --cyrius`) gives the same diplotype and the depth check passed. Alone, it is shown in the reports as a probable call and goes no further. Step 27 also lists the genes PharmCAT could not call but pypgx did; for the other genes consult [CPIC guidelines](https://cpicpgx.org/guidelines/).
- **The image and the bundle must be the same release.** The pinned pair is image `pypgx:0.26.0--pyh7e72e81_0` with the `0.26.0` branch of pypgx-bundle; with that pair all 23 genes, including the four BAM-based ones, were called on a real 30x genome (see [lessons learned](lessons-learned.md)). The 0.27.0 image against the 0.26.0 bundle failed every gene. Bump both together and rerun a known sample before trusting the new calls.
- **GSTT1 needs an ALT contig, which the default reference does not have.** In GRCh38, GSTT1 lies on `chr22_KI270879v1_alt`. The pipeline's reference, the no-ALT analysis set, has no such contig (it leaves ALT contigs out so CYP2D6 and the other BAM genes keep their depth, see [realignment](realignment.md)). The step then prints a notice, leaves GSTT1 out of depth preparation (otherwise it fails for all four SV genes) and writes GSTT1 with Status `not callable on this reference`, not `FAILED`: nothing went wrong, the reference cannot hold the gene. Both reports list it as not assessable on this reference. The other 22 genes are called.
- **A copy-number call is only as good as the depth it reads.** On a reference with ALT contigs and an aligner that is not run ALT-aware, reads at CYP2D6 split between the primary and ALT copies, depth on the primary drops, and pypgx can report a deletion that is not there. The [CYP2D6 depth check](#cyp2d6-depth-check) marks such a call `Indeterminate`, and step 36 reports CYP2D6 to PharmCAT only when two callers agree.
- Individual gene failures do not stop the pipeline. However, if **all** genes fail, the script exits with status 1 before generating the summary TSV — this signals a systemic problem (e.g., wrong BAM path, corrupted index, missing pypgx-bundle). Rerun with verbose output to identify the root cause. For partial failures, check the summary TSV for "FAILED" entries.

## Maintenance

- pypgx is pinned to `0.26.0` in `versions.env`. Check [pypgx releases](https://github.com/sbslee/pypgx/releases) periodically for updates to star allele definitions or algorithm improvements.
- The pypgx-bundle must match the pypgx version. When updating pypgx, re-download the bundle: `cd ${GENOME_DIR}/reference && rm -rf pypgx-bundle && git clone --branch <new_version> --depth 1 https://github.com/sbslee/pypgx-bundle.git`
- If you update pypgx, rerun on a known sample and compare diplotype calls against the previous version before adopting the new results.
- The curated gene list should be reviewed against [CPIC guideline updates](https://cpicpgx.org/guidelines/) at least quarterly.

## Links

- [pypgx documentation](https://pypgx.readthedocs.io/)
- [pypgx GitHub](https://github.com/sbslee/pypgx)
- [PharmVar database](https://www.pharmvar.org/) (star allele definitions)
- [CPIC guidelines](https://cpicpgx.org/guidelines/)
- [Aldy](https://github.com/0xTCG/aldy) (academic-only alternative for CYP2D6)
