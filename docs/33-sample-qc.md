# Step 33: Sample Identity and Contamination

## What This Does
Checks two things about a sample before you read any result from it, and estimates a third:

- **Is it the person you think?** somalier reads the BAM at 17,766 known polymorphic sites and infers the sex from the reads: a male has no heterozygous sites on chrX outside the pseudoautosomal regions. A sex that differs from the one you declare stops the step. In a Nextflow run with several samples, somalier also compares every pair, so two rows that are the same person show up.
- **Is another person's DNA mixed in?** VerifyBamID2 estimates FREEMIX, the share of reads that come from someone else, from the allele balance at 100,000 markers of the 1000 Genomes panel. Above 0.03 the step warns.
- **Did some cells lose the Y, or one X?** From the depth step 16b measured and, for a woman, the VCF's chrX allele fractions, the step estimates the share of cells that lost the Y chromosome (men) or one X (women), with its 95% range. See [Mosaic loss of chrY or chrX](#mosaic-loss-of-chry-or-chrx).

## Why
A swapped sample gives a perfectly normal-looking report about someone else. Step 16 compares the declared sex with the X and Y copy numbers in the BAM index; this step checks again from the reads at known sites, and on a multi-sample run catches a duplicate or a swap between two people of the same sex, which a sex check cannot. Contamination is quieter: a few percent of foreign reads turn real homozygous sites into false heterozygous calls, so PharmCAT, ClinVar and the ROH steps all read a slightly wrong genotype. FREEMIX tells you how much to trust them.

Contamination only warns. A contaminated sample is still your sample; its calls are less reliable, above all the heterozygous ones, and you decide whether to resequence. A sex mismatch stops the step, because the steps that use the declared sex (DeepVariant's chrX and chrY ploidy, ExpansionHunter, Delly) would otherwise run with the wrong value.

## Tool
- **somalier** (Brent Pedersen): `extract` reads the sites, `relate --infer` infers the sex and the relatedness of every pair. It gets the declared sex as a one-line pedigree (`--ped`), so its own table, and MultiQC's somalier section, show the sex you declared beside the one from the reads
- **VerifyBamID2** (Fan Zhang and Hyun Min Kang): FREEMIX with ancestry-aware allele frequencies

## Docker Image
- `SOMALIER_IMAGE`
- `VERIFYBAMID2_IMAGE`
- `PYTHON_IMAGE` (the verdict, from `bin/collect_summary.py`, and the mosaic estimate, from `bin/mosaic_sex_loss.py`)

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags.

## Data
`setup.sh` installs both, each checked against a pinned sha256 (`./scripts/setup.sh --sample-qc-data ${GENOME_DIR}` installs only these):

| File | What it is |
|---|---|
| `reference/somalier/sites.hg38.vcf.gz` | somalier's GRCh38 sites, 17,766 common SNPs on chr1-22, chrX and chrY (265 kB) |
| `reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat.{UD,mu,bed}` | VerifyBamID2's panel of 100,000 1000 Genomes markers, from the VerifyBamID v2.0.3 release (10 MB) |

## Command
```bash
./scripts/33-sample-qc.sh your_name male      # or female: the sex check runs
./scripts/33-sample-qc.sh your_name           # no declared sex: report only
```

| Setting | Default | What it does |
|---|---|---|
| `SEX_CHECK` | `fail` | `warn` prints a mismatch and exits 0, as in step 16 |
| `FREEMIX_WARN` | `0.03` | FREEMIX above this is reported as possible contamination |
| `SOMALIER_SITES` | `reference/somalier/sites.hg38.vcf.gz` | Another sites VCF, inside `GENOME_DIR` |
| `VERIFYBAMID2_PANEL` | `reference/verifybamid2/1000g.phase3.100k.b38.vcf.gz.dat` | Another panel, as its prefix (the path without `.UD`) |

The script runs, in its containers:
```bash
somalier extract -d qc/somalier --sites "${SOMALIER_SITES}" -f "${REF_FASTA}" aligned/${SAMPLE}_sorted.bam
somalier relate --infer --ped qc/somalier/${SAMPLE}.declared.ped --sites "${SOMALIER_SITES}" -o qc/somalier/${SAMPLE} qc/somalier/<SM>.somalier
verifybamid2 --SVDPrefix "${VERIFYBAMID2_PANEL}" --Reference "${REF_FASTA}" \
  --BamFile aligned/${SAMPLE}_sorted.bam --Output qc/verifybamid2/${SAMPLE}
python3 bin/collect_summary.py sample-qc --sample ${SAMPLE} ...   # the verdict table
python3 bin/mosaic_sex_loss.py --table qc/${SAMPLE}_sample_qc.tsv \
  --regions mosdepth/${SAMPLE}.regions.bed.gz --vcf vcf/${SAMPLE}.vcf.gz   # the mosaic rows, each input when it is there
```

The pedigree is one line: the sample, itself as its family, no parents (`-9`), the declared sex (`1` male, `2` female, `-9` none) and no phenotype. With it somalier starts its `sex` column from the declared sex and changes it only when the reads tell; it then logs `setting sex to male for <sample>` (or female). The step's verdict does not read the declared sex back as a call: when chrX cannot tell (10 or fewer chrX sites, a chrX het/hom-ALT ratio between 0.05 and 0.4, or 6% or more sites with an allele balance outside 0.1 to 0.9, somalier's own rule), the sex from the reads is `unknown` and the check is `not_checked`, as without a pedigree.

The bash step lets somalier name the sample after the BAM's `@RG SM` tag; every BAM step 02 writes has one. In Nextflow, somalier names each sample by its samplesheet id instead, so two rows whose BAMs share an `SM` stay two samples, and a BAM with no read group works.

VerifyBamID2 refuses to estimate when fewer than 1,000 panel markers have reads ("Insufficient Available markers"), as on a targeted (WES) or sliced BAM. The step then runs it again with `--DisableSanityCheck` and records `verifybamid2_marker_check skipped`, so the report says FREEMIX rests on fewer than 1,000 markers. On the CI fixture, with about 320 markers that have reads, FREEMIX still read 0.0004 for the clean HG002 and 0.093 for HG002 with 9.3% of its reads from HG001.

In Nextflow, add `sample_qc` to `--tools` with `--somalier_sites` and `--verifybamid2_panel` (the folder that holds the `.UD`, `.mu` and `.bed` files); `--freemix_warn` and `--sex_check` work as above. `SOMALIER_RELATE` runs once over every sample of the run, with a pedigree of every sample's samplesheet sex. A sex mismatch stops the run, after `INDEXCOV` made the same check from the index. `SAMPLE_QC` gets `MOSDEPTH`'s windows when `mosdepth` is in `--tools` (it is in the default list), and the row's VCF. It waits for that VCF, so on a FASTQ row somalier's sex check reports after `DEEPVARIANT`; `INDEXCOV`'s check still stops the run before it.

## Mosaic loss of chrY or chrX
With age, some blood cells lose their Y chromosome (men) or one of their two X chromosomes (women), and the cells that lost it can grow into a clone. DNA from blood then holds less of that chromosome than one copy (Y) or two copies (X) per cell. Both become more common with age. Population studies link loss of Y with smoking and with higher rates of heart disease and some cancers; loss of X is less studied. That is an association across many people, not a diagnosis for one, so the step reports an estimate with its uncertainty and nothing more.

The sex is the one somalier infers from the reads, else the declared one; with neither, nothing is estimated. The step estimates the share of cells with the loss two ways:

- **Depth.** From step 16b's mosdepth windows (`mosdepth/${SAMPLE}.regions.bed.gz`; `coverage/` in a Nextflow run). The 500 bp windows are grouped in 10 kb blocks, each block the median of its windows with reads. The ratio is the median block of the sex chromosome over the median block of chr1 to chr22, counting only the single-copy part of the sex chromosome: chrY or chrX outside the pseudoautosomal regions, the X-transposed region, the ampliconic regions and the centromere, from GIAB's GRCh38 genome stratifications v3.3, in pieces of 200 kb or more (9.7 Mb of chrY, 143 Mb of chrX). That is the normalised ratio. A whole-chromosome mean, which mosdepth's summary gives, also counts chrY's gaps and repeats and the copies of chrX that chrY holds. A man with one Y in every cell reads 0.5 for chrY, a woman 1.0 for chrX. With a share f of the cells missing the copy, f = 1 - 2 × ratio for chrY and f = 2 × (1 - ratio) for chrX.
- **Allele fractions, women only.** From the VCF (`vcf/${SAMPLE}.vcf.gz`, its `GT` and `AD`), at heterozygous PASS SNVs on chrX's single-copy part with 10 or more reads. When one X is lost in a share f of the cells, the allele fraction at those sites moves from 0.5 to 0.5 ± d, with d = f / (2 (2 - f)). The spread a library has anyway is read from the heterozygous SNVs of chr1 to chr22, so a noisy library does not read as a shift. Allele fractions also show a loss in which the other X was copied (copy-neutral), which depth cannot see.

Each estimate has a 95% range, from the spread of the blocks or of the sites. Below about 10% of the cells nothing can be told at 30x: the depth of a sex chromosome against the autosomes varies by a few percent from one library to the next (GC content, mappability), and that is 5 to 10% of the estimate. So the call is `detected` only when the low end of a range is above 10%, and `not_detected` otherwise. An extra copy (47,XYY, 47,XXX) reads as `gain`. A constitutional mosaic such as 45,X/46,XX reads in blood the same as a loss with age; the step cannot tell them apart.

The step needs step 16b's windows for the depth estimate and, for a woman, the VCF for the allele fractions. Without them it writes `not_estimated`, or the one estimate it could make, and says what was missing in `mosaic_reason`. A missing or unreadable input never stops the step. A targeted (WES) or sliced BAM has too few 10 kb blocks on the sex chromosome (each estimate needs 10) and gets `not_estimated`.

## Output Files
| File | Description |
|---|---|
| `qc/${SAMPLE}_sample_qc.tsv` | The verdict, one `key` and `value` per line: `inferred_sex`, `sex_check` (`ok`, `mismatch` or `not_checked` with `sex_check_reason`), the chrX and chrY numbers it rests on, `freemix`, `contamination` (`ok` or `warn`), `panel_markers`, `verifybamid2_marker_check` (`passed`, or `skipped` when fewer than 1,000 markers had reads), `same_person_as`; then the mosaic rows below |
| `qc/somalier/${SAMPLE}.samples.tsv` | somalier's per-sample table: depth, genotype counts, chrX and chrY counts, `original_pedigree_sex` (the declared sex: `male`, `female`, or `unknown` when none was declared), `sex` (the pedigree's code, changed when the reads tell: 1 male, 2 female, -9 unknown, -2 when it is female but chrY has reads) |
| `qc/somalier/${SAMPLE}.declared.ped` | The pedigree somalier got: the declared sex |
| `qc/somalier/${SAMPLE}.pairs.tsv`, `.html` | Relatedness of every pair (one sample here; every sample of the run in Nextflow, under `somalier/`) |
| `qc/verifybamid2/${SAMPLE}.selfSM` | VerifyBamID2's result: `#SNPS` the size of the panel, `AVG_DP`, `FREEMIX` |

The mosaic rows, every one always present (empty when it does not apply):

| Key | Value |
|---|---|
| `mosaic_sex`, `mosaic_sex_from` | The sex the estimate assumes, and where it came from: `reads` (somalier) or `declared` |
| `mosaic_chrom` | `chrY` for a man, `chrX` for a woman |
| `mosaic_depth_ratio`, `mosaic_depth_ratio_expected` | The normalised depth ratio, and what one Y (0.5) or two X (1.0) in every cell read |
| `mosaic_depth_blocks` | The 10 kb blocks of the sex chromosome it rests on |
| `mosaic_loss_depth`, `_low`, `_high` | The share of cells with the loss, from depth, and its 95% range (a fraction: 0.25 is 25%) |
| `mosaic_x_het_sites` | Women: the heterozygous chrX SNVs read |
| `mosaic_loss_af`, `_low`, `_high` | Women: the share of cells with the loss, from allele fractions, and its 95% range |
| `mosaic_floor` | 0.100: below this, nothing can be told at 30x |
| `mosaic_call` | `detected`, `not_detected`, `gain` or `not_estimated` |
| `mosaic_reason` | Why, in words, and any input that was missing |

The HTML report (step 24) shows the inferred sex, FREEMIX and any duplicate in its Quality Control card, and the mosaic estimate in its own card, with the same sentence on what it means as the text report.

## Interpretation
| Result | Meaning | What to do |
|---|---|---|
| Sex check `ok` | somalier's sex matches yours | Nothing |
| Sex check `mismatch` | The reads say the other sex | Check the sample's origin first. A sex-chromosome aneuploidy (47,XXY, 45,X) is the other explanation; then `SEX_CHECK=warn` |
| Sex `unknown` | Fewer than 11 chrX sites have a genotype, or allele balance looks like more than one person | Nothing on a targeted or sliced BAM; on a whole genome, look at FREEMIX |
| FREEMIX below 0.01 | Clean | Nothing |
| FREEMIX 0.01 to 0.03 | Trace contamination, usual for a library | Nothing |
| FREEMIX above 0.03 | About that share of the reads are someone else's | Treat heterozygous calls with care; ask the lab whether it can resequence |
| `same_person_as` set | Two samples of the run share their genome (relatedness 0.9 or more) | A duplicate row, a resequenced sample, an identical twin, or a swap |
| `mosaic_call not_detected` | No loss of the Y (or of one X) above about 10% of the cells | Nothing; a smaller loss cannot be told from noise at 30x |
| `mosaic_call detected` | The low end of the 95% range is above 10%: a share of the blood cells lost the Y, or one X | Nothing follows from it alone. It is a population-level marker, not a diagnosis; a doctor would test blood cells directly if it ever mattered |
| `mosaic_call gain` | More depth than one Y or two X | An extra copy (47,XYY, 47,XXX) or a gain; it is not a loss |

## Runtime
Both tools read the BAM only at their sites (17,766 for somalier, 100,000 for VerifyBamID2), not the whole file. The mosaic estimate reads step 16b's windows file and, for a woman, the VCF once each; it has not been timed on a whole genome. On CI's fixture slices one run of the step takes a few seconds (six runs of it and the HG001 stream took 58 seconds on one e2e run); a whole genome has not been timed here.

## Notes
- The CI fixture holds only small slices of HG002, with 2 of somalier's chrX sites and 320 of VerifyBamID2's markers: somalier reports `unknown` there, and VerifyBamID2 runs without its marker check. The end-to-end case adds sites at the slice's own chrX calls to show the sex check stop a female-declared HG002, and mixes about 10% of HG001's reads into HG002 to show FREEMIX rise above 0.03.
- somalier stops with `sequence chr11 not found in fasta` when the reference lacks a contig its sites file names. Every full GRCh38 FASTA has chr1 to chr22, chrX and chrY; for a cut-down reference, give `SOMALIER_SITES` a copy of the sites without the missing contigs, as the end-to-end case does for its 14-contig fixture reference.
- VerifyBamID2's panel is GRCh38 (`b38`); for a GRCh37 BAM the pipeline would need the `b37` files, which it does not install.
- somalier compares samples only inside one `somalier relate` call: a Nextflow run relates every sample of its samplesheet. To compare samples run at different times with the bash step, run `somalier relate --infer --sites <sites> -o out */qc/somalier/*.somalier` in `SOMALIER_IMAGE` over their `.somalier` files.
