# Step 33: Sample Identity and Contamination

## What This Does
Checks two things about a sample before you read any result from it:

- **Is it the person you think?** somalier reads the BAM at 17,766 known polymorphic sites and infers the sex from the reads: a male has no heterozygous sites on chrX outside the pseudoautosomal regions. A sex that differs from the one you declare stops the step. In a Nextflow run with several samples, somalier also compares every pair, so two rows that are the same person show up.
- **Is another person's DNA mixed in?** VerifyBamID2 estimates FREEMIX, the share of reads that come from someone else, from the allele balance at 100,000 markers of the 1000 Genomes panel. Above 0.03 the step warns.

## Why
A swapped sample gives a perfectly normal-looking report about someone else. Step 16 compares the declared sex with the X and Y copy numbers in the BAM index; this step checks again from the reads at known sites, and on a multi-sample run catches a duplicate or a swap between two people of the same sex, which a sex check cannot. Contamination is quieter: a few percent of foreign reads turn real homozygous sites into false heterozygous calls, so PharmCAT, ClinVar and the ROH steps all read a slightly wrong genotype. FREEMIX tells you how much to trust them.

Contamination only warns. A contaminated sample is still your sample; its calls are less reliable, above all the heterozygous ones, and you decide whether to resequence. A sex mismatch stops the step, because the steps that use the declared sex (DeepVariant's chrX and chrY ploidy, ExpansionHunter) would otherwise run with the wrong value.

## Tool
- **somalier** (Brent Pedersen): `extract` reads the sites, `relate --infer` infers the sex and the relatedness of every pair
- **VerifyBamID2** (Fan Zhang and Hyun Min Kang): FREEMIX with ancestry-aware allele frequencies

## Docker Image
- `SOMALIER_IMAGE`
- `VERIFYBAMID2_IMAGE`
- `PYTHON_IMAGE` (the verdict, from `bin/collect_summary.py`)

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
somalier relate --infer --sites "${SOMALIER_SITES}" -o qc/somalier/${SAMPLE} qc/somalier/<SM>.somalier
verifybamid2 --SVDPrefix "${VERIFYBAMID2_PANEL}" --Reference "${REF_FASTA}" \
  --BamFile aligned/${SAMPLE}_sorted.bam --Output qc/verifybamid2/${SAMPLE}
python3 bin/collect_summary.py sample-qc --sample ${SAMPLE} ...   # the verdict table
```

somalier names the sample after the BAM's `@RG SM` tag; every BAM step 02 writes has one.

VerifyBamID2 refuses to estimate when fewer than 1,000 panel markers have reads ("Insufficient Available markers"), as on a targeted (WES) or sliced BAM. The step then runs it again with `--DisableSanityCheck` and records `verifybamid2_marker_check skipped`, so the report says FREEMIX rests on fewer than 1,000 markers. On the CI fixture, with about 320 markers that have reads, FREEMIX still read 0.0004 for the clean HG002 and 0.093 for HG002 with 9.3% of its reads from HG001.

In Nextflow, add `sample_qc` to `--tools` with `--somalier_sites` and `--verifybamid2_panel` (the folder that holds the `.UD`, `.mu` and `.bed` files); `--freemix_warn` and `--sex_check` work as above. `SOMALIER_RELATE` runs once over every sample of the run. A sex mismatch stops the run, after `INDEXCOV` made the same check from the index.

## Output Files
| File | Description |
|---|---|
| `qc/${SAMPLE}_sample_qc.tsv` | The verdict, one `key` and `value` per line: `inferred_sex`, `sex_check` (`ok`, `mismatch` or `not_checked` with `sex_check_reason`), the chrX and chrY numbers it rests on, `freemix`, `contamination` (`ok` or `warn`), `panel_markers`, `verifybamid2_marker_check` (`passed`, or `skipped` when fewer than 1,000 markers had reads), `same_person_as` |
| `qc/somalier/${SAMPLE}.samples.tsv` | somalier's per-sample table: depth, genotype counts, chrX and chrY counts, inferred `sex` (1 male, 2 female, -9 unknown) |
| `qc/somalier/${SAMPLE}.pairs.tsv`, `.html` | Relatedness of every pair (one sample here; every sample of the run in Nextflow, under `somalier/`) |
| `qc/verifybamid2/${SAMPLE}.selfSM` | VerifyBamID2's result: `#SNPS` markers used, `AVG_DP`, `FREEMIX` |

The HTML report (step 24) shows the inferred sex, FREEMIX and any duplicate in its Quality Control card.

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

## Runtime
Both tools read the BAM only at their sites (17,766 for somalier, 100,000 for VerifyBamID2), not the whole file. On CI's fixture slice the step takes about a minute; a whole genome has not been timed here.

## Notes
- The CI fixture holds only small slices of HG002, with 2 of somalier's chrX sites and 320 of VerifyBamID2's markers: somalier reports `unknown` there, and VerifyBamID2 runs without its marker check. The end-to-end case adds sites at the slice's own chrX calls to show the sex check stop a female-declared HG002, and mixes about 10% of HG001's reads into HG002 to show FREEMIX rise above 0.03.
- somalier stops with `sequence chr11 not found in fasta` when the reference lacks a contig its sites file names. Every full GRCh38 FASTA has chr1 to chr22, chrX and chrY; for a cut-down reference, give `SOMALIER_SITES` a copy of the sites without the missing contigs, as the end-to-end case does for its 14-contig fixture reference.
- VerifyBamID2's panel is GRCh38 (`b38`); for a GRCh37 BAM the pipeline would need the `b37` files, which it does not install.
- somalier compares samples only inside one `somalier relate` call: a Nextflow run relates every sample of its samplesheet. To compare samples run at different times with the bash step, run `somalier relate --infer --sites <sites> -o out */qc/somalier/*.somalier` in `SOMALIER_IMAGE` over their `.somalier` files.
