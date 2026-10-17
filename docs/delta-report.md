# What changed since the previous run

## What This Does

When you run the pipeline again on the same sample, for example after a ClinVar update or a new tool version, the new reports replace the old ones. The delta report says what moved between the two runs, so you do not have to compare two long reports by hand.

`scripts/generate-report.sh` writes it after the text report. `run-all.sh` runs that script last, so every `run-all.sh` run gets one.

## How it works

1. Each time `generate-report.sh` runs, it keeps a copy of the run's `summary.json` and `run_manifest.tsv` in `${SAMPLE}/history/<UTC time>/`, for example `history/2026-10-18T091500Z/`.
2. Before it keeps that copy, `bin/run_delta.py` compares the current `summary.json` and `run_manifest.tsv` with the newest folder in `history/`.
3. It writes `${SAMPLE}/${SAMPLE}_delta.txt`, `${SAMPLE}_delta.html` and `${SAMPLE}_delta.json`, and prints the text version at the end of the script's output.

The first run with no `history/` folder says there is no previous run and becomes the baseline. A run whose inputs and results match the previous one says "No change".

Runs made before this feature kept no history, so the first run after it has nothing to compare with.

## What it compares

| Area | Compared |
|---|---|
| Inputs | Every row of the run manifest: each `versions.env` value, the image each `*_IMAGE` names with the digest Docker resolved it to, the ClinVar file date, the VEP cache, the PCGR and pypgx bundles, the ancestry panel, the HLA database release, the header of each PGS scoring file, the pipeline commit and the declared sex. The time the manifest was written and the script that wrote it are left out |
| Report sections | A section that ran in one run and not in the other (for example `ok` to `missing` or `failed`) |
| ClinVar (step 06) | Hits added or dropped, a hit whose significance changed (reclassified), whose review status or stars changed (re-starred), or whose genotype changed. A hit is matched by chromosome, position, REF and ALT |
| CPSR (step 17) | The number of variants in each classification |
| Pharmacogenomics | Each CPIC gene's diplotype and phenotype (step 27), pypgx's CYP2D6 call (step 32) and the CYP2D6 consensus (step 36) |
| Polygenic risk scores (step 25) | Scores added or dropped, and a percentile that moved by 5 points or more, or appeared or went away. A smaller move is left out |

The other sections (structural variants, repeat expansions, HLA alleles, haplogroups, QC metrics) are compared only by whether they ran.

For each change in a result, the delta lists the changed inputs that feed that section: for ClinVar, the ClinVar file date and the fastp, aligner, DeepVariant, GATK and bcftools images; for CPSR, the PCGR image and data, the VEP image and cache, and DeepVariant; for pharmacogenomics, the PharmCAT, pypgx and Cyrius versions and the pypgx bundle (for the CPIC genes, DeepVariant and bcftools too); for PRS, the PGS headers, the pgsc_calc images and panel, plink2 and DeepVariant. A new pipeline commit counts for every section. When none of those changed, it says "the inputs this delta tracks for it are the same; the call itself differs". That line is a pointer, not a proof of cause.

## Command

```bash
export GENOME_DIR=/path/to/your/data
./scripts/generate-report.sh your_sample
cat ${GENOME_DIR}/your_sample/your_sample_delta.txt
```

`run-all.sh` rewrites `run_manifest.tsv` before the reports. After you rerun a single step by hand, write the manifest again before the report, so the delta sees the new inputs:

```bash
GENOME_DIR=/path/to/your/data bin/write_manifest.sh your_sample manual
./scripts/generate-report.sh your_sample
```

Each run of `generate-report.sh` adds a folder to `history/`, so running it twice in a row gives an empty delta the second time. The folders are small (a summary and a manifest); delete old ones whenever you like. The newest one is always the baseline.

## Nextflow

A bare `nextflow run main.nf` writes no text report, so it writes no delta and keeps no history either. `run-all.sh` writes both after the pipeline. On the output of a bare Nextflow run, run `24-html-report.sh` and then `generate-report.sh` with `GENOME_DIR` set to the `--outdir`.

## When to rerun

| Resource | How often it changes | What to rerun |
|---|---|---|
| ClinVar | Monthly | `./scripts/setup.sh --refresh clinvar <genome_dir>`, then the same `run-all.sh` command |
| PharmCAT | With each release | Steps 07 and 27 |
| PGS Catalog | When a score you use gets a new version | Step 25 |

`run-all.sh` runs Nextflow with `-resume`, so running the same command again redoes only the processes whose inputs or code changed. A new ClinVar file reruns the ClinVar screen and what reads it; the alignment and the variant calls come from the cache. The delta then shows what the new release moved.
