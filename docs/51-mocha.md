# Step 51: Mosaic Chromosomal Alterations (MoChA)

## What This Does

Looks for large mosaic chromosomal alterations (mCAs) on chr1 to chr22: a gain, a loss or a copy-neutral loss of heterozygosity (CN-LOH) of millions of bases that a share of the blood cells carries, not every cell. It reports each one with its type, where it is, the chromosome arm it covers and the share of cells that carry it (the cell fraction). Opt-in.

Mosaic loss of Y and of X are [step 33](33-sample-qc.md)'s. This step reads the autosomes only.

## Why

mCAs in blood become more common with age. [Loh et al. 2018](https://pubmed.ncbi.nlm.nih.gov/29995854/) (PMID 29995854) called 8,342 mCAs in 151,202 UK Biobank participants, and people who carried one had a higher risk of a later blood cancer. Most carriers never get one. The same paper introduced the method this step runs: in the cells that carry the event one parental copy of the stretch is lost or gained, so at every heterozygous site of the stretch the allele balance leans the same way along one phased haplotype. Following the phase makes a small lean over many sites visible.

## Tool

- [MoChA](https://github.com/freeseek/mocha) (`bcftools +mocha` and `bcftools +mochatools`, MIT licence), the authors' implementation of that method, version `MOCHA_VERSION` (2025-08-19, commit `MOCHA_COMMIT`)
- Beagle to phase the VCF against the 1000 Genomes high-coverage panel, as the MoChA README does
- `bin/mocha_report.py` (standard library Python) for the floor and the calls reported

## Docker Image

- `MOCHA_IMAGE` (bcftools 1.22 with the MoChA plugins, the authors' image, pinned by digest), `BEAGLE_IMAGE` and `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tags. MoChA has no tagged release, so `MOCHA_COMMIT` and `MOCHA_VERSION` move with the image, by hand. The e2e case checks that the image prints `MOCHA_VERSION` and that `mocha.c` at `MOCHA_COMMIT` defines it.

## Reference data

`./scripts/setup.sh --mocha-data <genome_dir>` installs `reference/mocha-<MOCHA_PANEL_VERSION>/`:

- `panel/chrN.vcf.gz`: the 1000 Genomes high-coverage phased panel of 3,202 people (release 20220422, EBI), one chromosome at a time. Each file is checked against the md5 of the release's manifest, and the manifest against `MOCHA_PANEL_MANIFEST_SHA256`. Only the biallelic SNVs seen twice or more are kept, with genotypes only; the downloaded file is then deleted.
- `maps/chrN.map`: the PLINK GRCh38 genetic map of step 26b (`ANCESTRY_FINE_MAP_URL`, its mirror and its sha256).
- `cnps.bed`: the panel's common germline deletions and duplications, MoChA's `--cnp` list. Each chromosome file gives its structural variants seen twice or more and over 10 kb (`SVTYPE` CNV, DEL or DUP) before it is deleted, the recipe of the MoChA README (GRCh38, "List of common germline duplications and deletions"). A panel installed before this list existed is rebuilt by a rerun of `--mocha-data`.

The 22 downloads come to about 30 GB (29.8 GB by the server's sizes on 2026-10-27). They are fetched one at a time and each is deleted once cut, so what stays on disk is much smaller. `MOCHA_REGIONS=<bed>` keeps only the regions of a BED file (the e2e case cuts the panel to one 500 kb slice of chr20 this way). `MOCHA_DATA` overrides the folder.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.vcf.gz` with its `.tbi`, with GT, AD, DP and GQ (DeepVariant writes all four)
- the panel above
- `MOCHA_PHASED_VCF=<file inside GENOME_DIR>`: a VCF already phased (GT with `|`), read instead; no Beagle run and no panel needed. The bash step only (`scripts/51-mocha.sh`); the Nextflow processes always phase with Beagle. Without the panel folder's `cnps.bed`, MoChA then runs with no `--cnp` list and says so

## Command

```bash
./scripts/setup.sh --mocha-data /path/to/genome_dir
./scripts/51-mocha.sh your_name
```

With `run-all.sh`: `TOOLS=...,mocha`; a run lists step 51 only when `TOOLS` names it. In the Nextflow pipeline: `mocha` in `--tools` with `--mocha_data reference/mocha-20220422` (processes `MOCHA_PREP`, `MOCHA_PHASE`, `MOCHA_CALL` and `MOCHA_REPORT`).

## How it works

1. The VCF's biallelic PASS SNVs on chr1 to chr22. A genotype under 10 reads or with GQ under 20 is set missing, as the MoChA README does for WGS. `mochatools` adds the GC content around each site from the reference.
2. Beagle phases the genotypes against the panel, one chromosome at a time (`impute=false`, `seed=1`). A chromosome the panel lacks is left out.
3. The phased genotypes go back onto the sites of step 1, and MoChA runs in its WGS mode (allele depths from AD) with `--genome GRCh38` and `cnps.bed` as `--cnp`. The MHC and KIR regions are left out of its model.
4. `bin/mocha_report.py` reads MoChA's tables and decides which calls are reported.

A call is reported when all of these hold:

- it is on chr1 to chr22;
- MoChA's own filter of the UK Biobank callset (MoChA README, "Filter callset") keeps it: not a germline copy-number polymorphism, not a likely constitutional event (a BAF deviation of 0.1 or more with sparse heterozygous sites), not a likely germline duplication;
- MoChA estimated its cell fraction, and it is at or above the floor.

**The floor** is the larger of two numbers:

- 0.05. This one is our choice: below a cell fraction of 0.05 the pipeline reports nothing, whatever the depth.
- the read-sampling floor at the sample's depth, `6 / sqrt(D x H)`, with `D` MoChA's median depth and `H` the heterozygous sites per 10 Mb. This is the cell fraction at which a 10 Mb loss or gain of one copy moves the allele balance by three standard errors of read sampling alone. A loss or gain of one copy in a fraction `c` of the cells moves the balance by about `c / 4`; CN-LOH moves it by `c / 2`, so it shows at half that floor.

At 30x with about 3,000 heterozygous sites per 10 Mb the read-sampling floor is about 0.02, so the floor is 0.05. At 5x with 300 sites per 10 Mb it is about 0.15. The floor is a lower bound: phasing errors and read bias are not in it, so an event above it is not always found.

The sample is **not assessed** when MoChA's sample filter of the same README fails (a call rate under 0.97, or a phased BAF autocorrelation over 0.03), when MoChA gives no value for one of the two, or when no heterozygous site is phased. Then no call is reported, and the reports say why.

## Output

| File | Contents |
|---|---|
| `mosaic/${SAMPLE}_mocha.tsv` | one row per MoChA call: `chrom`, `start`, `end`, `length_mb`, `arms` (for example "20p, to the p telomere"), `type` (loss, gain, CN-LOH, undetermined), `cell_fraction`, `bdev`, `rel_cov`, `n_hets`, `lod_baf_phase`, `lod_lrr_baf`, `reported` (yes or no) and `reason` |
| `mosaic/${SAMPLE}_mocha_summary.tsv` | `metric`, `value`: status (assessed or not assessed) and why, depth, heterozygous sites (all and phased), sites per 10 Mb, the floor and its two parts, the calls MoChA made and the calls reported, and how the VCF was phased |
| `mosaic/${SAMPLE}_mocha_calls_raw.tsv`, `_mocha_stats.tsv` | MoChA's own `--calls` and `--stats` tables |
| `mosaic/${SAMPLE}_mocha_sites.tsv` | per chromosome: first and last position, sites, heterozygous sites, phased heterozygous sites |
| `mosaic/${SAMPLE}_mocha.log` | MoChA's log |

Both reports show a Mosaic Chromosomal Alterations (MoChA) section with the reported calls and the floor. Without step 51 the report lists autosomal mosaic alterations as not assessed.

## Runtime

Not measured on a whole genome. Beagle's phasing of the 22 chromosomes takes most of the time; it runs with `THREADS` threads and a 14 GB heap (16 GB for the container). MoChA itself reads one phased BCF.

## Notes

- Research result, not a diagnosis. An mCA in blood is common with age and most carriers never get a blood cancer. A clinical test for clonal blood cells is done in a clinical laboratory.
- Blood only. The method assumes the DNA comes from blood cells; saliva holds a mix of cells, which changes what a cell fraction means.
- A constitutional event (one every cell carries, such as a germline deletion) has a cell fraction near 1 and is left out by MoChA's filter; this step is not a test for those.
- A germline copy-number polymorphism is left out by type: MoChA types a call over a region of `cnps.bed` as `CNP_Loss` or `CNP_Gain`. A germline CNV missing from that list (rare, or under 10 kb) gets no such type. Short events (under a few Mb) are hard to tell from these and are reported only when MoChA's filter keeps them.
- Calls are labelled by chromosome arm, not by gene: an mCA spans millions of bases.
- The e2e case (`tests/e2e/mocha-1.sh`) runs the step on HG002's chr20 slice (MoChA makes no call; with about 500 heterozygous sites the sample filter can mark the slice not assessed, since the BAF autocorrelation it reads is meant for a whole genome), and on a synthetic phased VCF with a 20 Mb CN-LOH on 20p in 15% of the cells (one CN-LOH call, cell fraction 0.10 to 0.20) beside the same VCF with no shift (no call). Nothing in the case is a finding about the fixture's sample.
