# Step 40: Raw-Data Export (23andMe Layout)

## What This Does

Writes your genotypes at the sites of the 23andMe v5 chip into a file with the layout of a 23andMe raw-data download, and the same sites in the layout of an AncestryDNA download. Both are zipped. Opt-in.

## Why

Genealogy and health services that match relatives or read traits take the raw-data file a chip company gives its customers. A whole genome holds those sites too, but a VCF lists only the places where you differ from the reference, so a file made from it by hand misses most of them. This step reads the gVCF, which also records where you match the reference, so a site where you carry two reference bases gets them, and only a site the calls hold nothing for is `--`.

## Uploading is your decision

The step writes files on this machine and sends nothing anywhere. Uploading one to a service shares your genotypes at up to about 620,000 sites with it. A matching service also shows you to the relatives it finds, and them to you, and those relatives did not choose to be found. Read the service's terms first: who can see your matches, whether law enforcement can search the database, whether the file is used for research, and how to delete it. A genotype file identifies you and your close relatives whatever name you give it. Keep the file as private as the genome itself.

## Tool

- `bin/export_consumer.py` (this repository), standard-library Python. It walks the gVCF once.

## Docker Image

- `PYTHON_IMAGE`

Pinned in `versions.env`; [Image versions](versions.md) lists the current tag.

## Data: the site table

`setup.sh --consumer-sites` builds it once into `reference/consumer/23andme_v5_sites_dbsnp<build>.tsv.gz`:

1. The 23andMe v5 site list from [23andme2vcf](https://github.com/arrogantrobot/23andme2vcf) (MIT licence), at the commit `CONSUMER_SITES_COMMIT` and checked against `CONSUMER_SITES_SHA256`: 638,463 sites with their GRCh37 position and reference base.
2. dbSNP's GRCh38 VCF of build `CONSUMER_DBSNP_BUILD` (about 30 GB), checked against `CONSUMER_DBSNP_MD5`. Each rsID takes its GRCh38 position from it: no liftover chain is used. The download resumes when cut, and is deleted once the table is built.

A site is dropped from the table, with its reason in `reference/consumer/23andme_v5_dropped_dbsnp<build>.tsv`, when it has no rsID (23andMe's own `i` numbers, about 17,800 sites), dbSNP has no SNV for it, dbSNP places it at two positions or on another chromosome, or neither strand of its GRCh37 base is one of dbSNP's alleles. A site whose GRCh37 base is on the other strand in GRCh38 is written on the GRCh37 strand. A pseudoautosomal site, which dbSNP lists on chrX and chrY, keeps the chromosome of the list. The list numbers its MT sites the way hg19's chrM does, which is one or two bases off the rCRS; the table gives them dbSNP's rCRS position, the same in GRCh37 and GRCh38.

## Input

- `${SAMPLE}/vcf/${SAMPLE}.g.vcf.gz` from step 03. A sample with a VCF alone has no gVCF: run step 03 again from the BAM, or give the gVCF of the same calls.

## Command

```bash
./scripts/setup.sh --consumer-sites /path/to/genome_dir   # once: about 30 GB download, deleted afterwards
./scripts/40-export-consumer.sh your_name male            # or female; the sex is optional
```

With `run-all.sh`: `TOOLS=...,consumer_export` (it is opt-in, so a default run lists it as skipped; without the table it is skipped with the reason, and a run that starts from a VCF with no gVCF beside it is skipped too). In the Nextflow pipeline: `consumer_export` in `--tools` with `--consumer_sites`; the samplesheet's sex is the sex below, and a row without a gVCF is skipped with a warning.

## How each site is read

| The gVCF at the site | Written |
|---|---|
| a reference block (`GT 0/0`) | the reference base twice, for example `AA` |
| a variant record at the site, `PASS` | its two bases in alphabetical order, for example `AG` |
| a variant record that did not pass its filter, with a non-reference allele | `--` |
| a no-call (`./.`), no record, or a record that starts before the site and is not a block | `--` |
| inside the bases of a called deletion or a longer allele | `--` |

Then the chip's ploidy rules. With `male`, a chrX or chrY site outside the pseudoautosomal regions (`assets/par_grch38.bed`) gets one letter, and a heterozygous call there is `--`. With `female`, chrY sites are `--`. Without a sex each site keeps the ploidy of its call. A chrM site always gets one letter, and a heteroplasmic call (two bases) is `--`. A site whose strand is `-` in the table gets the complementary bases.

Before it writes, the step checks the reference base of every gVCF record that starts at a site against the table. When more than 1% of them (and more than 10) differ, the gVCF is not on GRCh38 or the table is wrong, and the step stops without writing a file.

## Output

| File | Contents |
|---|---|
| `consumer_export/${SAMPLE}_23andme_v5.zip` | `${SAMPLE}_23andme_v5.txt`: a comment header, then `rsid`, `chromosome` (1-22, X, Y, MT), `position` (GRCh37) and `genotype`, tab-separated |
| `consumer_export/${SAMPLE}_ancestrydna_layout.zip` | the same sites in the AncestryDNA layout: `rsid`, `chromosome` (1-22, 23 X, 24 Y, 25 the X pseudoautosomal regions, 26 MT), `position`, `allele1`, `allele2`; a no-call is `0 0`, a one-letter site has the letter twice |
| `consumer_export/${SAMPLE}_consumer_export.tsv` | what was written: the site count, how many got a genotype, how many came from reference blocks, and each reason for a `--` with its count |

The header says the file was written by this pipeline from sequencing calls, not by 23andMe or AncestryDNA. Before zipping, the step reads both files back with a strict reader (`export_consumer.py check`): the column header, the column count, an rsID, a chromosome and a position on every row, a genotype of one or two of `A`, `C`, `G`, `T` or `--`, one letter only on X, Y and MT, rows in chromosome and position order, each rsID once, plain ASCII with Unix line ends. A file that fails stops the step.

Both reports (step 24 and `generate-report.sh`) show the counts and the file names.

## Runtime

On the CI fixture (the gVCF of the fixture's 15 regions, 7,623 sites of the CI table) the step took under one second, container start included (e2e run 38649573454). It is not measured on a whole genome. The step reads the gVCF from start to end once in one Python process, so the time grows with the size of the gVCF; memory holds the site table, about 620,000 rows.

## Notes

- The AncestryDNA file has the 23andMe v5 sites, not the sites of an AncestryDNA chip; the two chips share part of their sites. A service that matches on the sites it expects uses the ones the files have in common.
- The 23andMe file has no internal `i` sites and none of the dropped sites, so it has fewer rows than a file from 23andMe.
- In a duplicated region the list's GRCh37 position and dbSNP's GRCh38 position can sit on two different copies. The file then carries the list's position with the genotype read at dbSNP's. On the CI table this is one site of 7,338 off chrM (in SMN1/SMN2), found by the round trip below.
- `scripts/chip-to-vcf.sh` reads both files back. The e2e case `consumer-export-1` does that on the CI fixture, lifts the result to GRCh38 and checks it against the fixture's calls; chrM is left out there, because the hg19 chain maps hg19's chrM, not the rCRS positions the file uses.
