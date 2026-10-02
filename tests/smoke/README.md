# Image tests

Every container image in `versions.env` runs a real command on the e2e fixture
here, and the test checks what the tool wrote. The **Container Test** workflow
(`.github/workflows/container-test.yml`) runs the rows of the images a pull
request changes, and every row once a month.

- `commands.tsv`: the rows, one or more per `*_IMAGE` variable.
- `*.sh`, `*.vcf.in`, `*.toml`: files a row needs, mounted at `/smoke` (`*.vcf` is git-ignored, so a VCF here is named `.vcf.in`).
- `scripts/ci/image-smoke.sh`: prepares the inputs, runs the rows, checks them.
- `scripts/ci/changed-images.sh`: picks the images a change must run.

An image with no row fails the workflow, so a new image cannot be added
untested.

## Run it

```bash
scripts/ci/image-smoke.sh --check            # table only, no docker
scripts/ci/image-smoke.sh PYPGX_IMAGE        # the rows of one image (docker, gh)
scripts/ci/image-smoke.sh --full ANNOTSV_IMAGE
```

To run images on a branch in CI, start the Container Test workflow by hand
with `images` set to `all`, `changed` or a list of variable names.

## Columns

Five columns, tab separated. A line starting with `#` is a comment.

| Column | Meaning |
|---|---|
| var | The `*_IMAGE` variable of `versions.env`. |
| opts | `-`, or a comma list: `root` (run as root), `net` (allow network), `full` (monthly and dispatch only), `also=VAR` (a coupled data variable: a change to it runs this row too), `rw=INPUT` (that input is writable). |
| needs | `-`, or a comma list of inputs (below). |
| command | Runs as `sh -c` in the image, in an empty `/out`, inputs at `/in` (read-only), this folder at `/smoke`, every `versions.env` variable in the environment, no network unless `net`. |
| expect | Bash, run on the host in the row's output folder after the command exits 0. Every check prints `[PASS]` or `[FAIL]`. |

Checks: `at_least NAME VALUE MIN`, `between NAME VALUE LO HI`,
`same NAME VALUE WANT`, `has NAME VALUE ERE`, `contains NAME FILE ERE`,
`nonempty NAME FILE`, `py NAME CODE`. Values: `vcf_records`, `mapped_pct`,
`snv_recall` (share of the GIAB truth SNVs a VCF calls), `sv_near` (records
at the planted deletion), `vcf_field`, `tsv_cell`, `tsv_count`, `csv_cell`,
`pharmcat_called`. An exit code alone is never the check.

## Inputs

| Name | What it is |
|---|---|
| ref | `mini.fa`: chr20:10.0-10.5 Mb of the fixture reference, with 3,000 random bases inserted after position 450,000, so the sample shows a 3 kb deletion there. `orig.fa` is the slice without them. |
| reads | Read pairs of that slice, `reads_R1.fq.gz` and `reads_R2.fq.gz`. |
| bam | Those reads aligned to `mini.fa` with minimap2: `mini.bam`. |
| truth | The GIAB v4.2.1 calls of the slice up to position 449,000 (`truth.vcf.gz`, `truth.bed`), and `query90.vcf.gz`, the truth without every tenth SNV. |
| longreads, longbam | Synthetic HiFi reads from both truth haplotypes, and their alignment. |
| fullref | The whole fixture reference, `ref.fa`. |
| slice | The fixture's GIAB BAM, `HG002_slice.bam`. |
| vcf, vcf50 | The fixture's VEP-annotated VCF, and its first 50 records without CSQ. |
| revel, sv, cyrius | Fixture files: synthetic REVEL scores, the ten-record SV VCF, the Cyrius BAM. |
| dels | Two deletions for duphold: the planted one and a control. |
| bundle | The pypgx-bundle tag `PYPGX_BUNDLE_VERSION` names. |
| mito | chrM calls from the slice BAM. |
| hlareads | Read pairs of the HLA slice (chr6:29.9-33.1 Mb). |
| qc | samtools flagstat and stats of `mini.bam`. |
| ehcatalog | A one-locus ExpansionHunter catalog for a CA repeat in `mini.fa`. |
| chain | A chain file that maps `mini.fa` onto itself. |
| annotsv | AnnotSV's annotation data (5.3 GB download). |

## What the runner cannot test

- VEP runs with `--database` on 50 variants. The offline cache (about 26 GB)
  does not fit a runner.
- PCGR: the CPSR report needs the PCGR data bundle and a VEP 113 cache. The row
  runs what the image carries: `cpsr --version`, its VEP on 50 variants with
  `--database`, and the R package that writes the report.
- AnnotSV annotates only in the monthly and dispatched runs (5.3 GB of data).
  A pull request that changes it gets a warning to start the workflow by hand.
