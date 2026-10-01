# Testing

This page is for contributors. It says what CI runs, how the end-to-end test works, and how to change its data.

The stub run in `nextflow.yml` checks that the Nextflow wiring holds together, and `container-test.yml` checks that each image starts. Neither runs a tool on real reads, so a step can be broken in real use while both stay green. The **E2E** workflow ([`.github/workflows/e2e.yml`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/workflows/e2e.yml)) closes that gap: it runs the real tools, in the pinned containers, on a small slice of a public genome, and checks what they write.

## The fixture

The test data is a slice of **HG002**, the Genome in a Bottle (GIAB) son of the Ashkenazi trio. HG002 is a public, consented reference sample, so nothing in the fixture is personal data.

[`scripts/ci/build-fixture.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/build-fixture.sh) builds it. It reads only the regions below from GIAB's 60x GRCh38 BAM over HTTPS (samtools fetches byte ranges through the `.bai`, so the 126 GB file is never downloaded) and samples them down to about 30x.

| Region (GRCh38) | Why it is there |
|---|---|
| chr20:10,000,000-10,500,000 | small variants; the planted ClinVar record (in SNAP25); the GIAB truth slice |
| chr22:42,000,000-42,300,000 | CYP2D6 and CYP2D7 with flanks |
| chr12:47,800,000-47,950,000 | VDR, the control gene pypgx normalises depth with |
| chr10:94,700,000-95,000,000 | CYP2C19 and CYP2C9 |
| chr6:29,900,000-33,100,000 | HLA |
| chr5:69,900,000-71,100,000 | SMN1 and SMN2 |
| chrX:73,700,000-74,000,000 | non-PAR chrX |
| chrX:1,000,000-1,200,000 | PAR1 |
| chrY:2,700,000-3,000,000 | non-PAR chrY (HG002 is male) |
| chrM | the whole mitochondrial genome |

The chr5, chrX and chrY slices are not used by every step yet. They are there so the paralog, ploidy, sample-QC and Y-haplogroup work needs no new fixture.

What the release holds:

| File | Content |
|---|---|
| `HG002_R1.fastq.gz`, `HG002_R2.fastq.gz` | the sliced reads as name-sorted pairs, input for step 02 |
| `HG002_slice.bam` (+ `.bai`) | the same reads as GIAB aligned them |
| `fixture_ref.fa.gz` (+ `.fai`, `.gzi`, `.dict`) | whole chr5, chr6, chr10, chr12, chr20, chr22, chrX, chrY and chrM from the NCBI GRCh38 no-alt analysis set, so every coordinate is real |
| `clinvar.vcf.gz`, `clinvar_chr.vcf.gz`, `clinvar_pathogenic_chr.vcf.gz` (+ `.tbi`) | ClinVar records inside the regions, built the way `setup.sh` builds the full files |
| `planted.tsv` | one synthetic ClinVar record (CLNSIG Pathogenic, gene SNAP25, ID 900000001) at a SNV HG002 is homozygous for, so step 06 always has a hit with a known gene |
| `HG002_vep.vcf` | up to 200 GIAB truth variants (HLA-A, -B, -C, CYP2C19, CYP2C9, SNAP25) annotated by VEP `--database --everything`; it stands in for step 13, whose offline cache does not fit a runner |
| `revel_synthetic.tsv.gz` (+ `.tbi`) | a score file in REVEL's layout for the SNVs of the VEP subset. The scores are made up; they only show that a score track is applied |
| `HG002_sv_manta_style.vcf.gz` (+ `.tbi`) | ten Manta-style SV records, input for AnnotSV and the SV readers |
| `HG002_truth_chr20.vcf.gz` (+ `.tbi`), `HG002_truth_chr20.bed` | GIAB v4.2.1 truth for the chr20 slice |
| `regions.bed`, `MANIFEST.txt`, `SHA256SUMS` | the slices, how this build was made (sources, depth, ClinVar date, build-script checksum), checksums |

The whole release stays under 1.5 GB. The build checks that the BAM passes `samtools quickcheck`, that every `.gz` file passes `gzip -t`, and that `samtools idxstats` shows reads on all nine contigs. The e2e job repeats those checks after download.

### Where it lives and how to change it

The fixture is published as assets of a GitHub release, a prerelease that is never marked latest. Its tag is the one line in [`tests/fixtures/VERSION`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/tests/fixtures/VERSION) (`fixture-v1`). The e2e job reads the same file, so the workflow never changes when the data does.

To change the data:

1. Edit `scripts/ci/build-fixture.sh`.
2. Bump `tests/fixtures/VERSION` (for example to `fixture-v2`).
3. Push the branch. The `build-fixture` job runs on any push that changes either file, builds the data on a GitHub runner (about 30 minutes, most of it VEP querying Ensembl's database) and publishes the new release. The e2e job of your pull request waits up to 45 minutes for it.

A push that changes the build script but keeps the old version fails on purpose: the existing release was built by different code, and replacing its files would change the data under every open pull request. To rebuild a release in place anyway (for example after a failed upload), run the E2E workflow by hand with `job: build-fixture` and `rebuild: true`.

You can build it yourself on Linux with Docker, `bgzip` and `tabix` and about 10 GB of free disk: `scripts/ci/build-fixture.sh /path/to/out`.

## The e2e job

[`scripts/ci/e2e-run.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/e2e-run.sh) downloads the fixture, checks `SHA256SUMS`, lays out a `GENOME_DIR` the way `setup.sh` and step 13 would leave it, and runs every case file in [`tests/e2e/`](https://github.com/GeiserX/Personal-Genome-Pipeline/tree/main/tests/e2e). It runs all of them even when one fails, so one run lists every broken step. The job summary shows a table of case, result, time and log, plus the failed checks of each failed case; the full logs are in the `e2e-logs` artifact.

It runs on pull requests that touch `scripts/`, `modules/`, `workflows/`, `bin/`, `conf/`, `tests/e2e/`, `tests/fixtures/VERSION`, `main.nf`, `nextflow.config`, `versions.env` or the workflow itself; once a month; and by hand. It uses a standard GitHub-hosted runner (4 CPUs, 16 GB of RAM) after deleting the preinstalled Android, .NET and Haskell toolchains to make disk room. Pulled images are cached as one compressed tar keyed on `versions.env` and the module files.

What it covers:

- **Bash steps, in order:** `validate-setup.sh`, 02 (alignment from FASTQ), 03 (DeepVariant), 03a (GATK HaplotypeCaller on the chr20 slice), 06, 07 (PharmCAT), 11, 12, 16, 16b, 20 (Mutect2 on chrM), 21 (Cyrius), 32 (pypgx with its bundle), 27, then 30, 23 and 31 on the VEP subset with the synthetic score file, then 24 and `generate-report.sh`.
- **Nextflow:** `nextflow run main.nf -profile docker` with real containers on the VCF and BAM the bash steps produced, through today's VCF+BAM samplesheet, with `--tools clinvar,mosdepth,delly,vcfanno,roh,pharmcat,cpic`, `--max_cpus 4 --max_memory 14.GB`.
- **Report assets:** a last case lists every external `http(s)` address in a `src=` or `href=` attribute of the generated HTML reports, so the remote files a report loads when opened are known. The list goes to the job summary.

Each check is a count or a column: a VCF with records, the sample name in the VCF header, a gene in the ClinVar hit, at least one called gene in PharmCAT's `report.json`, a CYP2D6 row from pypgx, no `.|.` cell in the HTML report. An exit code alone never passes a case.

Many steps ask Docker for `--cpus 8`, which Docker refuses on a 4-CPU machine. The job puts a small shim first on `PATH` ([`tests/e2e/bin/docker`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/tests/e2e/bin/docker)) that lowers `--cpus` to the CPUs the machine has and passes everything else through.

### Adding a case

A case is one executable file in `tests/e2e/`. It sources `lib.sh`, runs one step, checks its output and ends with `finish`:

```bash
#!/usr/bin/env bash
# What this case shows, in one or two lines.
. "$(dirname "$0")/lib.sh"

run_step 11-roh-analysis.sh "$SAMPLE"
check_step_exit 11-roh-analysis.sh
check_ge "per-site (ST) lines in the ROH output" \
  "$(grep -c '^ST' "${GENOME_DIR}/${SAMPLE}/vcf/${SAMPLE}_roh.txt" || true)" 100

finish
```

Cases run in name order in the C locale: the numbered base cases first, then files named `<package-key>-<what>.sh`, which can read everything the base cases wrote. Add your own file rather than editing the runner, the workflow or another package's case. `lib.sh` has the helpers: `check`, `check_eq`, `check_ge`, `has`, `lacks`, `bcf` and `sam` (bcftools and samtools from the pinned images, paths relative to `GENOME_DIR`), `vcf_count`, `planted`, `sample_side_hits` and `pharmcat_called`.

A new check has to be seen failing once before it is trusted: run it against the code before your fix, or against a deliberately broken input, and link that red run in the pull request.

### What it cannot cover

These need data or hardware a GitHub runner does not have, so no CI job runs them:

| Step | Why not |
|---|---|
| 13 (VEP, offline) | the cache is 26 GB; the fixture's VEP subset was annotated once with `--database` instead |
| 17 (CPSR) | the PCGR bundle (about 8 GB) plus a second VEP cache |
| 05 (AnnotSV) | the annotation data is 5.3 GB; not part of the default e2e run |
| 04b (GRIDSS) | needs a 31 GB Java heap |
| 02a (BWA-MEM2 index build) | about 90 GB of RAM for GRCh38 |

## The settle-doubts job

[`scripts/ci/settle-doubts.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/settle-doubts.sh) answers questions a review could not settle by reading the code: whether GATK accepts a BAM without a read group, which Clair3 models the image ships, what TIDDIT does with only a BWA-MEM2 index, whether fastp trimming changes DeepVariant's accuracy against the GIAB truth, and so on. It prints one table row per question to the job summary, with the command and what was observed. It is not a gate. Run it by hand with `job: settle-doubts`; it also runs on a push that changes the script.

## Running the e2e job yourself

The e2e run needs Linux, Docker, Nextflow 25.10 with Java 17, the `gh` CLI and about 25 GB of disk, and takes most of an hour. Run it on a machine you do not need for anything else:

```bash
E2E_WORK=/path/with/space scripts/ci/e2e-run.sh          # every case
E2E_WORK=/path/with/space scripts/ci/e2e-run.sh '3*'     # only cases whose file name matches
```

A partial run is for debugging: later cases read what earlier cases wrote.
