# Testing

This page is for contributors. It says what CI runs, how the end-to-end test works, and how to change its data.

The stub runs in `nextflow.yml` check that the Nextflow wiring holds together, and `container-test.yml` runs each changed image's own tool on the fixture. Neither runs a pipeline step as the scripts and modules call it, so a step can be broken in real use while both stay green. The **E2E** workflow ([`.github/workflows/e2e.yml`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/workflows/e2e.yml)) closes that gap: it runs the real tools, in the pinned containers, on a small slice of a public genome, and checks what they write.

## The fixture

The test data is a slice of **HG002**, the Genome in a Bottle (GIAB) son of the Ashkenazi trio. HG002 is a public, consented reference sample, so nothing in the fixture is personal data.

[`scripts/ci/build-fixture.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/build-fixture.sh) builds it. It reads only the regions below from GIAB's 60x GRCh38 BAM on GIAB's S3 mirror (`https://giab.s3.amazonaws.com/`) over HTTPS (samtools fetches byte ranges through the `.bai`, so the 126 GB file is never downloaded) and samples them down to about 30x.

| Region (GRCh38) | Why it is there |
|---|---|
| chr1:109,600,000-109,800,000 | GSTM1 (pypgx reads depth over it) |
| chr2:233,600,000-233,800,000 | UGT1A1 and UGT1A4 |
| chr4:68,500,000-68,700,000 | UGT2B15 and UGT2B17 |
| chr16:28,580,000-28,630,000 | SULT1A1 |
| chr19:40,800,000-41,050,000 | CYP2A6, CYP2A7 and CYP2B6 |
| chr20:10,000,000-10,500,000 | small variants; the planted ClinVar record (in SNAP25); the GIAB truth slice |
| chr22:42,000,000-42,300,000 | CYP2D6 and CYP2D7 with flanks |
| chr12:47,800,000-47,950,000 | VDR, the control gene pypgx normalises depth with |
| chr10:94,700,000-95,000,000 | CYP2C19 and CYP2C9 |
| chr6:29,900,000-33,100,000 | HLA |
| chr5:69,900,000-71,100,000 | SMN1 and SMN2 |
| chrX:73,700,000-74,000,000 | non-PAR chrX |
| chrX:1,000,000-1,200,000 | PAR1 |
| chrY:2,700,000-3,000,000 | non-PAR chrY (HG002 is male) |
| chrM | the whole mitochondrial genome, sampled to about 500x (the source has thousands of x there) |

The chr5, chrX and chrY slices are not used by every step yet. They are there so the paralog, ploidy, sample-QC and Y-haplogroup work needs no new fixture.

What the release holds:

| File | Content |
|---|---|
| `HG002_R1.fastq.gz`, `HG002_R2.fastq.gz` | the sliced reads as name-sorted pairs, input for step 02 |
| `HG002_slice.bam` (+ `.bai`) | the same reads as GIAB aligned them |
| `HG002_cyrius.bam` (+ `.bai`) | GIAB's alignment of the regions Cyrius (step 21) reads, sampled the same way: CYP2D6, CYP2D7 and its 3,000 depth-normalisation bins on chr1 to chr22, plus the two 50 kb flanks (chr22:42.05-42.10 and 42.20-42.25 Mb) step 21's depth check compares CYP2D6 with. Cyrius stops on the first bin whose contig the BAM lacks, so it cannot run on the BAM step 02 writes against the fixture reference |
| `fixture_ref.fa.gz` (+ `.fai`, `.gzi`, `.dict`) | whole chr1, chr2, chr4, chr5, chr6, chr10, chr12, chr16, chr19, chr20, chr22, chrX, chrY and chrM from the NCBI GRCh38 no-ALT analysis set, the pipeline's default reference, so every coordinate and contig name is real |
| `clinvar.vcf.gz`, `clinvar_chr.vcf.gz`, `clinvar_pathogenic_chr.vcf.gz` (+ `.tbi`) | ClinVar records inside the regions, built the way `setup.sh` builds the full files |
| `planted.tsv` | one synthetic ClinVar record (CLNSIG Pathogenic, gene SNAP25, ID 900000001) at a SNV HG002 is homozygous for, so step 06 always has a hit with a known gene |
| `HG002_vep.vcf` | up to 200 GIAB truth variants (HLA-A, -B, -C, CYP2C19, CYP2C9, SNAP25) annotated by VEP `--database --everything`; it stands in for step 13, whose offline cache does not fit a runner. `--database` gives no gnomAD frequencies, so steps 23 and 31 only run their no-gnomAD branch |
| `revel_synthetic.tsv.gz` (+ `.tbi`) | a score file in REVEL's layout for the SNVs of the VEP subset. The scores are made up; they only show that a score track is applied |
| `HG002_sv_manta_style.vcf.gz` (+ `.tbi`) | ten Manta-style SV records, input for AnnotSV and the SV readers |
| `HG002_truth_chr20.vcf.gz` (+ `.tbi`), `HG002_truth_chr20.bed` | GIAB v4.2.1 truth for the chr20 slice |
| `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (+ `.tbi`), `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`, `HG002_GRCh38_v5.0q_smvar.vcf.gz` (+ `.tbi`), `HG002_GRCh38_v5.0q_smvar.benchmark.bed` | GIAB's two HG002 truth sets, whole and byte for byte as GIAB publishes them, for `benchmark-variants.sh --giab` (its md5 check accepts nothing else). The GIAB e2e case reads them from here, so it needs no download. About 210 MB. The build checks each file against the md5 the step accepts; MANIFEST.txt lists the md5s and sources |
| `regions.bed`, `MANIFEST.txt`, `SHA256SUMS` | the slices, how this build was made (sources, depth, ClinVar date, build-script checksum), checksums |

Why chr2, chr4 and chr16: pypgx (step 32) reads depth over the region of every gene it can call copy number for before it calls any of them, and samtools refuses a region on a contig the BAM does not have. One missing contig and step 32 calls no BAM-based gene at all, CYP2D6 included. Its GRCh38 region for GSTT1 is on an ALT contig, `chr22_KI270879v1_alt`, which the default reference does not have; step 32 leaves GSTT1 out on the fixture as it does on a real sample. Up to `fixture-v4` the fixture carried that contig; `fixture-v5` dropped it, so the e2e job runs on the contigs a default setup has.

The whole release stays under 1.5 GB. The build checks that both BAMs pass `samtools quickcheck`, that every `.gz` file passes `gzip -t`, that `samtools idxstats` shows reads on every primary contig of the reference, that the reference has no ALT or HLA contig, that the Cyrius BAM has reads on chr1 to chr22 and reads with MAPQ >= 1 in both CYP2D6 flanks, and that each GIAB truth file has the md5 `benchmark-variants.sh` accepts. After download the e2e job repeats the BAM, `gzip -t`, contig and size checks (`tests/e2e/01-fixture.sh`); the GIAB case checks the truth files by md5 and `pgx-1-depth-check.sh` checks the flank depth.

### Where it lives and how to change it

The fixture is published as assets of a GitHub release, a prerelease that is never marked latest. Its tag is the one line in [`tests/fixtures/VERSION`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/tests/fixtures/VERSION) (now `fixture-v6`). The e2e job reads the same file, so the workflow never changes when the data does.

To change the data:

1. Edit `scripts/ci/build-fixture.sh`.
2. Bump `tests/fixtures/VERSION` (for example from `fixture-v6` to `fixture-v7`).
3. Push the branch. The `build-fixture` job runs on any push that changes either file, builds the data on a GitHub runner (30 to 70 minutes, most of it VEP querying Ensembl's public database) and publishes the new release. The e2e job of your pull request waits up to 75 minutes for it.

A push that changes the build script but keeps the old version fails on purpose: the existing release was built by different code, and replacing its files would change the data under every open pull request. To rebuild a release in place anyway (for example after a failed upload), run the E2E workflow by hand with `job: build-fixture` and `rebuild: true`.

You can build it yourself on Linux with Docker, `bgzip` and `tabix` and about 10 GB of free disk: `scripts/ci/build-fixture.sh /path/to/out`.

## The e2e job

[`scripts/ci/e2e-run.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/e2e-run.sh) downloads the fixture, checks `SHA256SUMS`, lays out a `GENOME_DIR` the way `setup.sh` and step 13 would leave it, and runs every case file in [`tests/e2e/`](https://github.com/GeiserX/Personal-Genome-Pipeline/tree/main/tests/e2e). It runs all of them even when one fails, so one run lists every broken step. The job summary shows a table of case, result, time and log, plus the failed checks of each failed case; the full logs are in the `e2e-logs` artifact.

It runs on pull requests that touch `scripts/`, `modules/`, `workflows/`, `bin/`, `conf/`, `assets/`, `tests/e2e/`, `tests/fixtures/VERSION`, `main.nf`, `nextflow.config`, `versions.env` or the workflow itself; once a month; and by hand (`job: e2e`). Every pull request starts the workflow, and its `e2e-scope` job checks those paths: on a pull request that touches none of them the e2e job is skipped, which GitHub counts as passed, so main can require the e2e check without blocking a docs-only change. It takes about 75 minutes (76 on the run that added `run-all-launcher.sh`, the third run from FASTQ after the bash steps and the pipeline, about 12 minutes of it); the job's limit is 150 minutes. The DeepVariant case passes the fixture slices as `INTERVALS`, so DeepVariant calls only those slices instead of walking all 1.8 Gb of the reference, which alone took 27 to 48 minutes. It uses a standard GitHub-hosted runner (4 CPUs, 16 GB of RAM) after deleting preinstalled toolchains it does not use (Android, .NET, Haskell, CodeQL, Boost) to make disk room. Pulled images are cached as one compressed tar keyed on `versions.env` and the module files.

What it covers:

- **Bash steps, in order:** `validate-setup.sh`, 02 (alignment from FASTQ), 03 (DeepVariant), 03a (GATK HaplotypeCaller on the chr20 slice), 06, 07 (PharmCAT), 11, 12, 16, 16b, 20 (Mutect2 on chrM), 21 (Cyrius, on `HG002_cyrius.bam`), 32 (pypgx with its bundle), 27, then 30, 23 and 31 on the VEP subset with the synthetic score file, then 24 and `generate-report.sh`.
- **Nextflow:** `nextflow run main.nf -profile docker` with real containers on the VCF and BAM the bash steps produced, through today's VCF+BAM samplesheet, with `--tools clinvar,mosdepth,delly,vcfanno,roh,pharmcat,cpic`, `--max_cpus 4 --max_memory 14.GB`.
- **Parity and the launcher:** `nextflow-from-fastq-*.sh` run the bash steps 01b, 02, 03, 06, 07, 11 and 25 and the pipeline on the same FASTQ pair, and [`scripts/ci/parity-diff.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/parity-diff.sh) compares what they wrote ([Bash vs Nextflow parity](nextflow.md#bash-vs-nextflow-parity)). `run-all-launcher.sh` runs `scripts/run-all.sh` from the same reads in a `GENOME_DIR` of its own: it must publish the same files as the direct `nextflow run`, pass the same parity check against the bash steps, write both reports, and on a second run take every task from the cache in under ten minutes. Each comparison has a negative control in the same run.
- **Manta and the Nextflow hardening** (`nextflow-hardening-*.sh`): step 04 on a planted Manta-style VCF with two inversion breakend pairs (they must come out as two `SVTYPE=INV` records) and on the fixture BAM inside the fixture's regions (`MANTA_CALL_REGIONS`); then a Nextflow run with `--tools clinvar,pharmcat,mosdepth,delly,manta,vcfanno,pypgx,telomere_hunter,mito_haplogroup`, which checks that every container ran with `--network none` and that no task built a `.fai` of the reference.
- **Report assets:** a last case lists every external `http(s)` address in a `src=` or `href=` attribute of the generated HTML reports, so the remote files a report loads when opened are known. The list goes to the job summary.

Each check is a count or a column: a VCF with records, the sample name in the VCF header, a gene in the ClinVar hit, at least one called gene in PharmCAT's `report.json`, a CYP2D6 row from pypgx, the planted ClinVar row in the HTML report with its gene, its significance and no empty cell. An exit code alone never passes a case.

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

The e2e run needs Linux, Docker, Nextflow 25.10.8 with Java 17, the `gh` CLI and about 25 GB of disk, and takes over an hour. Run it on a machine you do not need for anything else:

```bash
E2E_WORK=/path/with/space scripts/ci/e2e-run.sh          # every case
E2E_WORK=/path/with/space scripts/ci/e2e-run.sh '3*'     # only cases whose file name matches
```

A partial run is for debugging: later cases read what earlier cases wrote.
