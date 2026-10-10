# Nextflow Execution

The pipeline is a [Nextflow](https://www.nextflow.io/) DSL2 pipeline, `main.nf`. Each samplesheet row starts from one of three places:

- **FASTQ**: the reads are trimmed (fastp), aligned (minimap2, read group, duplicates marked) and called (DeepVariant, a VCF and a gVCF, chrX and chrY haploid for a male sample);
- **a BAM or a CRAM** without a VCF: it is called;
- **a VCF** from any caller (nf-core/sarek, DRAGEN, a provider), with an optional BAM or CRAM.

Every BAM then goes through a sex check (indexcov), and the pipeline runs pharmacogenomics, variant annotation, clinical screening, BAM analyses, structural variant calling and reporting: 7 workflows, 62 processes in 39 module files under `modules/local/`. A VCF given in the samplesheet needs FILTER=PASS records and GRCh38 contig names with chr; see [FILTER=PASS required](#filterpass-required) and [Contig names and gVCF input](#contig-names-and-gvcf-input). Starting from a provider's VCF: [Starting from a Vendor VCF](vcf-first.md).

> **Nextflow is the pipeline; the scripts are single steps.** Each numbered script in `scripts/` runs one step on its own and takes its image tags from the same `versions.env` and its helpers from `scripts/lib/common.sh`; `run-all.sh` is a launcher: it writes a one-row samplesheet and starts this pipeline with `-resume`, then runs the few script-only steps you ask for and the reports. CI runs the scripts and the pipeline on the same reads and fails when their results differ (see [Bash vs Nextflow parity](#bash-vs-nextflow-parity)). The Singularity profile is untested (see [Profiles](#profiles)).

---

## Quick Start

### Prerequisites

1. **Docker** (already required for the bash pipeline)
2. **Java 17 or later** (Nextflow 26.04 runtime requirement; CI runs Java 17)
3. **Nextflow 26.04.7**, the version CI validates (`NEXTFLOW_VERSION` in versions.env). Pin it when installing, because the plain installer fetches the newest release:
   ```bash
   curl -s https://get.nextflow.io | NXF_VER=26.04.7 bash
   sudo mv nextflow /usr/local/bin/
   ```

### Run the Pipeline

```bash
# 1. Create a samplesheet CSV: one row per sample, from FASTQ, a BAM, or a VCF
cat > samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,bam,bam_index,vcf,vcf_index,sex
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,,,,,female
sample2,,,/path/to/sample2_sorted.bam,/path/to/sample2_sorted.bam.bai,/path/to/sample2.vcf.gz,/path/to/sample2.vcf.gz.tbi,male
EOF

# 2. Run (default tools need no external databases; prs and vcfanno are
#    skipped with a warning until --pgs_scoring or a score file is set)
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /path/to/GRCh38_no_alt_analysis_set.fasta \
    --outdir ./results \
    -profile docker

# 3. To enable database-requiring tools, add them to --tools with their flags:
#    --tools '...,vep,slivar,clinical_filter'  + --vep_cache /path/to/vep_cache
#    --tools '...,cpsr'                        + --pcgr_data + --vep_cache_cpsr
#    --tools '...,clinvar'                     + --clinvar + --clinvar_index
#    --tools '...,expansion_hunter'            + --expansion_catalog (and a sex column)
#    --tools '...,hla_typing'                  + --hla_dat + --hla_genes (hla.dat and the GENCODE gene lines setup.sh installs)
#    --tools '...,annotsv'                     + --annotsv_annotations
#    --tools '...,cnvpytor'                    + --cnvpytor_resources
#    --tools '...,delly'                       (optional --delly_exclude <excl.tsv>, passed as delly sr -x)
#    --tools '...,manta'                       (optional --manta_call_regions <regions.bed.gz>, its .tbi beside it)
#    --tools '...,sample_qc'                   + --somalier_sites + --verifybamid2_panel (setup.sh --sample-qc-data
#                                                installs both; see docs/33-sample-qc.md)
#    --tools '...,cram_archive'                writes a checked CRAM beside each BAM (docs/34-cram-archive.md)
#    --tools '...,cyrius'                      + --cyrius_install (setup.sh --cyrius; non-commercial licence)
#    --tools '...,parascopy'                   + --parascopy_data (setup.sh --parascopy-data; docs/35-paralogs.md)
#    --kir true (with hla_typing)              + --kir_dat (setup.sh --kir-data; KIR genes, a second T1K pass)
#    With pharmcat, PGX_CONSENSUS gives PharmCAT the HLA types of hla_typing and a CYP2D6 call only when
#    pypgx and cyrius agree (docs/36-pgx-consensus.md); PharmCAT then waits for those BAM steps.
#    telomere_hunter (a default tool) takes --cytoband <cytoBand.hg38.txt> for GRCh38 bands; without it,
#    TelomereHunter uses its own hg19 bands and the run logs a warning
#    An unknown name in --tools stops the run.
```

`--reference` is the FASTA the reads are aligned to and the BAMs were aligned to, the GRCh38 no-ALT analysis set that `setup.sh` installs as `reference/GRCh38_no_alt_analysis_set.fasta`. The Nextflow pipeline does not compare a given BAM's header with it; `./scripts/validate-setup.sh <sample>` does, and [Realigning after a reference change](realignment.md) covers a BAM aligned to another reference.

### From reads

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_trim` | false | Align the raw reads, without fastp (`SKIP_TRIM=true` of step 01b) |
| `--minimap2_index` | `<reference base>.sr.mmi` beside the FASTA | The minimap2 index built with `-x sr`, the file step 02 builds. When it is not there, the run builds one (about 30 minutes for GRCh38, cached by `-resume`) |
| `--intervals` | whole genome | Space-separated regions DeepVariant calls (`INTERVALS` of step 03) |
| `--sex_check` | `fail` | When the sex indexcov infers from a BAM differs from the samplesheet's: `fail` stops the run before any BAM step starts, `warn` logs both and goes on with the declared sex. With `sample_qc`, the same applies to the sex somalier infers from the reads |

The trimmed reads stay in the work directory; the fastp reports are published. A run from FASTQ needs the same memory as step 02 for alignment (minimap2 peaks near 10 GB on a 1.8 Gb reference) and DeepVariant's (32 GB by default); `--max_memory` caps both.

### Resume After Failure

`-resume` reuses every task that finished. If a step fails, fix the issue and resume:

```bash
nextflow run main.nf -resume [same params as before]
```

The failed tasks and the ones downstream of them run again. A task that was running when the run stopped starts over: DeepVariant is a single task, so a run stopped during it loses that task's progress. Run the pipeline inside tmux or screen, or with nohup, so closing the terminal does not stop it. Nextflow keys an input file by its path, size and modification time, not by its content: a file rewritten in place with the same size and time is not seen as changed.

One exception, once: Nextflow 26.04 hashes a map input by its keys and values where 25.10 hashed it differently, so the first `-resume` of a work directory written by 25.10 reruns every process that takes the `meta` map, which is every process here. Start a fresh work directory for the first run after the upgrade instead of waiting on a resume that caches nothing.

---

## Samplesheet Format

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier (used as output directory name) |
| `fastq_1`, `fastq_2` | One of three* | Paired gzipped FASTQ: trimmed, aligned and called |
| `bam`, `bam_index` | One of three* | Aligned BAM and its index (`.bam.bai`): called when the row has no VCF |
| `cram`, `crai` | Instead of a BAM | Aligned CRAM and its index, read with `--reference` (the FASTA it was written with). `CRAM_TO_BAM` writes it out as a BAM in the work directory, checked against it, and the row goes on as a BAM row |
| `vcf`, `vcf_index` | One of three* | Bgzipped VCF (`.vcf.gz`) and its tabix index from any caller; a BAM on the same row is optional |
| `gvcf`, `gvcf_index` | No | The gVCF of the row's VCF (step 03 writes one) and its index. PharmCAT and PRS read it as they read the gVCF DEEPVARIANT writes; the row is not called again. Needs `vcf` on the row; `run-all.sh` fills it from `vcf/<sample>.g.vcf.gz` |
| `sex` | On called rows** | `male` or `female` |

\* A row starts from FASTQ, or from a BAM or CRAM, or from a VCF (with or without a BAM or CRAM); a row with FASTQ and a BAM, CRAM or VCF stops the run. A VCF-only row is valid for annotation and PGx, but most default tools (mosdepth, telomere_hunter, mito_variants) and opt-in tools (expansion_hunter, hla_typing, pypgx, cyrius, parascopy) need a BAM. **Provide reads or a BAM for full analysis.** A row the pipeline calls gets a VCF and a gVCF, and a VCF row can give its gVCF in the `gvcf` column; PharmCAT and PRS then read the sites where the sample matches the reference from the gVCF.

The VCF must name its contigs the GRCh38 way with chr (`chr1` to `chr22`, `chrX`, `chrY`, `chrM`); a VCF named `1`, `MT` stops the run with the rename command. A gVCF given in the `vcf` column stops the run with `pharmcat` selected: the pipeline expands the reference blocks of the gVCF DeepVariant writes for a called row, or of the one in the `gvcf` column, not of a gVCF given as the VCF, and a variants-only VCF leaves about half of PharmCAT's genes Unknown. [Starting from a Vendor VCF](vcf-first.md) has the commands for both.

\*\* `sex` is required on every row the pipeline calls (FASTQ, or a BAM or CRAM without a VCF): for a male sample DeepVariant calls chrX and chrY haploid outside the pseudoautosomal regions. It is required on every row with a BAM or CRAM when `expansion_hunter` is in `--tools`, where it sets the chrX ploidy (ExpansionHunter's default is female). A row that needs it and lacks it stops the run at parse time.

**Sex check.** With `sample_qc` in `--tools`, somalier also infers the sex from the reads (chrX heterozygosity) and the same rule and `--sex_check` apply; see [Step 33](33-sample-qc.md). `INDEXCOV` (goleft indexcov, seconds per sample: it reads only the `.bai`) infers each BAM's sex from the chrX and chrY copy numbers and writes it to `<sample>_sex_check.tsv`. When the row declares a sex and indexcov infers another (or cannot tell), the run stops before any step reads the BAM, and the message gives both values and the copy numbers: the sample is not the one you think, the declared sex is wrong, or the sample has a sex-chromosome aneuploidy. `--sex_check warn` logs it and goes on with the declared sex. On a small region slice, like the test fixture, indexcov's call is not reliable.

Each `sample` value must appear once; a repeated id stops the run, because the id names the output directory and keys every per-sample join.

### Using Sarek Output

If you ran [nf-core/sarek](https://nf-co.re/sarek) for alignment and variant calling, point the samplesheet at sarek's output files. Sarek 3.x writes CRAM by default: give it in the `cram` and `crai` columns, with sarek's reference as `--reference`, or run sarek with `--save_output_as_bam` to get the `.recal.bam` files below:

```csv
sample,vcf,vcf_index,bam,bam_index
sample1,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz.tbi,results/preprocessing/recalibrated/sample1/sample1.recal.bam,results/preprocessing/recalibrated/sample1/sample1.recal.bam.bai
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (default for local). Every container runs with `--network none`. |
| `singularity` | Singularity/Apptainer. **Untested.** One module writes inside its image and needs a writable container: `cnvpytor` (copies resources into its `site-packages`). This profile does not cut the network. |
| `test` | Minimal test with reduced resources |
| `test_full` | Full-size test with real WGS data |

Combine profiles: `-profile docker,test`

---

## Resource Configuration

Default resource limits (tuned for 16-core consumer desktop). Each process label asks for a fixed CPU count, a memory and a time limit; after an exit that usually means a kill for memory (137 and the others `conf/base.config` lists), the task is retried once with memory and time doubled. On Nextflow 26.04.7, a task that runs past its time limit can stop the whole run instead of being retried ([nextflow-io/nextflow#7569](https://github.com/nextflow-io/nextflow/issues/7569), fixed upstream after 26.04.7); we read this in the source and have not reproduced it. These limits cap each request:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | 16 | Maximum CPUs per process |
| `--max_memory` | 64.GB | Maximum memory per process |
| `--max_time` | 48.h | Maximum wall time per process |

`--max_cpus` and `--max_memory` cap each task. They do not limit how much runs at once: Nextflow fills the machine's CPUs and RAM with as many tasks as fit.

Override for smaller machines, with `--max_memory` below the RAM the machine reports. On a 32 GB Linux machine the kernel reports a little less than 32 GiB, so a 32 GB cap leaves every 32 GB task (DeepVariant, VEP and the other `process_high` steps) larger than the machine, and Nextflow refuses to start it:

```bash
nextflow run main.nf --max_cpus 8 --max_memory 30.GB [other params]
```

---

## Output Structure

```
results/
├── sample1/
│   ├── fastq_trimmed/      # fastp reports, JSON + HTML (FASTQ rows)
│   ├── aligned/            # <sample>_sorted.bam + .bai: minimap2, duplicates marked (FASTQ rows);
│   │                       #   <sample>_sorted.cram + .crai, checked against the BAM (cram_archive)
│   ├── vcf/                # <sample>.vcf.gz and <sample>.g.vcf.gz + .tbi: DeepVariant (called rows)
│   ├── indexcov/           # goleft indexcov coverage plots and .ped (rows with a BAM)
│   ├── <sample>_sex_check.tsv  # the sex indexcov infers, with CNchrX and CNchrY
│   ├── qc/                 # <sample>_sample_qc.tsv: somalier's sex, FREEMIX, same person as (sample_qc);
│   │                       #   somalier/ and verifybamid2/ hold each tool's own files
│   ├── pharmcat/           # PharmCAT PGx reports (HTML + JSON)
│   ├── clinvar/            # ClinVar pathogenic variant screen: hits as VCF and TSV
│   ├── pypgx/              # pypgx star allele calling (optional)
│   ├── cpic/               # CPIC drug-gene recommendations (optional)
│   ├── vep/                # VEP VCF, and the vcfanno-enriched VCF when a score file is set
│   ├── slivar/             # Prioritized variants + compound hets
│   ├── clinical/           # Clinically relevant variant subset
│   ├── cpsr/               # Cancer predisposition report
│   ├── roh/                # Runs of homozygosity
│   ├── prs/                # Polygenic risk scores
│   ├── ancestry/           # Ancestry: projection onto the panel (with --ancestry_ref)
│   ├── mito/               # Mitochondrial variant calls, the haplogroup from them and haplocheck's contamination check
│   ├── y_haplogroup/       # Y-chromosome haplogroup of male samples (Yleaf, opt-in: y_haplogroup)
│   ├── hla/                # HLA typing
│   ├── expansion_hunter/   # Repeat expansion calls
│   ├── telomere/           # Telomere length estimation
│   ├── coverage/           # Coverage statistics (mosdepth)
│   ├── cyrius/             # CYP2D6 star allele (Cyrius, opt-in) and the CYP2D6 depth check
│   ├── pgx_consensus/      # PharmCAT's outside calls (HLA-A/B, an agreed CYP2D6) and who said what
│   ├── kir/                # KIR genotypes and the IPD-KIR release (--kir)
│   ├── paralogs/           # SMN1/SMN2 copy number (Parascopy, opt-in)
│   ├── manta/              # SV calling (optional): diploidSV.vcf.gz with inversions as SVTYPE=INV;
│   │                       #   diploidSV.raw.vcf.gz as Manta wrote it (an inversion is two BND records)
│   ├── sv_duphold/         # Manta SVs with duphold depth tags (optional)
│   ├── sv_filtered/        # Manta SVs after the duphold depth filter (optional)
│   ├── annotsv/            # AnnotSV ACMG classification of the filtered SVs (optional)
│   ├── delly/              # SV calling (optional)
│   ├── cnvpytor/           # CNV calling (optional)
│   ├── sv_merged/          # SV consensus of two or more callers, SURVIVOR merge (optional)
│   ├── summary.json        # The numbers the report is rendered from (html_report)
│   └── *_report.html       # Summary HTML report (published to sample root): QC, ClinVar, PharmCAT, CPIC,
│                           #   CPSR, clinical filter (with its ACMG SF tier), slivar, ROH, mito haplogroup and haplocheck,
│                           #   Y haplogroup; "Not run" for a tool not selected
├── somalier/               # somalier relate over every sample of the run: samples, pairs, HTML (sample_qc)
├── multiqc/                # MultiQC report across samples (reads mosdepth: needs a BAM; a VCF-only run logs the skip)
└── pipeline_info/
    ├── timeline_*.html
    ├── report_*.html
    ├── trace_*.txt
    └── dag_*.svg
```

---

## Before you share outputs

The pipeline does not anonymise anything: the outputs carry whatever identified you in the input. Before you share a file, know what it holds:

- **The VCF's sample name** (the last column of its `#CHROM` line) is repeated in the ClinVar VCFs (`clinvar/<sample>_clinvar_hits.vcf`, `<sample>_pass.vcf.gz`), on every line of `roh/<sample>_roh.txt`, in `mito/<sample>_haplogroup.txt` and as `sampleId` in the PharmCAT JSON files.
- **The input's header lines** pass through into the ClinVar VCFs, including the provider's and bcftools' command lines, which often name the sample or a file.
- **A BAM the pipeline aligned** names the sample label in its `@RG` line, and its `@PG` lines hold the command lines with your file paths.
- **The input file name** is in the `##bcftools_viewCommand` and `##bcftools_normCommand` lines of the ClinVar VCFs, and can be in the command lines other tools print into their outputs.
- **`pipeline_info/`** (report, timeline, trace) holds absolute paths of your machine.
- **The samplesheet's `sample` label** names every output folder and file.

A neutral input removes the metadata that names you. Use a label that does not name you, name the file after it, and keep only the header lines the tools read. The recipe below lists the lines to keep and renames the sample column. A list of keys to delete would miss lines, because callers name them differently.

It does not anonymise the data. Genotype files (VCF, gVCF, BAM, CRAM) can identify you and your relatives through genealogy databases: Erlich et al. 2018 (Science 362:690) projected a third-cousin-or-closer match for about 60% of searches for people of European descent, and Gymrek et al. 2013 (Science 339:321) recovered surnames from Y-chromosome STRs.

```bash
bcftools view --no-version -h in.vcf.gz | grep -E '^(##(fileformat|FILTER|INFO|FORMAT|ALT|contig)=|#CHROM)' > h.txt
echo SAMPLE > names.txt
bcftools reheader -h h.txt -s names.txt -o SAMPLE.vcf.gz in.vcf.gz
bcftools index -t SAMPLE.vcf.gz
```

`bcftools annotate -x` cannot remove these lines (it exits with "No matching tag"). The same recipe, with bcftools from the pinned image, is step 4 of [Starting from a Vendor VCF](vcf-first.md). Do not share `pipeline_info/`.

---

## The pipeline and the single-step scripts

Run the pipeline for a sample: it runs every step a default run selects, in parallel where the inputs allow, and `-resume` reruns only what changed. Run a numbered script to run or rerun one step by hand, on the layout `run-all.sh` and `setup.sh` use; the scripts and the modules run the same commands from the same images (see the parity section below). The pipeline needs Java 17 or later and Nextflow beside Docker; the scripts need Docker only. The HPC/Singularity profile exists but is untested.

---

## Known Limitations & Design Decisions

### Scope: from reads to report

The pipeline starts from paired short-read FASTQ, from a BAM or from a VCF. The reads go through fastp, minimap2 with the `sr` preset (`samtools fixmate`, `sort`, `markdup`) and DeepVariant's WGS model, with the flags of steps 01b, 02 and 03. Other starting points are scripts (below): ORA files, long reads, other aligners and callers, chip data.

### Bash vs Nextflow parity

The pipeline and the scripts run the same commands from the same images, so on the same input they must give the same answers. The E2E workflow checks this on the test fixture: the scripts 01b, 02, 03, 06, 07, 11 and 25, and the pipeline from the same FASTQ pair, then [`scripts/ci/parity-diff.sh`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/scripts/ci/parity-diff.sh) compares, item by item:

| Item | Compared |
|---|---|
| alignment | `samtools flagstat` of the two BAMs |
| variants | `bcftools isec` of the two VCFs: records in one only, and shared records with another FILTER or GT |
| gvcf | every gVCF record |
| clinvar | the rows of `<sample>_clinvar_hits.tsv` |
| roh | the `bcftools roh` segments and the 5 Mb summary |
| pharmcat | the diplotype of each gene |
| prs | each score's sum, matched variant count and input (gVCF or VCF) |
| sv_merged | each SV consensus record (position, end, type, length, support and each caller's GT); on the Nextflow hardening run, against steps 22 and 37 on the same calls and BAM |
| y_haplogroup | Yleaf's prediction table; on the same hardening run |

A difference fails the check unless the script lists it with its reason; none is listed today. Output file names and folders differ (the scripts write under `$GENOME_DIR/<sample>/`, the pipeline under `<outdir>/<sample>/`, mapped in the script).

Where a module and its script differ on purpose:

| Step | Script | Module |
|---|---|---|
| Alignment (02) | pipes minimap2 into samtools, each in its own image | `ALIGN_MINIMAP2` writes the SAM compressed with `gzip -1` and `ALIGN_MARKDUP` reads it: a task runs in one image, and no image in `versions.env` holds both tools. Same commands and the same BAM, plus a temporary file in the work directory, larger than the gzipped FASTQ: one 30x run wrote 119 GB |
| Sex check (16) | runs beside the other steps and stops itself on a mismatch | `INDEXCOV` runs before every BAM step, and a mismatch stops the run before DeepVariant starts |
| DeepVariant (03) | `MODEL_TYPE` picks WGS, WES, PACBIO or ONT_R104 | the WGS model only: the pipeline takes paired short reads |
| HLA typing (08) | keeps the T1K index under `t1k_idx/`, named after the T1K version, the IPD-IMGT/HLA release and the GENCODE release | `T1K_BUILD` builds it once per run for every sample; the task hash covers the same three, and `-resume` reuses it |
| PRS (25) | scores the list in `assets/pgs_scores.tsv`, downloading a missing file | scores every file `--pgs_scoring` holds, labelled from `assets/pgs_scores.tsv` (an id not in it is labelled with its file's `trait_reported`); `PRS` runs pgsc_calc on the host, see [No network inside the containers](#no-network-inside-the-containers) |
| ExpansionHunter (09) | uses the GRCh38 catalog inside the image, or `EH_CATALOG` | needs `--expansion_catalog` |
| Mito variants (20) | extracts the chrM reads with `samtools view` | with GATK `PrintReads`. Mutect2 applies its own read filters to either, so the calls are expected to match; CI does not compare them. Both mark possible NuMTs at the median autosomal depth from mosdepth (step 16b, or `MOSDEPTH` when `mosdepth` is in `--tools`) |
| AnnotSV (05) | annotates step 15's duphold-filtered calls; without them, or when Manta's calls are newer, Manta's calls, and says so | `ANNOTSV` always reads `DUPHOLD_FILTER`'s output: `annotsv` in `--tools` needs `duphold` and `manta` |
| Y haplogroup (37) | reads the sex step 16 (indexcov) infers | runs on the rows whose samplesheet sex is male, which `INDEXCOV` has checked against the reads |
| HTML report (24) | renders every section from `bin/collect_summary.py`'s summary | `HTML_REPORT` runs the same code on the outputs of this run's QC, ClinVar, PharmCAT, CPIC, CPSR, clinical filter, slivar, ROH, mito haplogroup, haplocheck and Y haplogroup steps; the clinical filter and slivar cards show counts only (the module gets their VCFs, not their tables). For every section, run `GENOME_DIR=<outdir> scripts/24-html-report.sh <sample>` on the Nextflow output |
| CNVpytor (18) | mounts each resource file over the image's data folder | copies the files into the image's `site-packages`, so it needs a writable container |
| Cyrius (21) | runs Cyrius from the install `setup.sh --cyrius` made | the same install, given as `--cyrius_install` (see [Cyrius, opt-in](#cyrius-opt-in)) |

Scripts with no module, and why:

| Script | Why it stays bash-only |
|---|---|
| `01-ora-to-fastq.sh` | `orad` is a native binary that runs on the host, on the raw reads |
| `02a-alignment-bwamem2.sh`, `03a-gatk-haplotypecaller.sh`, `03b-freebayes.sh`, `03c-strelka2-germline.sh`, `03d-octopus.sh`, `03e-clair3.sh` | alternative and legacy aligners and callers, kept for benchmarking: the pipeline aligns with minimap2 and calls with DeepVariant, and takes the VCF of any other caller in the `vcf` column |
| `02b-alignment-longread.sh`, `04c-sniffles2.sh` | the long-read path: the samplesheet holds one short-read BAM per sample |
| `04a-tiddit.sh` | an alternative to Manta; its local assembly needs the classic BWA index |
| `04b-gridss.sh` | needs a 31 GB Java heap and the BWA index, more than the SV modules are sized for |
| `14-imputation-prep.sh` | writes per-chromosome files for an upload to an imputation server, a manual step |
| `29-mutect2-somatic.sh` | experimental tumor-only somatic calling |
| `chip-to-vcf.sh` | turns consumer array data into a VCF: an input for the pipeline, not a step on one |
| `benchmark-variants.sh` | compares callers with each other or with a truth set: a check of the calling, not an analysis of the sample |
| `generate-report.sh` | the text report. Run `GENOME_DIR=<outdir> scripts/generate-report.sh <sample>` on the Nextflow output |
| `run-all.sh`, `setup.sh`, `validate-setup.sh` | running, installing and checking the bash path |

### Reference databases not auto-downloaded

Several tools require large reference databases that are **not automatically downloaded** by the pipeline. You must obtain and provide paths for these yourself:

| Parameter | Required by | Size |
|-----------|------------|------|
| `--vep_cache` | VEP annotation | ~28 GB download (release 116) |
| `--pcgr_data` | CPSR cancer predisposition | ~7 GB bundle (20260620) |
| `--vep_cache_cpsr` | CPSR | ~24 GB download (release 115) |
| `--pypgx_bundle` | PyPGx star allele calling | ~370 MB |
| `--annotsv_annotations` | AnnotSV SV classification | ~5.3 GB download |
| `--cadd_snv`, `--spliceai_snv`, etc. | vcfanno score annotation | ~100 GB total |
| `--gnomad_constraint` | Slivar gene constraint | ~95 MB |
| `--pgs_scoring` | Polygenic risk scores | ~400 MB (`setup.sh` fills `<genome_dir>/prs_scores`) |
| `--ancestry_ref` | PRS percentiles and ancestry (step 26) | ~7 GB (`setup.sh --ancestry-panel`) |

Tools that require external databases (VEP, slivar, clinvar, CPSR, ExpansionHunter, HLA typing, pypgx, AnnotSV, CNVpytor) will **fail at startup** if enabled in `--tools` without their required parameters. vcfanno and prs are in the default tools but are skipped, with a warning, until a score file or `--pgs_scoring` is set. The gnomAD constraint table is optional for slivar; when it is set and no gene matches it, the task fails.

### Ancestry reference panel

`--ancestry_ref` is pgsc_calc's panel, `<genome_dir>/reference/pgsc_calc/pgsc_1000G_v1.tar.zst`, with the `pgsc_1000G_v1_GRCh38_sites.tsv` that `setup.sh --ancestry-panel` writes beside it (the run stops when the list is missing). With it, `prs` genotypes the panel's SNVs from the gVCF too, pgsc_calc projects the sample onto the panel, each score gets a percentile among the most similar population, and `PRS_SUMMARY` publishes step 26's table in `ancestry/`. `ancestry` in `--tools` needs `prs` and `--ancestry_ref`; without them it is skipped with a warning. See [step 26](26-ancestry.md).

### SV consensus merge

`survivor_merge` runs `SURVIVOR merge` (step 22's parameters: breakpoints within 1 kb, same type and strands, at least 50 bp, two or more callers) over the callers this run selected among `manta`, `delly` and `cnvpytor`, each cut to its PASS records first (`SURVIVOR_PREP`, `SURVIVOR_MERGE`, `SURVIVOR_SORT`). CNVpytor's depth-only calls count as one caller like the others; their coarse breakpoints often lie more than 1 kb from the paired-end callers'. TIDDIT and Sniffles2 have no module, so the consensus of the script (step 22) can hold those two callers as well. GRIDSS is in neither: its breakend (BND) records never match the DEL, DUP and INV records of the other callers. See [step 22](22-survivor-merge.md).

### FILTER=PASS required

ClinVar screen, clinical filter and slivar keep only records with FILTER=PASS. Before any analysis, `VCF_PRECHECK` counts the FILTER values of each sample. A VCF with no PASS record at all (for example unfiltered GATK HaplotypeCaller or FreeBayes output, where FILTER is `.`) stops the run with a message naming the sample, because every PASS-only step would report zero hits. Filter it with your caller's recommended filters, or add `--allow_unfiltered` to treat FILTER `.` as PASS for that file. A VCF with any PASS record is used as given.

### One sample per VCF

One samplesheet row is one sample, so a VCF given on a row must hold one sample column. `VCF_PRECHECK` stops the run, before any analysis, on a VCF with more than one (a joint-called family file, for example): the message names the sample count and the first five sample names, and prints the `bcftools view -s <name> -a -c 1` command that keeps one sample's column. That command writes a variant-only VCF, which is what the `vcf` column takes. A `gvcf` for the same row is split with `-s <name>` alone, so its reference blocks stay. `validate-setup.sh <sample>` applies the same rule to `vcf/<sample>.vcf.gz`, and fails when it cannot read that file's header.

### Contig names and gVCF input

`VCF_PRECHECK` also stops the run, before any analysis, in two cases, and the message names the sample and the fix:

- No contig that holds records is named the chr way (`1`, `MT` instead of `chr1`, `chrM`). Without this stop the mito haplogroup file comes out empty and chrX segments leak into the autosomal ROH summary, with exit 0. The message prints the `bcftools annotate --rename-chrs` command with its 25-line map.
- `pharmcat` is selected and the VCF is a gVCF (a `##GVCFBlock` header line, or reference-block records: ALT `<*>`, `<NON_REF>` or `.` with `INFO/END`), or only its name says so (`.g.vcf`, `.genomic.vcf`). PharmCAT refuses both. Without `pharmcat`, a gVCF runs: ROH, the mito haplogroup and the ClinVar screen give the same results as on the matching variants-only file.

[Starting from a Vendor VCF](vcf-first.md) has the commands that fix both.

### Security model

This pipeline is designed for **personal, single-user use** on trusted data. Sample labels are restricted to `[A-Za-z0-9._-]` (they name folders and go into shell commands); this is not anonymisation, see [Before you share outputs](#before-you-share-outputs). HTML report fields from VCF INFO are escaped to prevent XSS. However, it is **not hardened for multi-tenant or untrusted-input scenarios**. Do not expose the pipeline or its outputs as a web service without additional security review.

### Failed commands fail the task

Every task script runs under `bash -euo pipefail` (`process.shell` in `conf/base.config`). In a pipe such as `bcftools view -f PASS in.vcf.gz | bcftools +split-vep ...`, a failure of the first command, a truncated input for example, fails the task. Under Nextflow's default `bash -ue` only the last command counted, and the task wrote a short, valid file and passed.

### No network inside the containers

With `-profile docker` every container runs with `--network none` (`process.containerOptions` in `nextflow.config`). The steps read only their inputs, so a tool that tries to download something at run time fails instead of fetching an unpinned file. The `singularity` profile does not cut the network.

One process runs on the host instead of in a container: `PRS` starts pgsc_calc, a Nextflow pipeline of its own that starts its own containers. It gets the images of `versions.env` through `conf/containers.config` (`ext.pipeline_images`, written by `scripts/ci/gen-containers-config.sh`), and gives every one of its containers `--network none`. With `--pgsc_calc <genome_dir>/tools/pgsc_calc-<release>` (setup.sh installs it with its nf-schema plugin) it runs offline; without, the task first fetches GitHub's archive of the pinned release, checked against `PGSC_CALC_SHA256`, and Nextflow the plugin.

### Cyrius, opt-in

Cyrius (CYP2D6 star alleles) is under a non-commercial licence and is not in the default tools. `scripts/setup.sh --cyrius <genome_dir>` installs it once into `tools/cyrius-1.1.1/`, with `pip --require-hashes --no-deps --only-binary :all:` from `scripts/cyrius-constraints.txt`, which pins every wheel by its sha256. The `CYRIUS` module runs it from there (`--tools ...,cyrius --cyrius_install <genome_dir>/tools/cyrius-1.1.1`) with no network, like every other task. See [step 21](21-cyrius.md).

### CI validation scope

`nextflow.yml` lints the pipeline and stub-runs every process, on the pinned Nextflow release and the newest 26.04.x, with one samplesheet row per starting point; it also checks that `T1K_BUILD` ran once for all samples, and that a VCF+BAM-only samplesheet against a reference with no `.sr.mmi` beside it finishes without building an index or aligning anything. The E2E workflow runs the pipeline with real containers on a slice of the public HG002 genome: from the FASTQ pair through alignment, the sex check and DeepVariant to the default tools but cyrius (the fixture's BAM lacks the autosomal bins Cyrius normalises on), plus clinvar and hla_typing; through the VCF+BAM entry, over several cases, with clinvar, mosdepth, delly, manta, vcfanno, roh, pharmcat, cpic, pypgx, telomere_hunter, mito_haplogroup, survivor_merge, y_haplogroup, sample_qc, prs, cram_archive and html_report (the prs case's row names a gVCF too); and the sex-check stop (see [Testing](testing.md#the-e2e-job)). It then runs the parity check above. It cannot run the tools whose databases do not fit a CI runner (vep's offline cache, cpsr) or the ones it has no data for. Before trusting results from a new installation, run the pipeline on a known sample and compare key outputs (PharmCAT star alleles, ClinVar hit counts) against expected values.

---

## Relationship to nf-core

This pipeline uses [nf-core](https://nf-co.re/) template patterns and tooling for code quality, but is **not an official nf-core pipeline** (it uses a GPL-3.0 license; nf-core requires MIT).

Individual modules (PharmCAT, pypgx, slivar) will be contributed to [nf-core/modules](https://github.com/nf-core/modules) under MIT license for use by the broader community.

### Acknowledgement

> This pipeline was created using tools and best practices from the nf-core community (Ewels et al., 2020, Nat Biotechnol). nf-core components used here are released under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
