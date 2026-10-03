# Nextflow Execution

The pipeline has a [Nextflow](https://www.nextflow.io/) DSL2 execution path for **post-calling interpretation and clinical analysis**. It accepts VCF + BAM from any upstream caller (e.g. nf-core/sarek, DRAGEN, the bash alignment scripts) and runs pharmacogenomics, variant annotation, clinical screening, structural variant analysis, and reporting across 6 workflows: 35 processes in 29 module files under `modules/local/`. The VCF needs FILTER=PASS records and GRCh38 contig names with chr; see [FILTER=PASS required](#filterpass-required) and [Contig names and gVCF input](#contig-names-and-gvcf-input). Starting from a provider's VCF: [Starting from a Vendor VCF](vcf-first.md).

> **Both execution paths are maintained.** The bash scripts (`run-all.sh`) remain the simpler option for single-machine use. Nextflow adds automatic parallelism and content-hash resume. It has a Singularity profile, but that profile is untested (see [Profiles](#profiles)). Both paths produce biologically equivalent results, though output file names and report scope may differ.

---

## Quick Start

### Prerequisites

1. **Docker** (already required for the bash pipeline)
2. **Java 17 or later** (Nextflow 25.10 runtime requirement; CI runs Java 17)
3. **Nextflow 25.10.4**, the version CI validates. Pin it when installing, because the plain installer fetches the newest release:
   ```bash
   curl -s https://get.nextflow.io | NXF_VER=25.10.4 bash
   sudo mv nextflow /usr/local/bin/
   ```

### Run the Pipeline

```bash
# 1. Create a samplesheet CSV
cat > samplesheet.csv << 'EOF'
sample,vcf,vcf_index,bam,bam_index
sample1,/path/to/sample1.vcf.gz,/path/to/sample1.vcf.gz.tbi,/path/to/sample1_sorted.bam,/path/to/sample1_sorted.bam.bai
EOF

# 2. Run (default tools need no external databases; prs and vcfanno are
#    skipped with a warning until --pgs_scoring or a score file is set)
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /path/to/Homo_sapiens_assembly38.fasta \
    --outdir ./results \
    -profile docker

# 3. To enable database-requiring tools, add them to --tools with their flags:
#    --tools '...,vep,slivar,clinical_filter'  + --vep_cache /path/to/vep_cache
#    --tools '...,cpsr'                        + --pcgr_data + --vep_cache_cpsr
#    --tools '...,clinvar'                     + --clinvar + --clinvar_index
#    --tools '...,expansion_hunter'            + --expansion_catalog (and a sex column)
#    --tools '...,annotsv'                     + --annotsv_annotations
#    --tools '...,cnvpytor'                    + --cnvpytor_resources
#    --tools '...,delly'                       (optional --delly_exclude <excl.tsv>, passed as delly call -x)
#    An unknown name in --tools stops the run.
```

### Resume After Failure

Nextflow caches completed steps using content hashes. If a step fails, fix the issue and resume:

```bash
nextflow run main.nf -resume [same params as before]
```

Only the failed and downstream steps re-run.

---

## Samplesheet Format

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier (used as output directory name) |
| `vcf` | Yes | Path to bgzipped VCF (`.vcf.gz`) |
| `vcf_index` | Yes | Path to tabix index (`.vcf.gz.tbi`) |
| `bam` | No* | Path to aligned BAM (needed for BAM-based steps like pypgx) |
| `bam_index` | No* | Path to BAM index (`.bam.bai`) |
| `sex` | No** | `male` or `female` |

\* BAM is technically optional (VCF-only runs are valid for annotation and PGx), but most default tools (mosdepth, telomere_hunter, cyrius, mito_variants) and opt-in tools (expansion_hunter, hla_typing, pypgx) require BAM input. **Provide BAM for full analysis.**

The VCF must name its contigs the GRCh38 way with chr (`chr1` to `chr22`, `chrX`, `chrY`, `chrM`); a VCF named `1`, `MT` stops the run with the rename command. A gVCF is the better PharmCAT input, but the pipeline does not expand its reference blocks yet, so with `pharmcat` selected a gVCF stops the run; a variants-only VCF leaves about half of PharmCAT's genes Unknown. [Starting from a Vendor VCF](vcf-first.md) has the commands for both.

\*\* `sex` is required on every row that has a BAM when `expansion_hunter` is in `--tools`: it sets the chrX ploidy, and ExpansionHunter's default is female. A BAM row without it stops the run at parse time. VCF-only rows never reach ExpansionHunter, so they need no `sex`.

Each `sample` value must appear once; a repeated id stops the run, because the id names the output directory and keys every per-sample join.

### Using Sarek Output

If you ran [nf-core/sarek](https://nf-co.re/sarek) for alignment and variant calling, point the samplesheet at sarek's output files. Sarek 3.x writes CRAM by default, and this pipeline reads BAM only, so run sarek with `--save_output_as_bam` to get the `.recal.bam` files below (or convert the CRAM with `samtools view -b -T <reference>`):

```csv
sample,vcf,vcf_index,bam,bam_index
sample1,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz.tbi,results/preprocessing/recalibrated/sample1/sample1.recal.bam,results/preprocessing/recalibrated/sample1/sample1.recal.bam.bai
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (default for local) |
| `singularity` | Singularity/Apptainer. **Untested.** Three modules write inside their image and need a writable container: `cyrius` (pip-installs at run time), `pypgx` (links its bundle under `/root`) and `cnvpytor` (copies resources into its `site-packages`). |
| `test` | Minimal test with reduced resources |
| `test_full` | Full-size test with real WGS data |

Combine profiles: `-profile docker,test`

---

## Resource Configuration

Default resource limits (tuned for 16-core consumer desktop). Each process label asks for a fixed CPU count and a memory and time that double on the one retry after an out-of-memory or time-limit exit; these limits cap every request:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_cpus` | 16 | Maximum CPUs per process |
| `--max_memory` | 64.GB | Maximum memory per process |
| `--max_time` | 48.h | Maximum wall time per process |

Override for smaller machines:

```bash
nextflow run main.nf --max_cpus 8 --max_memory 32.GB [other params]
```

---

## Output Structure

```
results/
├── sample1/
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
│   ├── ancestry/           # Ancestry PCA (optional)
│   ├── mito/               # Mitochondrial haplogroup and mitochondrial variant calls
│   ├── hla/                # HLA typing
│   ├── expansion_hunter/   # Repeat expansion calls
│   ├── telomere/           # Telomere length estimation
│   ├── coverage/           # Coverage statistics (mosdepth)
│   ├── cyrius/             # CYP2D6 star allele (Cyrius)
│   ├── manta/              # SV calling (optional)
│   ├── sv_duphold/         # Manta SVs with duphold depth tags (optional)
│   ├── sv_filtered/        # Manta SVs after the duphold depth filter (optional)
│   ├── annotsv/            # AnnotSV ACMG classification of the filtered SVs (optional)
│   ├── delly/              # SV calling (optional)
│   ├── cnvpytor/           # CNV calling (optional)
│   ├── sv_merged/          # SV consensus of two or more callers (optional)
│   └── *_report.html       # Summary HTML report (published to sample root): ClinVar, PharmCAT, CPIC,
│                           #   CPSR, clinical filter, slivar, ROH, mito haplogroup; "Not run" for a tool not selected
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
- **The input file name** is in the `##bcftools_viewCommand` and `##bcftools_normCommand` lines of the ClinVar VCFs, and can be in the command lines other tools print into their outputs.
- **`pipeline_info/`** (report, timeline, trace) holds absolute paths of your machine.
- **The samplesheet's `sample` label** names every output folder and file.

The fix is a neutral input. Use a label that does not name you, name the file after it, and keep only the header lines the tools read. The recipe below lists the lines to keep and renames the sample column. A list of keys to delete would miss lines, because callers name them differently.

```bash
bcftools view --no-version -h in.vcf.gz | grep -E '^(##(fileformat|FILTER|INFO|FORMAT|ALT|contig)=|#CHROM)' > h.txt
echo SAMPLE > names.txt
bcftools reheader -h h.txt -s names.txt -o SAMPLE.vcf.gz in.vcf.gz
bcftools index -t SAMPLE.vcf.gz
```

`bcftools annotate -x` cannot remove these lines (it exits with "No matching tag"). The same recipe, with bcftools from the pinned image, is step 4 of [Starting from a Vendor VCF](vcf-first.md). Do not share `pipeline_info/`.

---

## Nextflow vs Bash: Which Should I Use?

| Feature | Bash (`run-all.sh`) | Nextflow (`main.nf`) |
|---------|---------------------|----------------------|
| Setup complexity | Just Docker | Docker + Java + Nextflow |
| Resume on failure | File-existence checks | Content-hash caching (more robust) |
| Parallelism | Manual (`wait`, throttle) | Automatic DAG-based |
| HPC / Singularity | Not supported | Profile exists, untested |
| Learning curve | Shell scripting | Nextflow DSL2 + Groovy |
| Target audience | Non-bioinformaticians | Bioinformaticians, HPC users |

**Recommendation:** If you're comfortable with bash and running on a single machine, use the bash scripts. If you want automatic parallelism or robust resume, use Nextflow.

---

## Known Limitations & Design Decisions

### Post-calling scope

This Nextflow pipeline is a **post-calling interpretation pipeline**, not a FASTQ-to-results pipeline. It accepts VCF + BAM from any upstream caller (e.g. nf-core/sarek, DRAGEN, the bash alignment scripts) and runs pharmacogenomics, annotation, clinical screening, structural variant calling, and reporting. Alignment and primary variant calling are handled upstream.

### Bash vs Nextflow parity

Both execution paths (bash `run-all.sh` and Nextflow `main.nf`) aim for **biologically equivalent results** — the same clinical conclusions, gene calls, and risk assessments. However, they are **not output-identical**: file names, directory structure, report formatting, and intermediate files may differ. When in doubt, the bash scripts are the reference implementation.

### Reference databases not auto-downloaded

Several tools require large reference databases that are **not automatically downloaded** by the pipeline. You must obtain and provide paths for these yourself:

| Parameter | Required by | Size |
|-----------|------------|------|
| `--vep_cache` | VEP annotation | ~28 GB download (release 116) |
| `--pcgr_data` | CPSR cancer predisposition | ~5 GB bundle (20250314) |
| `--vep_cache_cpsr` | CPSR | ~25 GB download (release 113) |
| `--pypgx_bundle` | PyPGx star allele calling | ~370 MB |
| `--annotsv_annotations` | AnnotSV SV classification | ~5.3 GB download |
| `--cadd_snv`, `--spliceai_snv`, etc. | vcfanno score annotation | ~100 GB total |
| `--gnomad_constraint` | Slivar gene constraint | ~95 MB |
| `--pgs_scoring` | Polygenic risk scores | varies |

Tools that require external databases (VEP, slivar, clinvar, CPSR, ExpansionHunter, HLA typing, pypgx, AnnotSV, CNVpytor) will **fail at startup** if enabled in `--tools` without their required parameters. vcfanno and prs are in the default tools but are skipped, with a warning, until a score file or `--pgs_scoring` is set. The gnomAD constraint table is optional for slivar; when it is set and no gene matches it, the task fails.

### Ancestry reference panel

The `--ancestry_ref` parameter expects a **single VCF file** (not a directory). The step reads the panel only for its list of variant ids: it never merges the panel's genotypes with the sample or projects the sample onto the panel. On one sample it therefore produces a SNP overlap count and `pca_status: skipped_single_sample`, whatever panel is given.

### SV consensus merge (experimental)

The `survivor_merge` module uses a simplified bcftools-based heuristic (1kb position binning) rather than the full SURVIVOR or Jasmine algorithm. CNVpytor calls (depth-based, no PASS/FAIL marking) are treated equally with paired-end callers in the "2+ callers" consensus. For production SV analysis, consider running SURVIVOR or Jasmine externally.

### FILTER=PASS required

ClinVar screen, clinical filter and slivar keep only records with FILTER=PASS. Before any analysis, `VCF_PRECHECK` counts the FILTER values of each sample. A VCF with no PASS record at all (for example unfiltered GATK HaplotypeCaller or FreeBayes output, where FILTER is `.`) stops the run with a message naming the sample, because every PASS-only step would report zero hits. Filter it with your caller's recommended filters, or add `--allow_unfiltered` to treat FILTER `.` as PASS for that file. A VCF with any PASS record is used as given.

### Contig names and gVCF input

`VCF_PRECHECK` also stops the run, before any analysis, in two cases, and the message names the sample and the fix:

- No contig that holds records is named the chr way (`1`, `MT` instead of `chr1`, `chrM`). Without this stop the mito haplogroup file comes out empty and chrX segments leak into the autosomal ROH summary, with exit 0. The message prints the `bcftools annotate --rename-chrs` command with its 25-line map.
- `pharmcat` is selected and the VCF is a gVCF (a `##GVCFBlock` header line, or reference-block records: ALT `<*>`, `<NON_REF>` or `.` with `INFO/END`), or only its name says so (`.g.vcf`, `.genomic.vcf`). PharmCAT refuses both. Without `pharmcat`, a gVCF runs: ROH, the mito haplogroup and the ClinVar screen give the same results as on the matching variants-only file.

[Starting from a Vendor VCF](vcf-first.md) has the commands that fix both.

### Security model

This pipeline is designed for **personal, single-user use** on trusted data. Sample labels are restricted to `[A-Za-z0-9._-]` (they name folders and go into shell commands); this is not anonymisation, see [Before you share outputs](#before-you-share-outputs). HTML report fields from VCF INFO are escaped to prevent XSS. However, it is **not hardened for multi-tenant or untrusted-input scenarios**. Do not expose the pipeline or its outputs as a web service without additional security review.

### Cyrius runtime installation

The Cyrius module (CYP2D6 star allele calling) installs `cyrius==1.1.1` via pip at runtime because no pre-built container image exists. This requires **network access on every run** and means Nextflow's container-only reproducibility guarantee does not apply to this module. Only Cyrius itself is pinned here: its dependencies (pysam, numpy, scipy, statsmodels) and the `PYTHON_IMAGE` base tag are not, so they resolve to whatever is newest on the day. The bash script (`scripts/21-cyrius.sh`) also pins `cyrius==1.1.1` and holds the dependencies to the versions in `scripts/cyrius-constraints.txt`; the base tag moves there too.

### CI validation scope

The CI test suite validates the stub-testable subset of modules using `-stub` dry runs (tools that do not require external databases). It does **not** cover database-dependent tools (vep, cpsr, clinvar, expansion_hunter) or run real bioinformatics tools on real data. Before trusting results from a new installation, run the pipeline on a known sample and compare key outputs (PharmCAT star alleles, ClinVar hit counts, PCA eigenvectors) against expected values.

---

## Relationship to nf-core

This pipeline uses [nf-core](https://nf-co.re/) template patterns and tooling for code quality, but is **not an official nf-core pipeline** (it uses a GPL-3.0 license; nf-core requires MIT).

Individual modules (PharmCAT, pypgx, slivar) will be contributed to [nf-core/modules](https://github.com/nf-core/modules) under MIT license for use by the broader community.

### Acknowledgement

> This pipeline was created using tools and best practices from the nf-core community (Ewels et al., 2020, Nat Biotechnol). nf-core components used here are released under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
