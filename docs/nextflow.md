# Nextflow Execution (v0.5.0)

The pipeline has a [Nextflow](https://www.nextflow.io/) DSL2 execution path for **post-calling interpretation and clinical analysis**. It accepts VCF + BAM from any upstream caller (e.g. nf-core/sarek, DRAGEN, the bash alignment scripts) and runs pharmacogenomics, variant annotation, clinical screening, structural variant analysis, and reporting across 6 workflows with 27 modules.

> **Both execution paths are maintained.** The bash scripts (`run-all.sh`) remain the simpler option for single-machine use. Nextflow adds automatic parallelism, content-hash resume, and HPC/Singularity support. Both paths produce biologically equivalent results, though output file names and report scope may differ.

---

## Quick Start

### Prerequisites

1. **Docker** (already required for the bash pipeline)
2. **Java 11-21** (Nextflow runtime requirement)
3. **Nextflow** — install with:
   ```bash
   curl -s https://get.nextflow.io | bash
   sudo mv nextflow /usr/local/bin/
   ```

### Run the Pipeline

```bash
# 1. Create a samplesheet CSV
cat > samplesheet.csv << 'EOF'
sample,vcf,vcf_index,bam,bam_index
sergio,/path/to/sergio.vcf.gz,/path/to/sergio.vcf.gz.tbi,/path/to/sergio_sorted.bam,/path/to/sergio_sorted.bam.bai
EOF

# 2. Run (default tools need no external databases)
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /path/to/Homo_sapiens_assembly38.fasta \
    --outdir ./results \
    -profile docker

# 3. To enable database-requiring tools, add them to --tools with their flags:
#    --tools '...,vep,slivar,clinical_filter'  + --vep_cache /path/to/vep_cache
#    --tools '...,cpsr'                        + --pcgr_data + --vep_cache_cpsr
#    --tools '...,clinvar'                     + --clinvar + --clinvar_index
#    --tools '...,expansion_hunter'            + --expansion_catalog
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

\* BAM is technically optional (VCF-only runs are valid for annotation and PGx), but most default tools (mosdepth, telomere_hunter, cyrius, mito_variants) and opt-in tools (expansion_hunter, hla_typing, pypgx) require BAM input. **Provide BAM for full analysis.**

### Using Sarek Output

If you ran [nf-core/sarek](https://nf-co.re/sarek) for alignment and variant calling, point the samplesheet at sarek's output files:

```csv
sample,vcf,vcf_index,bam,bam_index
sergio,results/variant_calling/deepvariant/sergio/sergio.deepvariant.vcf.gz,results/variant_calling/deepvariant/sergio/sergio.deepvariant.vcf.gz.tbi,results/preprocessing/recalibrated/sergio/sergio.recal.bam,results/preprocessing/recalibrated/sergio/sergio.recal.bam.bai
```

---

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (default for local) |
| `singularity` | Run with Singularity/Apptainer (HPC clusters) |
| `test` | Minimal test with reduced resources |
| `test_full` | Full-size test with real WGS data |

Combine profiles: `-profile docker,test`

---

## Resource Configuration

Default resource limits (tuned for 16-core consumer desktop):

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
├── sergio/
│   ├── pharmcat/           # PharmCAT PGx reports (HTML + JSON)
│   ├── clinvar/            # ClinVar pathogenic variant screen
│   ├── pypgx/              # pypgx star allele calling (optional)
│   ├── cpic/               # CPIC drug-gene recommendations (optional)
│   ├── vep/                # VEP + vcfanno enriched VCF (CADD, SpliceAI, REVEL, AlphaMissense)
│   ├── slivar/             # Prioritized variants + compound hets
│   ├── clinical/           # Clinically relevant variant subset
│   ├── cpsr/               # Cancer predisposition report
│   ├── roh/                # Runs of homozygosity
│   ├── prs/                # Polygenic risk scores
│   ├── ancestry/           # Ancestry PCA (optional)
│   ├── mito/               # Mitochondrial haplogroup
│   ├── hla/                # HLA typing
│   ├── expansion_hunter/   # Repeat expansion calls
│   ├── telomere/           # Telomere length estimation
│   ├── coverage/           # Coverage statistics (mosdepth)
│   ├── mito_variants/      # Mitochondrial variant calling
│   ├── cyrius/             # CYP2D6 star allele (Cyrius)
│   ├── manta/              # SV calling (optional)
│   ├── delly/              # SV calling (optional)
│   ├── cnvnator/           # CNV calling (optional)
│   └── *_report.html       # Consolidated HTML report (published to sample root)
└── pipeline_info/
    ├── timeline_*.html
    ├── report_*.html
    ├── trace_*.txt
    └── dag_*.svg
```

---

## Nextflow vs Bash: Which Should I Use?

| Feature | Bash (`run-all.sh`) | Nextflow (`main.nf`) |
|---------|---------------------|----------------------|
| Setup complexity | Just Docker | Docker + Java + Nextflow |
| Resume on failure | File-existence checks | Content-hash caching (more robust) |
| Parallelism | Manual (`wait`, throttle) | Automatic DAG-based |
| HPC / Singularity | Not supported | Built-in |
| Learning curve | Shell scripting | Nextflow DSL2 + Groovy |
| Target audience | Non-bioinformaticians | Bioinformaticians, HPC users |

**Recommendation:** If you're comfortable with bash and running on a single machine, use the bash scripts. If you want automatic parallelism, robust resume, or HPC support, use Nextflow.

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
| `--vep_cache` | VEP annotation | ~15 GB |
| `--pcgr_data` | CPSR cancer predisposition | ~20 GB |
| `--pypgx_bundle` | PyPGx star allele calling | ~2 GB |
| `--cadd_snv`, `--spliceai_snv`, etc. | vcfanno score annotation | ~100 GB total |
| `--gnomad_constraint` | Slivar gene constraint | ~5 MB |
| `--pgs_scoring` | Polygenic risk scores | varies |

Tools that require external databases (VEP, slivar, clinvar, CPSR, ExpansionHunter, HLA typing, pypgx) will **fail at startup** if enabled in `--tools` without their required parameters. Annotation scores (CADD, SpliceAI, REVEL, AlphaMissense, gnomAD constraint) are optional enrichments — vcfanno and slivar degrade gracefully without them.

### Ancestry reference panel

The `--ancestry_ref` parameter expects a **single VCF file** (not a directory). Single-sample PCA without a multi-population reference panel produces mathematically limited results — the module will run but report `pca_status: skipped_single_sample`. For meaningful ancestry estimation, provide a reference panel VCF containing multiple population samples.

### SV consensus merge (experimental)

The `survivor_merge` module uses a simplified bcftools-based heuristic (1kb position binning) rather than the full SURVIVOR or Jasmine algorithm. CNVnator calls (depth-based, no PASS/FAIL marking) are treated equally with paired-end callers in the "2+ callers" consensus. For production SV analysis, consider running SURVIVOR or Jasmine externally.

### Security model

This pipeline is designed for **personal, single-user use** on trusted data. Sample names are sanitized (alphanumeric, `.`, `_`, `-` only), and HTML report fields from VCF INFO are escaped to prevent XSS. However, it is **not hardened for multi-tenant or untrusted-input scenarios**. Do not expose the pipeline or its outputs as a web service without additional security review.

### Cyrius runtime installation

The Cyrius module (CYP2D6 star allele calling) installs `cyrius==1.1.1` via pip at runtime because no pre-built container image exists. This requires **network access on first run** and means Nextflow's container-only reproducibility guarantee does not fully apply to this module. The version is pinned to avoid floating dependencies. The matching bash script (`scripts/21-cyrius.sh`) has the same limitation.

### CI validation scope

The CI test suite validates the stub-testable subset of modules using `-stub` dry runs (tools that do not require external databases). It does **not** cover database-dependent tools (vep, cpsr, clinvar, expansion_hunter) or run real bioinformatics tools on real data. Before trusting results from a new installation, run the pipeline on a known sample and compare key outputs (PharmCAT star alleles, ClinVar hit counts, PCA eigenvectors) against expected values.

---

## Relationship to nf-core

This pipeline uses [nf-core](https://nf-co.re/) template patterns and tooling for code quality, but is **not an official nf-core pipeline** (it uses a GPL-3.0 license; nf-core requires MIT).

Individual modules (PharmCAT, pypgx, slivar) will be contributed to [nf-core/modules](https://github.com/nf-core/modules) under MIT license for use by the broader community.

### Acknowledgement

> This pipeline was created using tools and best practices from the nf-core community (Ewels et al., 2020, Nat Biotechnol). nf-core components used here are released under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
