# Nextflow Execution (v0.5.0+)

Starting with v0.5.0, the pipeline can be run via [Nextflow](https://www.nextflow.io/) as an alternative to the bash scripts. Nextflow provides DAG-based parallelism, robust resume-on-failure, and container orchestration.

> **The bash scripts remain first-class.** If you prefer the simplicity of `./scripts/run-all.sh`, nothing has changed. Nextflow is an additional execution path for users who want workflow engine features.

---

## Quick Start

### Prerequisites

1. **Docker** (already required for the bash pipeline)
2. **Java 11+** (Nextflow runtime requirement)
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
sample1,/path/to/sample1.vcf.gz,/path/to/sample1.vcf.gz.tbi,/path/to/sample1_sorted.bam,/path/to/sample1_sorted.bam.bai
EOF

# 2. Run
nextflow run main.nf \
    --input samplesheet.csv \
    --reference /path/to/Homo_sapiens_assembly38.fasta \
    --clinvar /path/to/clinvar_pathogenic_chr.vcf.gz \
    --clinvar_index /path/to/clinvar_pathogenic_chr.vcf.gz.tbi \
    --outdir ./results \
    -profile docker
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
| `bam` | No | Path to aligned BAM (needed for BAM-based steps like pypgx) |
| `bam_index` | No | Path to BAM index (`.bam.bai`) |

### Using Sarek Output

If you ran [nf-core/sarek](https://nf-co.re/sarek) for alignment and variant calling, point the samplesheet at sarek's output files:

```csv
sample,vcf,vcf_index,bam,bam_index
sample1,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz,results/variant_calling/deepvariant/sample1/sample1.deepvariant.vcf.gz.tbi,results/preprocessing/recalibrated/sample1/sample1.recal.bam,results/preprocessing/recalibrated/sample1/sample1.recal.bam.bai
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
├── sample1/
│   ├── pharmcat/
│   │   ├── sample1.report.html
│   │   └── sample1.report.json
│   └── clinvar/
│       ├── isec/
│       │   └── 0002.vcf          # Shared pathogenic variants
│       └── sample1_pass.vcf.gz
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

## Relationship to nf-core

This pipeline uses [nf-core](https://nf-co.re/) template patterns and tooling for code quality, but is **not an official nf-core pipeline** (it uses a GPL-3.0 license; nf-core requires MIT).

Individual modules (PharmCAT, pypgx, slivar) are contributed to [nf-core/modules](https://github.com/nf-core/modules) under MIT license for use by the broader community.

### Acknowledgement

> This pipeline was created using tools and best practices from the nf-core community (Ewels et al., 2020, Nat Biotechnol). nf-core components used here are released under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
