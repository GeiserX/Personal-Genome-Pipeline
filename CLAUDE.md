# CLAUDE.md — Personal Genome Pipeline

## Overview
Whole genome sequencing (WGS) analysis pipeline for consumer hardware. Takes raw FASTQ/BAM/VCF data and runs 34 analysis steps locally in Docker containers: variant calling, pharmacogenomics, structural variants, cancer predisposition, polygenic risk scores, ancestry, telomere length, and more. Designed for non-bioinformaticians analyzing their own genome data.

## Tech Stack
- Bash (pipeline scripts, `set -euo pipefail`)
- Docker (all bioinformatics tools containerized)
- Key tools: DeepVariant, minimap2, BWA-MEM2, VEP, PharmCAT, GATK, FreeBayes, Strelka2, TIDDIT, Manta, PCGR/CPSR, plink2

## Development

```bash
# Validate setup
./scripts/validate-setup.sh

# Run all steps
GENOME_DIR=/path/to/data ./scripts/run-all.sh sample_name

# Run individual step
GENOME_DIR=/path/to/data ./scripts/03-deepvariant.sh sample_name

# Lint
shellcheck scripts/*.sh
```

Requirements: 16+ cores recommended, 500 GB disk per sample, Docker. Runs on Linux, macOS, WSL2.

### Testing Changes

After modifying any script, verify:
1. No personal paths remain (`grep -r '/mnt/user\|internal-host\|sample1\|sample2' scripts/ docs/`)
2. All scripts use `GENOME_DIR` not `GENOMA_DIR`
3. Docker mount is `:/genome` not `:/genoma`
4. `shellcheck` passes on all scripts

## Architecture

```
personal-genome-pipeline/
  README.md                    # Pipeline overview, quick start
  docs/
    00-reference-setup.md      # One-time reference data downloads
    01-ora-to-fastq.md         # Step docs (one per pipeline step)
    ...
    hardware-requirements.md   # Disk, RAM, CPU, runtime breakdown
    vendor-guide.md            # Data formats from each WGS vendor
    chip-data-guide.md         # Using 23andMe/MyHeritage/AncestryDNA chip data
    interpreting-results.md    # Plain-language guide for non-experts
    multi-sample.md            # Comparing two or more samples
    glossary.md                # Genomics terms
    quick-test.md              # Verify setup with public test data
    troubleshooting.md         # Comprehensive troubleshooting
    lessons-learned.md         # Every failure and fix (KEEP UPDATED)
  scripts/
    01-ora-to-fastq.sh         # Step scripts (one per pipeline step)
    ...
    27-cpic-lookup.sh
    chip-to-vcf.sh             # Chip data converter
    02a-alignment-bwamem2.sh   # Alternative aligner
    03a-gatk-haplotypecaller.sh # Alternative caller
    03b-freebayes.sh           # Alternative caller
    03c-strelka2-germline.sh   # Alternative caller
    04a-tiddit.sh              # Alternative SV caller
    benchmark-variants.sh      # Concordance benchmarking
    run-all.sh                 # Orchestrator
    validate-setup.sh          # Pre-flight check
    generate-report.sh         # Summary report
  .github/workflows/
    lint.yml                   # ShellCheck + markdownlint
    smoke-test.yml             # Dry-run validation
```

### Data Flow

```
User's FASTQ/BAM/VCF
  ├─ Step 2: minimap2 alignment (FASTQ -> BAM)
  ├─ Step 3: DeepVariant variant calling (BAM -> VCF)
  ├─ VCF-dependent steps: 6, 7, 9, 11, 12, 13, 14, 17, 25, 26
  ├─ BAM-dependent steps: 4, 10, 15, 16, 18, 19, 20, 21
  ├─ Post-VCF-analysis: 22 (SV merge), 23 (clinical filter), 24 (report), 27 (CPIC)
  └─ Both: 5 (needs Manta VCF from step 4)
```

## Key Rules

### No Personal Information
- NEVER commit personal paths, server hostnames, or IP addresses
- NEVER use specific sample names as defaults (use `your_name` or `$SAMPLE` placeholder)
- All environment variables must require user to set them: `${VAR:?Set VAR to...}`
- Docker mount point is always `:/genome`

### Script Conventions
- Shebang: `#!/usr/bin/env bash`
- Error handling: `set -euo pipefail`
- Parameters: `SAMPLE=${1:?Usage: $0 <sample_name>}`
- Environment: `GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}`
- Docker: always use `--cpus N --memory Xg` limits, `-v "${GENOME_DIR}:/genome"` mount, `--rm` flag
- Add `--user root` when the container needs write access to bind mounts
- Validate all input files exist before running Docker commands
- Print clear status messages: step name, input files, output location

### Documentation Conventions
- Each pipeline step has a matching doc in `docs/XX-name.md` and script in `scripts/XX-name.sh`
- Docs must include: What it does, Why, Tool name, Docker image, Command, Output, Runtime estimate, Notes
- README step table must stay in sync with actual docs and scripts
- All Docker images must include the exact tag (not floating)

### Lessons Learned
- **ALWAYS update `docs/lessons-learned.md`** when encountering a new failure, workaround, or non-obvious behavior
- Include: what failed, why it failed, and the fix

### Adding a New Step

1. Create `docs/NN-tool-name.md` following existing template
2. Create `scripts/NN-tool-name.sh` following script conventions
3. Update `README.md` step table
4. Update `scripts/run-all.sh` with the new step
5. Update `docs/00-reference-setup.md` if new reference data needed
6. Update `docs/interpreting-results.md` if output needs explanation
7. Add Docker image to pre-pull list in `docs/00-reference-setup.md`
8. Test on at least one sample before committing

- All processing is local; genomic data never leaves the machine
- Pin tool versions; never use floating tags
- Reference genome: GRCh38/hg38
- License: GPL-3.0

## Tool-Specific Gotchas

### PharmCAT 2.15.5
- **Two-step workflow**: Preprocessor (`pharmcat_vcf_preprocessor.py` with `-refFna`) then main jar (`pharmcat.jar`). The old `-refFasta` flag no longer exists.
- Preprocessor outputs `.preprocessed.vcf.bgz` (NOT `.vcf`).
- JSON output: `genes` is `{source -> {gene_name -> data}}` (dict of dicts, NOT a list). `sourceDiplotypes` contains `allele1`/`allele2` objects with `.name` field.
- Pipeline pinned to 2.15.5. Before bumping, revalidate steps 7 and 27 end-to-end — JSON structure and preprocessor flags change between versions.

### plink2 (PRS / Ancestry)
- **chrX requires sex info**: Use `--chr 1-22 --allow-extra-chr` for PRS/PCA (autosomal only).
- **`--output-chr chrM`** preserves `chr` prefix. Without it, prefix is stripped.
- **`--set-all-var-ids '@:#'`**: `@` includes full contig name. Do NOT use `chr@:#`.
- **Scoring file duplicates**: Large PGS files contain duplicate variant:allele pairs. Deduplicate before `--score`.
- **LD pruning requires >=50 samples**. PCA requires >=2. Single-sample ancestry is fundamentally limited.
- **PRS guardrail**: Raw scores are NOT percentiles or portable labels. Require ancestry-matched reference cohort.
- **Ancestry guardrail**: Single-sample step is a starting point, not a population-placement tool.

### Chip Data Conversion
- **NEVER use plink 1.9 for single-sample chip-to-VCF.** plink's `.bim` format encodes monomorphic sites with one allele. For single-sample data, ALL homozygous positions are monomorphic. `--ref-from-fa` cannot fix these. Result: all hom-ALT genotypes silently become hom-REF.
- **Use `bcftools convert --tsv2vcf -f <reference.fa>`** instead.
- **MyHeritage CSV needs pre-conversion** to TSV format.
- **hg19 VCF needs chr prefix** before liftover — use `bcftools annotate --rename-chrs`.
- **PharmCAT on chip data**: Misses CYP2C19 (25 positions), VKORC1 (1 position), miscalls CYP3A5.
- **ROH on chip data** requires `-G30` flag (no FORMAT/PL tags).
- **PRS on chip data** requires `no-mean-imputation` flag. Matches ~12% of large scoring files vs ~28% from WGS.

### Alternative Callers & Benchmarking
- **Output isolation**: Alternative tools write to separate directories to never overwrite defaults.
- **INTERVALS env var**: GATK and FreeBayes support `INTERVALS=chr22`. Strelka2 and TIDDIT do not.
- **Strelka2 is a small-variant caller** (SNVs + indels <=49bp), not an SV caller. Scoring model trained on BWA-MEM data; SNP precision drops with minimap2.
- **FreeBayes is single-threaded**: Full WGS ~9 hours. Needs `--memory 32g`.
- **GATK full-genome**: ~8.6 hours on i5-14500. Requires `.dict` file alongside FASTA.
- **BWA-MEM2 index files**: Created alongside FASTA. Check for `.bwt.2bit.64`.
- **ALIGN_DIR env var**: All alternative scripts accept `ALIGN_DIR=aligned_bwamem2`.
- **TIDDIT --skip_assembly auto-detected**: Checks for BWA index files.
- **benchmark-variants.sh**: Pairwise mode (auto-discovers vcf*/ dirs) and truth set mode (hap.py).

### bcftools
- **`bcftools sort` requires `##contig` headers** — fails on VCFs without them. Inject from reference `.fai`.
- **`set -euo pipefail` + `find | grep -q`**: If directory doesn't exist, `find` exits 1, poisoning pipefail. Use per-directory flag variables.

### VEP
- Running without `--af_gnomade` produces VCF lacking gnomAD frequencies. Clinical filter (step 23) then can't filter by population frequency.

### Cyrius (CYP2D6)
- Returns `None/None` — common limitation of short-read WGS due to CYP2D7 homology.

## Knowledge Base / Tool Update Cadence

| Resource / Tool | Update Frequency | Re-run Steps | Time |
|---|---|---|---|
| ClinVar | Monthly (first Thursday) | 6 (ClinVar screen) | ~5 min |
| Ensembl / VEP cache | ~6 months | 13, 23 | ~3 hr |
| PCGR/CPSR data | Annually | 17 | ~45 min |
| PharmCAT | Quarterly check | 7, 27 | ~15-30 min |
| CPIC / ClinPGx | Quarterly check | 27 | ~15 min |
| PGS Catalog | Quarterly check | 25 | ~30 min |

### Minimal Revalidation Before Publishing Updates

1. **ClinVar / VEP**: Run steps 6 and 23, compare pathogenic hits and filtered variants against previous run.
2. **PharmCAT / CPIC**: Run steps 7 and 27, diff diplotypes, phenotypes, and recommendations.
3. **PGS Catalog**: Rerun step 25, compare `variants_used/variants_total` and raw score deltas. New scoring file version = new baseline.
4. **Documentation**: Update pinned versions and interpretation guardrails before merging.

## Common Issues

- **Docker image not found**: Biocontainer tags change frequently. Check quay.io/biocontainers directly.
- **Permission denied in container**: Add `--user root`. Most bioinformatics images run as non-root.
- **0-byte output**: Usually wrong input path inside container. Double-check `:/genome` mount mapping.
- **PCGR/CPSR path confusion**: `--pcgr_dir` should point to PARENT of `data/`, not `data/` itself.
- **VEP cache download**: Use `wget -c` (resume-capable), not VEP's `INSTALL.pl` (can't resume 26 GB).

*Generated by [LynxPrompt](https://lynxprompt.com) CLI*
