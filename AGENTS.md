# AGENTS.md — Genomics Pipeline

Instructions for AI agents working on this repository.

## Project Context

This is a public, open-source WGS (Whole Genome Sequencing) analysis pipeline designed for consumer hardware. It must be:
- **Generic**: No personal data, no hardcoded paths, no user-specific defaults
- **Reproducible**: Every command must work on any Linux amd64 machine with Docker
- **Well-documented**: Target audience includes non-bioinformaticians analyzing their own genome data

## Critical Rules

### No Personal Information
- NEVER commit personal paths (e.g., `/mnt/user/Multimedia/`, server hostnames, IP addresses)
- NEVER use specific sample names as defaults (use `your_name` or `$SAMPLE` placeholder)
- All environment variables must require user to set them: `${VAR:?Set VAR to...}`
- Docker mount point is always `:/genome` (not locale-specific)

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
- README.md step table must stay in sync with actual docs and scripts
- All Docker images must include the exact tag (not just `:latest` unless no versioned tags exist)

### Lessons Learned
- **ALWAYS update `docs/lessons-learned.md`** when encountering a new failure, workaround, or non-obvious behavior
- Include: what failed, why it failed, and the fix
- This is the most valuable document for future users — every Docker image issue, permission error, path confusion, and tool quirk should be recorded here

### Testing Changes
- After modifying any script, verify:
  1. No personal paths remain (`grep -r '/mnt/user\|watchtower\|sergio\|annais' scripts/ docs/`)
  2. All scripts use `GENOME_DIR` not `GENOMA_DIR`
  3. Docker mount is `:/genome` not `:/genoma`
  4. `shellcheck` passes on all scripts (if available)

### Git Practices
- Commit messages: descriptive, multi-line for large changes
- Never force-push to main
- Keep commits atomic: docs + scripts for the same feature in one commit

## Architecture

```
genomics-pipeline/
  README.md                    # Main entry point, pipeline overview, quick start
  LICENSE                      # GPL-3.0
  AGENTS.md                    # This file
  .gitignore                   # Excludes BAM, VCF, tar.gz, etc.
  docs/
    00-reference-setup.md      # One-time reference data downloads
    01-ora-to-fastq.md         # Step docs (one per pipeline step)
    ...
    20-mtoolbox.md
    hardware-requirements.md   # Disk, RAM, CPU, runtime breakdown
    vendor-guide.md            # Data formats from each WGS vendor
    interpreting-results.md    # Plain-language guide for non-experts
    multi-sample.md            # Comparing two or more samples (partners, family)
    glossary.md                # Alphabetical glossary of genomics terms
    quick-test.md              # Verify setup with public test data
    resources.md               # Free courses, databases, and learning resources
    troubleshooting.md         # Comprehensive troubleshooting by symptom
    lessons-learned.md         # Every failure and fix (KEEP UPDATED)
  scripts/
    01-ora-to-fastq.sh         # Step scripts (one per pipeline step)
    ...
    27-cpic-lookup.sh
    run-all.sh                 # Orchestrator: runs all steps with parallelism
    validate-setup.sh          # Pre-flight check: Docker, refs, images, sample
    generate-report.sh         # Text summary report aggregating all outputs
  .github/workflows/
    lint.yml                   # ShellCheck + markdownlint
    smoke-test.yml             # Dry-run validation of all scripts
```

## Data Flow

```
User's FASTQ/BAM/VCF
  │
  ├─ Step 2: minimap2 alignment (FASTQ → BAM)
  ├─ Step 3: DeepVariant variant calling (BAM → VCF)
  │
  ├─ VCF-dependent steps: 6, 7, 9, 11, 12, 13, 14, 17, 25, 26
  ├─ BAM-dependent steps: 4, 10, 15, 16, 18, 19, 20, 21
  ├─ Post-VCF-analysis: 22 (SV merge), 23 (clinical filter), 24 (report), 27 (CPIC)
  └─ Both: 5 (needs Manta VCF from step 4)
```

## When Adding a New Step

1. Create `docs/NN-tool-name.md` following the template of existing docs
2. Create `scripts/NN-tool-name.sh` following the script conventions above
3. Update `README.md` step table with the new step
4. Update `scripts/run-all.sh` to include the new step in the appropriate phase
5. Update `docs/00-reference-setup.md` if new reference data or Docker images are needed
6. Update `docs/interpreting-results.md` if the output needs explanation
7. Add the Docker image to the pre-pull list in `docs/00-reference-setup.md`
8. Test on at least one sample before committing

## Tool-Specific Gotchas (Learned from Real Execution)

### PharmCAT 2.15.5
- **Two-step workflow**: Preprocessor (`pharmcat_vcf_preprocessor.py` with `-refFna`) → main jar (`pharmcat.jar`). The old `-refFasta` flag on the jar no longer exists.
- Preprocessor outputs `.preprocessed.vcf.bgz` (NOT `.vcf`).
- JSON output structure: `genes` is `{source → {gene_name → data}}` (dict of dicts), NOT a list. `sourceDiplotypes` contains `allele1`/`allele2` objects with `.name` field.
- Star allele calls may differ from other pipelines (e.g., Sanitas hg19 vs our hg38 DeepVariant). PharmCAT 2.15.5 definitions update frequently.
- **Pipeline pin vs upstream**: The pipeline is currently pinned to PharmCAT `2.15.5` for reproducibility, but upstream PharmCAT releases continue to ship new guideline content and parser-relevant format changes. Before bumping the Docker tag, revalidate both step 7 and step 27 end-to-end — the JSON structure and preprocessor flags have changed between major versions.

### plink2 (PRS / Ancestry)
- **chrX requires sex info**: Use `--chr 1-22 --allow-extra-chr` for PRS/PCA (autosomal only).
- **`--output-chr chrM`** preserves `chr` prefix in output. Without it, `--chr 1-22` strips prefix → variant IDs become `1:pos` instead of `chr1:pos`.
- **`--set-all-var-ids '@:#'`**: The `@` placeholder includes the full contig name (including `chr`). Do NOT use `chr@:#` or you get `chrchr1:pos`.
- **Scoring file duplicates**: Large PGS Catalog files (e.g., PGS000014 with 7M variants) contain duplicate variant:allele pairs. Deduplicate before `--score` or plink2 errors.
- **LD pruning requires >=50 samples**. PCA requires >=2. Single-sample ancestry is fundamentally limited.
- **PRS guardrail**: Raw PRS scores are NOT percentiles, absolute risks, or portable labels across tool versions. Never describe them that way unless you have an ancestry-matched reference cohort scored with the exact same PGS file and preprocessing.
- **Ancestry guardrail**: Treat the current single-sample ancestry step as overlap/QC plus a starting point for downstream projection work, not as a population-placement tool by itself.

### bcftools
- **`bcftools sort` requires `##contig` headers** — fails silently or errors on VCFs without them. Always inject contig headers from the reference `.fai` when building VCFs.
- **`set -euo pipefail` + `find | grep -q`**: If the directory doesn't exist, `find` exits 1, which poisons pipefail even with `2>/dev/null`. Use per-directory flag variables instead.

### VEP
- Running without `--af_gnomade` produces VCF lacking gnomAD frequencies. The clinical filter (step 23) then can't filter by population frequency, resulting in thousands of unfiltered MODERATE variants.

### Cyrius (CYP2D6)
- Returns `None/None` for both samples — common limitation of short-read WGS due to CYP2D7 homology and structural rearrangements.

## Knowledge Base / Tool Update Cadence

| Resource / Tool | Update Frequency | Re-run Steps | Time |
|---|---|---|---|
| ClinVar | Monthly | 6 (ClinVar screen) | ~5 min |
| Ensembl / VEP cache | Each Ensembl release (~4-6 months) | 13, 23 | ~3 hr |
| PCGR/CPSR data | Annually or when upstream bundle changes materially | 17 | ~45 min |
| PharmCAT upstream release | Check quarterly | 7, 27 | ~15-30 min validation |
| CPIC static lookup table | Check quarterly or when CPIC adds/updates guideline pairs | 27 | ~15 min code refresh |
| PGS Catalog | Check quarterly | 25 | ~30 min |

ClinVar is the highest-value update — new pathogenic classifications happen monthly.
Before bumping PharmCAT, validate the preprocessor flags, JSON parsing in step 27, and any phenotype/diplotype changes on a known test sample.
For a public pipeline, keep PGS IDs, PharmCAT Docker tags, and the CPIC lookup table explicitly versioned in git so result changes are auditable over time.

## Common Issues When Developing

- **Docker image not found**: Biocontainer tags change frequently. Use `docker search` or check quay.io/biocontainers directly.
- **Permission denied in container**: Add `--user root`. Most bioinformatics images run as non-root.
- **0-byte output**: Usually means the input path was wrong inside the container. Double-check the `:/genome` mount mapping.
- **PCGR/CPSR path confusion**: `--pcgr_dir` should point to the PARENT of `data/`, not `data/` itself. CPSR appends `/data` internally.
- **VEP cache download**: Use `wget -c` (resume-capable), not VEP's `INSTALL.pl` which can't resume 26 GB downloads.

## Audience Reminder

Every decision should be evaluated through the lens of: "Would a non-bioinformatician who just received their WGS data be able to follow this?" If the answer is no, add more documentation, clearer error messages, or a simpler default.
