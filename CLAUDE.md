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
1. No personal paths or identifiers remain — grep `scripts/` and `docs/` for your own home paths, server hostnames, and names (e.g. `grep -ri '/mnt/user\|/Users/\|/home/' scripts/ docs/`)
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
  scripts/lib/common.sh        # Sourced by every script: versions.env, run_in, fetch, validate_sample
  versions.env                 # Every image tag and coupled data version, one line each
  .github/workflows/
    lint.yml                   # ShellCheck, actionlint, gitleaks, doc links
    personal-data.yml          # Personal-data scan of every tracked text file
    guard.yml                  # image variables, image tags (check-images.sh), coupled versions, fake-docker suite, helper binaries, unit tests; ci-ok sums them up
    container-test.yml         # Runs each changed image on the fixture (tests/smoke/commands.tsv); container-ok sums it up
    e2e.yml                    # Real tools on a small fixture
    nextflow.yml               # Nextflow config, lint, schema and stub runs on the pinned release and the newest 26.04.x
    renovate-dry-run.yml       # Renovate lookups without an app, and their checks
    freshness.yml              # Monthly issue of pins that fell behind upstream
    docs.yml                   # mkdocs build
    release.yml, stale.yml     # Releases, stale issues
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
- Source the library right after reading `SAMPLE` and `GENOME_DIR`: `. "$(dirname "$0")/lib/common.sh"`, then `validate_sample "$SAMPLE"`
- Containers: `run_in --cpus N --memory Xg "${TOOL_IMAGE}" tool ...`. It mounts `GENOME_DIR` read-only at `/genome` and the sample directory writable, with no network, as the calling user
- Opt out at the call, with the reason in a comment: `--rw DIR` (shared index or database), `--net` (the step downloads), `--root` (the image cannot run unprivileged)
- Images come from `versions.env` as quoted variables; a script never spells an image name or tag
- Reference: `${REF_FASTA}` on the host, `${REF_FASTA_C}` inside a container; never spell the reference file name in a step script
- Downloads: `fetch URL DEST [md5|sha256|sum VALUE-or-URL]`
- Validate all input files exist before running Docker commands
- Print clear status messages: step name, input files, output location

### Documentation Conventions
- Each pipeline step has a matching doc in `docs/XX-name.md` and script in `scripts/XX-name.sh`
- Docs must include: What it does, Why, Tool name, Docker image, Command, Output, Runtime estimate, Notes
- The category table in `docs/pipeline-overview.md` and the `nav:` of `mkdocs.yml` must stay in sync with actual docs and scripts
- All Docker images must include the exact tag (not floating)

### Lessons Learned
- **ALWAYS update `docs/lessons-learned.md`** when encountering a new failure, workaround, or non-obvious behavior
- Include: what failed, why it failed, and the fix

### Adding a New Step

1. Add the image as one line in `versions.env` (`setup.sh` and `validate-setup.sh` read their list from it)
2. Create `scripts/NN-tool-name.sh` following script conventions
3. Create `modules/local/<tool>/main.nf` with no `container` line, add its process to the table in `scripts/ci/gen-containers-config.sh` and run it; or note in `docs/nextflow.md` why the step stays bash-only
4. Create `docs/NN-tool-name.md` following existing template, and add it to `nav:` in `mkdocs.yml`
5. Add the step to the category table in `docs/pipeline-overview.md`
6. Update `scripts/run-all.sh` with the new step
7. Add a row for the image to `tests/smoke/commands.tsv`: a real command on the fixture and a check on what it wrote (an image with no row fails `container-test.yml`)
8. Update `docs/00-reference-setup.md` if new reference data needed
9. Update `docs/interpreting-results.md` if output needs explanation
10. Test on at least one sample before committing

### Bumping a Tool

Change its line in `versions.env`, plus the coupled data variable its comment names. Run `scripts/ci/gen-containers-config.sh` and `scripts/ci/gen-versions-doc.sh` to rewrite `conf/containers.config` and `docs/versions.md` (CI fails until both match). `container-test.yml` then runs the new image on the fixture. No script or module names the tag. Lines marked `hold:` or `legacy:` say why a tool is pinned.

- All processing is local; genomic data never leaves the machine
- Pin tool versions; never use floating tags
- Reference genome: GRCh38, NCBI's no-ALT analysis set (`reference/GRCh38_no_alt_analysis_set.fasta`). No aligner here runs ALT-aware, so a reference with ALT contigs gives MAPQ 0 at CYP2D6, the MHC and KIR; `validate-setup.sh` refuses one (unless `ALLOW_ALT_REFERENCE=true`) and refuses a BAM whose @SQ lines differ from the `.fai`. A change of reference means realigning every sample (`docs/realignment.md`)
- License: GPL-3.0

## Tool-Specific Gotchas

### PharmCAT 3.x (pinned: 3.4.0)
- **Two-step workflow**: Preprocessor (`pharmcat_vcf_preprocessor` with `-refFna`) then main jar (`pharmcat.jar`). NOTE: since 3.0 the preprocessor script lost its `.py` extension and the Python package was renamed `preprocessor` → `pcat`. The old `-refFasta` flag is long gone.
- Preprocessor outputs `.preprocessed.vcf.bgz` (NOT `.vcf`).
- **3.x JSON changes vs 2.15.x**: `wildtypeAllele` → `referenceAllele`; the HTML report is **no longer emitted unless `-reporterHtml` is passed explicitly**. The `genes` map may be flat (`{gene -> data}`) or nested (`{source -> {gene -> data}}`). `sourceDiplotypes` (or `recommendationDiplotypes`) carry `allele1`/`allele2` objects with a `.name`. Both CPIC consumers (`scripts/27-cpic-lookup.sh` and `modules/local/cpic_lookup`) run one parser, `bin/pgx_parse.py`: it **auto-detects both shapes**, gives each gene a status (`normal`, `non-normal`, `ambiguous` when the possible diplotypes have different phenotypes, `not called`) and **fails loud**: a recognized report yielding zero genes is reported as a parse failure, never "all genes were successfully called". Guarded by `tests/test_cpic_parser.py` on real 3.2.0 and 3.4.0 reports of the HG002 fixture, `tests/fixtures/pharmcat/report-<version>.json`.
- Pipeline pinned to **3.4.0** (3.4.0 still bundles vcf-parser 0.3.1, so the `##` header rewrite in step 7 and the module stays). Before bumping, revalidate steps 7 and 27 end-to-end against a known sample — JSON structure and preprocessor flags change between major versions — and capture the new version's `report.json` as a parser fixture beside the others (a row in the test's `REAL_REPORTS`).

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
| Ensembl / VEP cache | ~6 months | 13, 30, 23, 31 | ~3 hr |
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
- **Permission denied in container**: the step writes outside the sample directory. Add `--rw DIR` to its `run_in` call; use `--root` only for an image that cannot run as an unprivileged user. A sample directory written by an older version, which ran every container as root, holds root-owned files: `run_in` then prints the `sudo chown -R` that gives it back.
- **0-byte output**: Usually wrong input path inside container. Double-check `:/genome` mount mapping.
- **PCGR/CPSR path confusion**: `--pcgr_dir` should point to PARENT of `data/`, not `data/` itself.
- **VEP cache**: step 13 installs it (resumable, checked against Ensembl's CHECKSUMS); never use VEP's `INSTALL.pl`.

*Generated by [LynxPrompt](https://lynxprompt.com) CLI*


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:6cd5cc61 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->

## Where the tracker syncs

This repo is public, so its tracker syncs only to the private remote named by `sync.remote` in `.beads/config.yaml`. The block above says sync uses "your git remote". Here that never means this GitHub repo. Don't add it as a Dolt remote and don't push `refs/dolt/*` to it.
