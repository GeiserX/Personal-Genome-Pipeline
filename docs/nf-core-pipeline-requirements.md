# nf-core Pipeline Requirements Research Report

> Researched 2026-04-08. Sources: nf-co.re docs, github.com/nf-core/proposals, github.com/nf-core/sarek

---

## 1. Pipeline Submission / Proposal Process

### Where to propose

All new pipeline proposals go through **github.com/nf-core/proposals/issues** using the title format `New pipeline: nf-core/<name>`.

### Pre-requisites before proposing

- Join the nf-core Slack (nf-co.re/join)
- Request GitHub org access via `#github-invitations` Slack channel
- Be familiar with Nextflow, git, and nf-core guidelines
- Ensure the workflow meets nf-core guidelines

### The process

1. **Discuss early** -- "no two pipelines should overlap too much in their purpose and results." If an existing pipeline covers similar ground, you will be asked to contribute to it instead.
2. **Open a proposal issue** at nf-core/proposals with label `new-pipeline`. Status labels: `proposed` (under discussion) then `accepted` (approved by core + maintainers team).
3. **Create the pipeline** using `nf-core pipelines create` under your personal GitHub account first.
4. **Develop with test data** -- configure test profiles.
5. **Transfer to nf-core org** -- after acceptance, the repo moves to github.com/nf-core/.
6. **Make first release** -- semantic versioned tag.
7. **Ongoing maintenance** -- community-owned from this point.

### Current state (Apr 2026)

19 open proposals, ~130 total historical proposals. Recent examples: nf-core/dartseq, nf-core/nanocirc, nf-core/shallowseq. The ecosystem has **147 pipelines** listed.

### If your pipeline already exists externally

- Major breaking changes may be required to meet guidelines
- A new nf-core-style name will likely be needed
- If a similar nf-core pipeline exists (released or in-development), you will be asked to contribute to it instead

### Pipelines NOT suitable for nf-core

Bespoke, proprietary, or highly specialized workflows. These can still use nf-core templates/tools under MIT license as "external" pipelines.

---

## 2. Mandatory Requirements (Must-Have)

These are non-negotiable for any official nf-core pipeline:

| # | Requirement | Detail |
|---|-------------|--------|
| 1 | **Nextflow** | Must be built using Nextflow (DSL2) |
| 2 | **nf-core template** | Must be created via `nf-core pipelines create` and kept in sync |
| 3 | **Community owned** | Pipelines belong to the community, not individuals |
| 4 | **nf-core org hosting** | Primary development on the nf-core GitHub organization |
| 5 | **No overlap** | Only one pipeline per data/analysis type |
| 6 | **Appropriate scope** | "Not too big, not too small" |
| 7 | **Naming** | Lowercase, no punctuation (e.g., `nf-core/mypipeline`) |
| 8 | **MIT license** | All nf-core pipelines must use MIT |
| 9 | **Documentation** | Must be hosted on nf-co.re website |
| 10 | **Docker** | All software bundled in Docker, versioned containers |
| 11 | **CI testing** | Must run CI tests (GitHub Actions) |
| 12 | **Semantic versioning** | Stable release tags (vX.Y.Z) |
| 13 | **Standard parameters** | Standardized parameter naming and usage |
| 14 | **Single command** | Pipeline runs with a single `nextflow run` command |
| 15 | **Minimum inputs** | Run with as little input as possible |
| 16 | **Pass lint** | Zero failures in `nf-core pipelines lint` |
| 17 | **Credits** | Properly acknowledge prior work |
| 18 | **Keywords** | GitHub repo keywords for discoverability |
| 19 | **Git branches** | Must use `master`, `dev`, and `TEMPLATE` branches |
| 20 | **RO-Crate** | Research Object metadata (ro-crate-metadata.json) |

---

## 3. Recommendations (Should-Have)

| Recommendation | Detail |
|----------------|--------|
| **Bioconda** | Package software via bioconda/biocontainers |
| **Modern file formats** | Use CRAM over BAM, etc. |
| **nf-test** | Test pipeline with nf-test using minimal examples |
| **DOIs** | Zenodo DOI for each release |
| **Cloud compatible** | Test on AWS/Azure/GCP |
| **Publication credit** | Acknowledge nf-core in publications |
| **Build with community** | Propose and build collaboratively, not submit finished work |
| **Custom containers** | Only when Bioconda unavailable |

---

## 4. Lint Checks -- All 29 Tests (v3.5.2)

`nf-core pipelines lint` runs these checks. ALL must pass for an official pipeline.

| Test | What it checks |
|------|----------------|
| `actions_awsfulltest` | AWS full test GitHub Actions workflow |
| `actions_awstest` | AWS test GitHub Actions workflow |
| `actions_nf_test` | nf-test CI workflow setup |
| `actions_schema_validation` | Schema validation in CI |
| `base_config` | Base configuration correctness |
| `files_exist` | Required files present |
| `files_unchanged` | Template files not modified |
| `included_configs` | Config files properly included |
| `local_component_structure` | Local modules/subworkflows structure |
| `merge_markers` | No unresolved git merge markers |
| `modules_config` | Module configuration correctness |
| `modules_json` | modules.json validity |
| `modules_structure` | Proper module directory structure |
| `multiqc_config` | MultiQC config file |
| `nextflow_config` | Nextflow config params and settings |
| `nf_test_content` | nf-test file content |
| `nfcore_yml` | .nf-core.yml validity |
| `pipeline_if_empty_null` | Proper ifEmpty/null patterns |
| `pipeline_name_conventions` | Naming conventions |
| `pipeline_todos` | TODO statements in code |
| `plugin_includes` | Nextflow plugin includes |
| `readme` | README content and format |
| `rocrate_readme_sync` | RO-Crate and README sync |
| `schema_description` | Schema param descriptions exist |
| `schema_lint` | JSON schema validity |
| `schema_params` | Schema params match pipeline params |
| `system_exit` | No System.exit calls |
| `template_strings` | No unresolved template strings |
| `version_consistency` | Versions consistent across files |

### Configuring lint exceptions

In `.nf-core.yml`:

```yaml
lint:
  # Disable entire tests
  actions_awsfulltest: False
  pipeline_todos: False
  # Skip specific files within a test
  files_exist:
    - CODE_OF_CONDUCT.md
  files_unchanged:
    - assets/email_template.html
```

### Running lint

```bash
nf-core pipelines lint                          # lint current directory
nf-core pipelines lint --dir <path>             # specific directory
nf-core pipelines lint -k files_exist           # specific test only
nf-core pipelines lint --fix                    # auto-fix (requires clean git)
nf-core pipelines lint --json                   # JSON output
nf-core pipelines lint --show-passed            # show all passed tests
```

---

## 5. Template Structure from `nf-core pipelines create`

### Generated directory tree

```
nf-core-mypipeline/
├── .devcontainer/              # GitHub Codespaces config
├── .github/                    # CI workflows, issue/PR templates
├── .vscode/                    # VSCode settings
├── assets/                     # methods_description_template.yml
├── bin/                        # Helper scripts
├── conf/                       # Config profiles (base, test, test_full)
├── docs/                       # usage.md, output.md, CONTRIBUTING.md
├── modules/                    # nf-core modules directory
│   ├── local/                  # Pipeline-specific modules
│   └── nf-core/                # Installed remote modules
├── subworkflows/               # Composed workflow units
│   ├── local/                  # Pipeline-specific subworkflows
│   └── nf-core/                # Installed remote subworkflows
├── tests/                      # nf-test tests
├── workflows/                  # Main workflow definition
│   └── mypipeline/
│       └── main.nf             # Primary workflow logic
├── .gitignore
├── .nf-core.yml                # nf-core tool configuration + lint exceptions
├── .pre-commit-config.yaml     # Pre-commit hooks
├── .prettierignore
├── .prettierrc.yml
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── LICENSE                     # MIT
├── README.md                   # With badges
├── main.nf                     # Entry point (delegates to workflows/)
├── modules.json                # Tracks installed remote modules
├── nextflow.config             # Main configuration
├── nextflow_schema.json        # Parameter schema (JSON)
├── nf-test.config              # nf-test configuration
├── ro-crate-metadata.json      # Research Object metadata
└── tower.yml                   # Seqera Platform config
```

### Skippable features during creation

The `nf-core pipelines create` command (interactive or via `--template-yaml`) allows skipping:

- `github` -- CI, issue templates (git/gitignore still created)
- `github_badges` -- README badges
- `ci` -- GitHub Actions CI tests
- `changelog` -- CHANGELOG.md
- `license` -- MIT license file
- `test_config` -- test/test_full profiles
- `nf-test` -- nf-test config and CI
- `igenomes` -- iGenomes reference config
- `modules` -- nf-core modules infrastructure
- `multiqc` -- MultiQC integration
- `fastqc` -- FastQC module
- `nf_schema` -- nf-schema validation
- `nf_core_configs` -- Institution config profiles
- `is_nfcore` -- nf-core branding (set false for external pipelines)
- `seqera_platform` -- tower.yml
- `gpu` -- GPU support
- `codespaces` -- Devcontainer
- `vscode` -- VSCode config
- `code_linters` -- Pre-commit + prettier
- `citations` -- CITATIONS.md
- `documentation` -- docs/ directory
- `rocrate` -- RO-Crate metadata
- `email` -- Email notifications

### The .nf-core.yml file

Central configuration for the nf-core tools. Controls:
- Lint exceptions (see section 4)
- Template feature toggles for `nf-core pipelines sync`
- Pipeline metadata

---

## 6. nf-core Modules: Local vs Remote

### Remote modules (from nf-core/modules)

- Hosted in the centralized `github.com/nf-core/modules` repository
- Installed via `nf-core modules install <tool_name>`
- Tracked in `modules.json` with specific git SHA pins
- Live under `modules/nf-core/` in the pipeline
- Available to all nf-core pipelines and the Nextflow community
- ~100+ modules in sarek alone

### Local modules (pipeline-specific)

- Live under `modules/local/` in the pipeline
- Custom processes specific to this pipeline
- Not shared with other pipelines
- Must follow `local_component_structure` lint rules

### Module file structure (per module)

```
modules/nf-core/toolname/
├── main.nf              # Process definition (DSL2)
├── meta.yml             # Metadata, input/output descriptions
└── tests/
    └── main.nf.test     # nf-test unit tests
```

### Key module requirements

- **DSL2 syntax** -- modular, importable processes
- **Meta maps** -- `[meta, file]` tuples for sample metadata tracking
- **Version reporting** -- emit tool version via `eval` output qualifier
- **Every module MUST have a test** (nf-test)
- **meta.yml validated** against JSON schema during linting
- **publishDir** configured in pipeline's `modules.config`, not in the module itself

### Subworkflows

Combine multiple modules into reusable units. Same structure (main.nf, meta.yml, tests). Created via `nf-core subworkflows create`. Live under `subworkflows/nf-core/` (remote) or `subworkflows/local/` (pipeline-specific).

### Creating a new module for nf-core/modules

1. Check it does not already exist (`nf-core modules list`)
2. Fork nf-core/modules, create feature branch
3. Run `nf-core modules create` for boilerplate
4. Fill in main.nf, meta.yml, tests
5. Run `nf-core modules test` (creates snapshots)
6. Run `nf-core modules lint`
7. PR to nf-core/modules, label "Ready for Review", request review from `nf-core/modules-team`

---

## 7. Technical Requirements: Nextflow Version, DSL2, Testing

### Nextflow version

- Sarek (latest, v3.8.1) requires Nextflow >= 25.10.2
- The template pins a minimum Nextflow version in nextflow.config
- Modules require Nextflow >= 21.04.0 (for DSL2 support)

### DSL2

- **Mandatory** -- all nf-core pipelines use DSL2
- DSL1 is no longer supported
- Enables: modular processes, imports, subworkflows, one container per process

### Testing requirements

**nf-test framework:**
- Default testing framework for nf-core (replaces pytest-workflow)
- Pipeline-level and module-level tests
- Snapshot-based validation
- CI runs tests via GitHub Actions

**Assertion patterns:**
- All assertions wrapped in `assertAll()`
- Minimum: process success + version.yml snapshot
- Best practice: snapshot complete `process.out`
- Fallbacks for unstable outputs: file existence, line counts, substring checks, file size
- Gzipped file support via `.linesGzip`
- Directory outputs need explicit handling

**Test profiles:**
- `test` -- minimal dataset, runs on GitHub Actions
- `test_full` -- larger dataset, runs on AWS (full-size benchmarking)

**CI pipeline:**
- Linting (prettier, nf-core lint)
- Pipeline tests (small dataset on GitHub)
- Full tests (larger dataset on AWS)
- Stale issue marking

### Container requirements

- One container per process (DSL2 pattern)
- Docker images must be versioned (no `latest` tag)
- Support Docker, Singularity, and Conda
- Prefer Bioconda + BioContainers for automatic multi-platform images

---

## 8. nf-core/sarek as Reference WGS Pipeline

### Overview

- **Purpose:** Germline + somatic variant calling for WGS/WES/targeted
- **Latest:** v3.8.1 (Feb 2026), Nextflow >= 25.10.2
- **Scale:** ~100+ nf-core modules, 8,515 commits, 555 stars
- **Species:** Any with reference genome (optimized for human/mouse)
- **License:** MIT

### Pipeline steps

1. **Pre-processing:** UMI consensus (fgbio), QC (FastQC), trimming (fastp), contamination removal (BBSplit)
2. **Alignment:** BWA-mem, BWA-mem2, dragmap, Sentieon BWA-mem, Parabricks GPU
3. **BAM processing:** MarkDuplicates, BQSR (GATK4 or Sentieon), stats (samtools, mosdepth)
4. **Variant calling (17 tools):** DeepVariant, GATK HaplotypeCaller, GATK Mutect2, freebayes, Strelka, Manta, TIDDIT, ASCAT, CNVkit, Control-FREEC, Lofreq, MuSE, MSIsensor, Sentieon, indexcov
5. **Filtering:** bcftools (view, norm, isec consensus >= 2 callers), Varlociraptor
6. **Annotation:** SnpEff, Ensembl VEP, BCFtools annotate, SnpSift
7. **QC reporting:** MultiQC

### Architecture patterns to follow

- Root `main.nf` is thin -- delegates to `workflows/sarek/main.nf`
- Subworkflows compose multiple modules (e.g., `vcf_annotate_snpeff`)
- All tools installed as remote nf-core modules (pinned by git SHA in modules.json)
- Local modules only for pipeline-specific logic
- Configuration layered: nextflow.config + conf/ profiles + nextflow_schema.json
- Input: CSV samplesheet (`patient, sample, lane, fastq_1, fastq_2`)
- Tool selection via `--tools` parameter (comma-separated list)

### Benchmarking

Three full-size tests on each release:
- `test_full` -- tumor-normal from SEQ2C consortium
- `test_full_germline` -- WGS 30X GIAB NA12878
- `test_full_germline_ncbench_agilent` -- WES evaluated against truth data via NCBench

---

## 9. Official vs Community/External Pipeline -- Summary

| Aspect | Official nf-core | External/Community |
|--------|------------------|--------------------|
| **Hosting** | github.com/nf-core/ | Any GitHub org |
| **Naming** | `nf-core/<name>` | Must NOT use `nf-core/` prefix |
| **License** | MIT (mandatory) | Any (MIT recommended) |
| **Template** | Required, kept in sync | Optional, can customize |
| **Lint** | Must pass all 29 tests | Can use with exceptions |
| **Docs** | Hosted on nf-co.re | Self-hosted |
| **Modules** | Full access to nf-core/modules | Full access (MIT) |
| **CI** | GitHub Actions (required) | Optional |
| **Logo** | nf-core branding | Must use own branding |
| **Slack** | nf-core Slack references | Must remove nf-core Slack refs |
| **Ownership** | Community-owned | Developer-owned |
| **Proposal** | Required via nf-core/proposals | Not needed |
| **iGenomes** | Included by default | Optional |
| **nf-core configs** | Institutional profiles | Can use (generic options) |

### For external pipelines

- Say the pipeline **"uses"** nf-core, not **"is"** nf-core
- Add acknowledgment citing nf-core Nat Biotechnol 2020 paper
- Can suppress nf-core-specific lint failures via `.nf-core.yml`
- Can still list on nextflow-io/awesome-nextflow

---

## 10. Key Implications for Personal-Genome-Pipeline

### If targeting official nf-core inclusion:

1. Must use MIT license (currently GPL-3.0 -- would need to change)
2. Must create via `nf-core pipelines create` and keep template in sync
3. Must propose at nf-core/proposals first -- risk of rejection if overlap with sarek
4. Pipeline becomes community-owned (loss of control)
5. Must pass all 29 lint tests
6. Must use Docker + Bioconda for all tools
7. Must provide test profiles with minimal test data
8. Must use nf-test for testing

### If building as external pipeline using nf-core tools:

1. Can keep GPL-3.0 license
2. Can use nf-core template as foundation
3. Can install and use all nf-core modules (MIT licensed)
4. Can use nf-core linting with exceptions
5. Must NOT use `nf-core/` prefix
6. Must replace nf-core branding
7. Can still follow all best practices
8. Developer retains full control

### Overlap concern with sarek

Sarek covers the full WGS germline+somatic pipeline. A personal genome pipeline focused on **personal/clinical interpretation** (PRS, pharmacogenomics, HLA typing, ancestry, clinical filtering, STR expansions) would differentiate from sarek's focus on **variant calling**. The post-calling analysis steps (PRS, PGx, HLA, STR, ancestry, clinical reporting) are NOT covered by any existing nf-core pipeline and would be the strongest basis for a proposal.
