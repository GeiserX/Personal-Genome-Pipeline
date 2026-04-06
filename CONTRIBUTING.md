# Contributing to Personal Genome Pipeline

Thank you for your interest in contributing. This pipeline aims to be the most accessible WGS analysis tool for non-bioinformaticians. Every contribution should be evaluated through that lens.

## How to Contribute

### Reporting Bugs

[Open an issue](https://github.com/GeiserX/Personal-Genome-Pipeline/issues/new) with:
- Which step failed (step number and script name)
- Full error message (copy-paste, not screenshot)
- Your platform (OS, Docker version, CPU architecture)
- Input data type (FASTQ, BAM, VCF) and vendor

### Suggesting New Analysis Steps

Before implementing a new step, open an issue to discuss it. Include:
- What the tool does and why it is useful for personal genomics
- A working Docker image (with exact tag) that is publicly available
- Whether it requires additional reference data
- Expected runtime and resource requirements on a 30X WGS sample
- Whether the output is interpretable by a non-expert

### Submitting Pull Requests

1. Fork the repository and create a feature branch
2. Follow the conventions below
3. Test your changes on at least one sample
4. Ensure `shellcheck` passes on all scripts
5. Ensure no personal data leaks: `grep -r '/mnt/user\|internal-host\|/home/' scripts/ docs/`
6. Open a PR with a clear description of what changed and why

## Conventions

### Scripts

Every script must:

```bash
#!/usr/bin/env bash
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
```

- Use `${GENOME_DIR}` for data paths, never hardcoded paths
- Mount as `-v "${GENOME_DIR}:/genome"` in Docker commands
- Always include `--rm --cpus N --memory Xg --user root` in Docker runs
- Validate input files exist before running Docker commands
- Print status messages showing what step is running and where output goes

### Documentation

Every pipeline step needs:
- `docs/NN-tool-name.md` — What it does, why, Docker image, command, output, runtime, notes
- `scripts/NN-tool-name.sh` — The executable script
- Entry in the README step table
- Section in `docs/interpreting-results.md` if the output needs explanation
- Entry in `scripts/validate-setup.sh` for pre-flight checks

### No Personal Data

This repository must never contain:
- Personal file paths (`/mnt/user/`, `/home/username/`, etc.)
- Server hostnames or IP addresses
- Specific sample names as defaults
- Any information that could identify a person's genome

The CI pipeline enforces this with automated scanning.

### Docker Images

- Always specify exact tags (e.g., `staphb/bcftools:1.21`, not `:latest`)
- Exception: images with no versioned tags (e.g., `lgalarno/telomerehunter:latest`)
- Verify the image exists and is publicly pullable before committing
- Document the image in `docs/lessons-learned.md` if there are any gotchas

## Adding a New Pipeline Step

1. **Choose a step number.** Steps 1-20 are taken. New steps should use 21+.
2. **Verify the Docker image works.** Pull it, run it manually on test data, confirm the output.
3. **Create the script** following the template of existing scripts.
4. **Create the documentation** following the template of existing docs.
5. **Update these files:**
   - `README.md` — step table
   - `scripts/run-all.sh` — add to appropriate phase
   - `scripts/validate-setup.sh` — add image check and reference data check
   - `docs/interpreting-results.md` — add output interpretation
   - `docs/00-reference-setup.md` — if new reference data is needed
   - `AGENTS.md` — if the architecture tree changes
6. **Document failures** in `docs/lessons-learned.md` if you hit any issues during development.
7. **Test** on at least one 30X WGS sample.
8. **Open a PR** with all changes in a single commit.

## Code of Conduct

Be respectful. Genomic data is deeply personal. This project exists to empower individuals with their own health data. Keep that mission in mind.

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0.
