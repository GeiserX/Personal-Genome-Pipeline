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

Every script starts like this:

```bash
#!/usr/bin/env bash
set -euo pipefail

SAMPLE=${1:?Usage: $0 <sample_name>}
GENOME_DIR=${GENOME_DIR:?Set GENOME_DIR to your data directory}
# shellcheck source=lib/common.sh
. "$(dirname "$0")/lib/common.sh"
validate_sample "$SAMPLE"
```

[`scripts/lib/common.sh`](scripts/lib/common.sh) reads [`versions.env`](versions.env) and gives the script its image variables, `REF_FASTA`, `THREADS` and the helpers below.

- Start containers with `run_in --cpus N --memory Xg "${TOOL_IMAGE}" tool ...`. It runs `docker run --rm` with no network, `GENOME_DIR` mounted read-only at `/genome`, the sample directory writable, and the calling user.
- Say so at the call when a step needs more, with the reason in a comment: `--rw DIR` to write a shared index or database, `--net` to download, `--root` for an image that cannot run as an unprivileged user.
- Name images by their `versions.env` variable, quoted. A script never spells an image name or tag.
- Use `${REF_FASTA}` on the host and `${REF_FASTA_C}` inside a container for the reference.
- Download with `fetch URL DEST [md5|sha256|sum VALUE-or-URL]`. It writes `DEST.part`, checks it, then renames it.
- Use `${GENOME_DIR}` for data paths, never hardcoded paths
- Validate input files exist before starting a container
- Print status messages showing what step is running and where output goes

### Documentation

Every pipeline step needs:
- `docs/NN-tool-name.md` — What it does, why, Docker image, command, output, runtime, notes
- `scripts/NN-tool-name.sh` — The executable script
- An entry in the `nav:` of [`mkdocs.yml`](mkdocs.yml), or the docs build fails on the orphan page
- A mention in the category table of [`docs/pipeline-overview.md`](docs/pipeline-overview.md)
- Section in `docs/interpreting-results.md` if the output needs explanation

### No Personal Data

This repository must never contain:
- Personal file paths (`/mnt/user/`, `/home/username/`, etc.)
- Server hostnames or IP addresses
- Specific sample names as defaults
- Any information that could identify a person's genome

The CI pipeline enforces this with automated scanning.

### Docker Images

- Every image is one line in [`versions.env`](versions.env), with an exact tag (e.g., `staphb/bcftools:1.21`, not `:latest`). Scripts, `setup.sh` and `validate-setup.sh` read it from there.
- Mark an image no default step runs with `# optional` on its line. `setup.sh` then leaves it for the step that uses it.
- When a publisher offers no versioned tags, pin by immutable digest (`name@sha256:<digest>`) — never a floating `:latest`. Resolve with `docker manifest inspect -v <name>:latest`.
- Verify the image exists and is publicly pullable before committing
- Document the image in `docs/lessons-learned.md` if there are any gotchas

## Adding a New Pipeline Step

1. **Choose a step number.** Steps 1-32 are taken. New steps should use 33+.
2. **Verify the Docker image works.** Pull it, run it manually on test data, confirm the output.
3. **Add the image to [`versions.env`](versions.env)** as `TOOL_IMAGE="name:tag"`. That one line makes `setup.sh` pull it and `validate-setup.sh` check it.
4. **Create the script** `scripts/NN-tool-name.sh` following the conventions above.
5. **Create the Nextflow module** `modules/local/<tool>/main.nf` and wire it into the workflow, or say in [`docs/nextflow.md`](docs/nextflow.md) why the step stays bash-only. A module has no `container` line: add its process to the table in `scripts/ci/gen-containers-config.sh` and run that script.
6. **Create the documentation** `docs/NN-tool-name.md` following the existing step docs.
7. **Update these files:**
   - [`mkdocs.yml`](mkdocs.yml): the page in `nav:`
   - [`docs/pipeline-overview.md`](docs/pipeline-overview.md): the category table
   - `scripts/run-all.sh`: add to the appropriate phase
   - [`.github/workflows/container-test.yml`](.github/workflows/container-test.yml): the image and a smoke command (`tool --version`) in the matrix
   - `scripts/validate-setup.sh`: only if the step needs reference data to check
   - `docs/interpreting-results.md`: add output interpretation
   - `docs/00-reference-setup.md`: if new reference data is needed
   - `CLAUDE.md` — if the architecture tree changes
8. **Document failures** in `docs/lessons-learned.md` if you hit any issues during development.
9. **Test** on at least one 30X WGS sample. CI runs the default steps on a small fixture, see [`docs/testing.md`](docs/testing.md).
10. **Open a PR** with all changes.

## Bumping a Tool

1. Change the tool's line in [`versions.env`](versions.env). No script or module names the tag.
2. If a comment next to the line couples it to a data version (the VEP cache release, the PCGR bundle, the pypgx bundle tag), change that variable in the same commit.
3. Run `scripts/ci/gen-containers-config.sh`. It rewrites [`conf/containers.config`](conf/containers.config), where the Nextflow modules take their image from, and CI fails until that file matches `versions.env`.
4. Update the matrix entry in [`.github/workflows/container-test.yml`](.github/workflows/container-test.yml) and any doc that prints the tag. The `version-consistency` check in CI names a stale tag in the docs.
5. Run `./scripts/setup.sh --pull-only` to pull the new image.
6. A line marked `hold:` or `legacy:` says why the tool is pinned and when the hold ends. Read it before bumping.

## Code of Conduct

Be respectful. Genomic data is deeply personal. This project exists to empower individuals with their own health data. Keep that mission in mind.

## License

By contributing, you agree that your contributions will be licensed under GPL-3.0.
