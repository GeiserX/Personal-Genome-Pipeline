# Keeping versions up to date

This page is for maintainers. It says how a pinned version gets bumped, what checks the bump, which pins are held on purpose and what still has to be checked by hand.

## Who updates what

| What | Where it is pinned | Who proposes the bump |
|---|---|---|
| Container images | [`versions.env`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/versions.env), copied into `conf/containers.config` and [Image versions](versions.md) | Renovate, one PR per tool |
| Nextflow | `NEXTFLOW_VERSION` in `versions.env`, the `setup-nextflow` pins in the workflows, the `NXF_VER=` lines in the docs | Renovate, after approval on the dependency dashboard |
| actionlint, gitleaks and ShellCheck | `ACTIONLINT_VERSION` and `GITLEAKS_VERSION` in `lint.yml`, and the `version:` of its action-shellcheck step | Renovate |
| Renovate CLI used by the dry run | `RENOVATE_VERSION` in `renovate-dry-run.yml` and the command at the end of this page | Renovate, one PR for both |
| IPD-IMGT/HLA release (step 8) | `HLA_DB_RELEASE` in `scripts/lib/common.sh` | Renovate, from the `v<x.y.z>-alpha` tags of `ANHIG/IMGTHLA` |
| pgsc_calc (steps 25 and 26) | `PGSC_CALC_VERSION` in `versions.env`, copied into `conf/containers.config` | Renovate, after approval on the dependency dashboard. `PGSC_CALC_SHA256`, `PGSC_CALC_NF_SCHEMA` and the `PGSC_*_IMAGE` lines move with it by hand ([PRS](25-prs.md#maintenance)) |
| Python of the docs build | `python-version` in `docs.yml` | Renovate, after approval on the dependency dashboard |
| Python packages of the screenshot script | the `==` pins in `tests/demo/retake-screenshots.sh` | Renovate |
| Aldy, the alternative CYP2D6 caller | its version in [Cyrius](21-cyrius.md), in the prose and the `pip install` line | Renovate |
| GitHub Actions (`uses:` lines) | the workflows | Dependabot, monthly |
| MkDocs and its plugins | `docs/requirements-docs.txt` | Dependabot, monthly |
| Data releases: VEP cache, PCGR bundle, pypgx bundle, AnnotSV annotations, Yleaf data, gnomAD constraint, the GRIDSS blacklist commit | the variables at the end of `versions.env` | nobody: bump them by hand, together with the image they belong to. Guard's `coupled-versions` job fails while the VEP cache, pypgx bundle, AnnotSV or Yleaf data disagrees with its image |
| GENCODE gene lines (step 8) | `GENCODE_RELEASE` in `scripts/lib/common.sh` | nobody: it moves with a VEP major, to the GENCODE release of that Ensembl release (the VEP major minus 66). The `coupled-versions` job fails until it does |
| VerifyBamID2 marker panel (step 33) | the tag in the `VB2_PANEL_FILES` URLs of `scripts/setup.sh`, with their sha256 | nobody: it moves with `VERIFYBAMID2_IMAGE`. The `coupled-versions` job fails while the tag and the image version differ |
| Delly exclude list (step 19) | `DELLY_EXCL_COMMIT` in `scripts/lib/common.sh` | nobody: a commit of dellytools/delly, moved by hand |
| IPD-KIR release (step 8 with `KIR=true`) | `KIR_DB_RELEASE` in `versions.env` | nobody: the monthly freshness report says when a newer release exists |
| Parascopy data (step 35) | `PARASCOPY_DATA_VERSION` in `versions.env`, with the Zenodo record and the md5s in `scripts/setup.sh` | nobody: the monthly freshness report says when the record has a newer version |
| Cyrius lock (step 21) | `CYRIUS_VERSION` in `versions.env` and every wheel in `scripts/cyrius-constraints.txt` | nobody: resolved again by hand in `PYTHON_IMAGE`, as the lock's header says. The freshness report says when PyPI has a newer Cyrius |
| The site table of step 40: the 23andMe v5 site list commit, the dbSNP build and the table's hash | `CONSUMER_SITES_COMMIT`, `CONSUMER_DBSNP_BUILD`, their checksums and `CONSUMER_SITES_TABLE_SHA256` in `versions.env` | nobody: dbSNP releases a build every year or two, and the site list has not changed since 2019. A new build, a new commit or a change to `bin/export_consumer.py build-sites` needs one full build (`setup.sh --consumer-sites`, about 30 GB), then `gzip -dc <table> \| sha256sum` into `CONSUMER_SITES_TABLE_SHA256`. CI never builds the full table. The monthly freshness report says when NCBI publishes a newer build |

Three pins are floors, not versions that drift: Node 24 in `renovate-dry-run.yml`, Java 17 in the workflows that run Nextflow, and `manifest.nextflowVersion` in `nextflow.config`. Nothing proposes a bump for them, on purpose.

Renovate runs only its regex managers ([`renovate.json`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/renovate.json)), so it never touches what Dependabot owns ([`.github/dependabot.yml`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/dependabot.yml)).

## Turning Renovate on

The Renovate GitHub App is installed on this repository only and runs in Interactive mode: it keeps a dependency dashboard issue and opens PRs on the schedule below. Majors (VEP) and every Python, Nextflow, pgsc_calc and PCGR update wait for a tick on the dashboard; ticking a box under "Awaiting Schedule" opens that PR at once. The mode and the app's repository list are set on [developer.mend.io](https://developer.mend.io/) and in the GitHub app settings; the dry run workflow below checks the config before Renovate reads it.

## What a Renovate PR looks like

- It is opened on a Monday, before 06:00 UTC, with the labels `dependencies` and `automated`.
- One tool per PR. No group rules.
- The title starts with `deps:` (`ci:` for actionlint, gitleaks, ShellCheck and the Renovate CLI).
- It changes the tool's line in `versions.env` and the same string in `conf/containers.config` and `docs/versions.md`. Those two generated files are written by `scripts/ci/gen-containers-config.sh` and `scripts/ci/gen-versions-doc.sh`, which copy the string from `versions.env` as it is, so the PR leaves them exactly as the scripts would. If a check still says one of them is stale, run both scripts on the branch and commit the result.
- For an image, Renovate does not edit version strings outside those three files. One is checked by CI, so these PRs go red until you change it on the branch by hand:
    - bcftools and VEP: the two tags written in the `--check-versions` self-test of `scripts/ci/gen-containers-config.sh`, which Guard runs against the real `versions.env`.
- A major update (VEP 116 to 117, for example), every PCGR update, every Python update, every Nextflow update and every pgsc_calc update waits on the dependency dashboard. Nothing is opened until you tick its box there.

Renovate reads biocontainer tags (`1.0.9--h5ca1c30_0`) as version, then build number, and ignores the conda build hash in between, which changes from one version to the next. It keeps the Python build (`py312` in `3.9.5--py312h6e8b409_0`): only tags for the same Python are proposed. A new build of the same version (mosdepth `0.3.14--h05c3d44_0` to `0.3.14--h87be163_2`) is proposed as a patch: same tool version, rebuilt by Bioconda. A biocontainer tag in any other shape is never proposed. VEP tags (`release_116.0`) are read as major and minor.

## What checks a bump

Every PR that changes `versions.env` runs:

- **Guard**: `conf/containers.config` matches `versions.env` and every Nextflow process has an image.
- **Renovate dry run**: every image is still found in every file it is copied to, and `docs/versions.md` is exactly what `scripts/ci/gen-versions-doc.sh` writes from `versions.env`. When it is not, the job prints the command to run.
- **Container test** (`container-test.yml`): runs each image whose line changed on the fixture, with its row of `tests/smoke/commands.tsv`, and checks what the tool wrote. Every image has a row; `tests/smoke/README.md` names the few that cannot run their real job on a runner (VEP's offline cache, the CPSR data bundle). Its `pinned-images` job runs `setup.sh --pull-only` on the images pinned by digest and checks that none is left dangling, so `docker image prune` keeps them, on Docker's classic image store and on its containerd image store.
- **E2E**: runs the real tools on a small slice of a public genome. [Testing](testing.md#the-e2e-job) lists the steps it covers and the ones no CI job can run (offline VEP, CPSR, AnnotSV, GRIDSS).

A bump whose real job neither the container test nor the E2E job can run (offline VEP, CPSR) is not tested by CI. Run that step by hand on a sample before merging.

## Held pins

Each hold is a rule in `renovate.json` with its reason in the rule's `description`. Most of them are also noted on the tool's line in `versions.env`.

| Pin | Rule | Why | When it ends |
|---|---|---|---|
| pypgx 0.26.0 | only `0.26.0--` tags | 0.27.0 pulled pandas 3.0 and broke every gene while its smoke test still passed | upstream releases the pandas fix and a pypgx-bundle tag of the same version exists; move `PYPGX_BUNDLE_VERSION` with it |
| plink2 2.00a5.10 | disabled | the only versioned tag `pgscatalog/plink2` ships, and the build pgsc_calc uses | a newer tag from that publisher |
| Manta, Strelka2, duphold, Octopus, GRIDSS | disabled | legacy: archived or quiet upstream | bump by hand after a run on the fixture |
| minimap2+samtools image (`MINIMAP2_SAMTOOLS_IMAGE`) | disabled | BioContainers tags it with a hash, not a version; the Nextflow alignment task needs both tools in one image | never by itself: it moves by hand with `MINIMAP2_IMAGE` and `SAMTOOLS_IMAGE`. `scripts/ci/check-images.sh` fails a bump of either until the versions comment of `MINIMAP2_SAMTOOLS_IMAGE` names the new release. Request the new pair in [BioContainers/multi-package-containers](https://github.com/BioContainers/multi-package-containers) (`combinations/hash.tsv`), pin the tag it publishes and update that comment; the container test runs both tools' version commands in it |
| Sniffles 2.8.0 | none needed | 2.8.1 has no biocontainer tag in the normal shape yet | Renovate proposes it once the tag exists |
| Yleaf 3.2.1 | none needed | the only Yleaf biocontainer; upstream is on 4.x. `YLEAF_DATA_VERSION` must equal the image version | a 4.x biocontainer, after step 37 runs on it; move `YLEAF_DATA_VERSION` and `YLEAF_DATA_SHA256` with it |
| hap.py 0.3.12 (`HAPPY_IMAGE`) | none needed | the hap.py biocontainer (0.3.15) has no RTG Tools, so it cannot run `--engine=vcfeval` | a biocontainer that ships RTG Tools |
| ADMIXTURE, RFMix, Beagle (step 26b) | disabled | their tags are not versions Renovate can read (`1.3.0--0`, a commit, a date) | bump by hand after step 26b runs on the e2e fixture; the freshness report shows a newer tag |
| pgsc_calc images (`ghcr.io/pgscatalog/...`) | disabled | they are the images pgsc_calc's `conf/modules.config` names at `PGSC_CALC_VERSION` | never by themselves: they move with a pgsc_calc release |

## What to check by hand

- **PharmCAT**: diff the diplotype table on the HG002 fixture between the old and the new image, and run the CPIC parser tests on a `report.json` from the new version.
- **VEP major**: move `VEP_CACHE_RELEASE` to the new major and `GENCODE_RELEASE` (`scripts/lib/common.sh`) to the new major minus 66 in the same PR, download the new cache (about 26 GB) and rerun step 13. CI cannot run the offline cache.
- **PCGR/CPSR**: move `PCGR_DATA_BUNDLE` and `PCGR_VEP_CACHE_RELEASE` with the image, and rerun step 17 on a sample with an earlier result. CI does not run CPSR on a sample. The image test checks that `cpsr --help` still lists every flag step 17 passes, that `cpsr --version` names the release line of the image tag and the `pcgrr` R package its full release, and runs the bundled VEP with `--database` on 50 variants.
- **Python**: resolve `scripts/cyrius-constraints.txt` again in the new image, as its header says, and rerun step 21 on a test BAM.
- **Nextflow**: the PR runs the stub and E2E jobs on the new version. The prose names `NEXTFLOW_VERSION` rather than a number, and the fake nextflow of the fake-docker suite reads it from `versions.env`, so a patch needs no hand edit. A new release line (26.10 to 27.04, say) also needs the line named in `nextflow.config`, `nextflow.yml` and the prose of [Nextflow Execution](nextflow.md) and [Lessons Learned](lessons-learned.md) moved by hand.
- **pypgx, when the hold ends**: check out the matching pypgx-bundle tag and compare the gene calls on the fixture.
- **Images pinned by digest** (`name:tag@sha256:...`): after the bump, run `./scripts/setup.sh --pull-only`. On Docker's classic image store an image pulled by digest has no tag, so it is dangling and `docker image prune` deletes it (the containerd image store keeps the digest reference as its name, so it is not dangling there); `setup.sh` gives it a local tag, `localhost/pgp-pinned/<name>:sha256-<digest>`, and a pruned one is restored the same way. Tasks still run the digest reference. The old image keeps its own tag, so a prune no longer removes it: list the tags with `docker images --format '{{.Repository}}:{{.Tag}}' | grep '^localhost/pgp-pinned/'` and `docker rmi` the ones whose digest is no longer in `versions.env`.

## Checking a change to the Renovate config

The [Renovate dry run](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/workflows/renovate-dry-run.yml) workflow runs on every PR that changes `renovate.json`, `versions.env`, `conf/containers.config`, a page in `docs/`, a workflow, `scripts/ci/gen-versions-doc.sh`, `scripts/lib/common.sh` or `tests/demo/retake-screenshots.sh`, and by hand. It runs Renovate's lookups with no app and no write access (`--platform=local --dry-run=lookup`), validates `renovate.json`, and fails when:

- an `*_IMAGE` line of `versions.env` is not a detected dependency;
- a package name carries a `:` or `@`;
- `docs/versions.md` does not hold exactly the images of `versions.env`, or `conf/containers.config` holds one it does not;
- a tool would get a different update in one file than in another;
- a pin that is not an image (Nextflow, actionlint, gitleaks, ShellCheck, the Renovate CLI, the HLA database, pgsc_calc, the screenshot script's packages, Aldy, the docs build's Python) is not detected in every place the text pins it;
- a lookup failed, or a dependency was skipped for any reason other than a rule that disables it.

The lookup job's summary lists every dependency with the update Renovate would propose, or why it proposes none. The full debug log is the `renovate-dry-run-log` artifact. A second job, `versions-doc`, fails when `docs/versions.md` is not what `scripts/ci/gen-versions-doc.sh` writes from `versions.env`.

To run the same lookups on your machine, with Node 24:

```bash
LOG_LEVEL=debug npx --yes renovate@44.166.0 --platform=local --dry-run=lookup
```
