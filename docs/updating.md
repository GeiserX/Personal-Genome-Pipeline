# Keeping versions up to date

This page is for maintainers. It says how a pinned version gets bumped, what checks the bump, which pins are held on purpose and what still has to be checked by hand.

## Who updates what

| What | Where it is pinned | Who proposes the bump |
|---|---|---|
| Container images | [`versions.env`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/versions.env), copied into `conf/containers.config`, [Image versions](versions.md) and the matrix of `container-test.yml` | Renovate, one PR per tool |
| Nextflow | `NEXTFLOW_VERSION` in `versions.env`, the `setup-nextflow` pins in the workflows, the `NXF_VER=` lines in the docs | Renovate, after approval on the dependency dashboard |
| actionlint and gitleaks | `ACTIONLINT_VERSION` and `GITLEAKS_VERSION` in `lint.yml` | Renovate |
| GitHub Actions (`uses:` lines) | the workflows | Dependabot, monthly |
| MkDocs and its plugins | `docs/requirements-docs.txt` | Dependabot, monthly |
| Data releases: VEP cache, PCGR bundle, pypgx bundle, AnnotSV annotations, gnomAD constraint, Cyrius, the GRIDSS blacklist commit | the variables at the end of `versions.env` | nobody: bump them by hand, together with the image they belong to |

Renovate runs only its regex managers ([`renovate.json`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/renovate.json)), so it never touches what Dependabot owns ([`.github/dependabot.yml`](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/dependabot.yml)).

## Turning Renovate on

The Renovate GitHub App is installed on the repository in Silent mode: it reads the config and opens nothing. The one owner step is to switch the repository to Interactive on [developer.mend.io](https://developer.mend.io/). From then on Renovate opens a dependency dashboard issue and its PRs.

Until that switch, the Renovate dry run workflow below is the only evidence that the config finds the pins and the updates.

## What a Renovate PR looks like

- It is opened on the first day of the month, before 06:00 UTC, with the labels `dependencies` and `automated`.
- One tool per PR. No group rules.
- The title starts with `deps:` (`ci:` for actionlint and gitleaks).
- It changes the tool's line in `versions.env` and the same string in `conf/containers.config`, `docs/versions.md` and the `container-test.yml` matrix. Those two generated files are written by `scripts/ci/gen-containers-config.sh` and `scripts/ci/gen-versions-doc.sh`, which copy the string from `versions.env` as it is, so the PR leaves them exactly as the scripts would. If a check still says one of them is stale, run both scripts on the branch and commit the result.
- A major update (VEP 116 to 117, for example), every PCGR update and every Nextflow update waits on the dependency dashboard. Nothing is opened until you tick its box there.

Renovate reads biocontainer tags (`1.0.9--h5ca1c30_0`) as version, then build number, and ignores the conda build hash in between, which changes from one version to the next. A biocontainer tag in any other shape is never proposed. VEP tags (`release_116.0`) are read as major and minor.

## What checks a bump

Every PR that changes `versions.env` runs:

- **Guard**: `conf/containers.config` matches `versions.env` and every Nextflow process has an image.
- **Renovate dry run**: every image is still found in every file it is copied to.
- **Container Smoke Test**: pulls each image in its matrix and runs a short command. The matrix is a subset of `versions.env`.
- **E2E**: runs the real tools on a small slice of a public genome. [Testing](testing.md#the-e2e-job) lists the steps it covers and the ones no CI job can run (offline VEP, CPSR, AnnotSV, GRIDSS).

A bump of a tool that neither the smoke test matrix nor the E2E job runs is not tested by CI. Run that step by hand on a sample before merging.

## Held pins

Each hold is a rule in `renovate.json` with its reason in the rule's `description`. Most of them are also noted on the tool's line in `versions.env`.

| Pin | Rule | Why | When it ends |
|---|---|---|---|
| pypgx 0.26.0 | only `0.26.0--` tags | 0.27.0 pulled pandas 3.0 and broke every gene while its smoke test still passed | upstream releases the pandas fix and a pypgx-bundle tag of the same version exists; move `PYPGX_BUNDLE_VERSION` with it |
| Delly 2.1.0 | disabled | 2.3.0 renamed `delly call` to `delly sr`; 2.6.0 has the pin's build hash, so a bump would look safe and break step 19 | the bump that changes `scripts/19-delly.sh` and the Delly module removes the rule |
| plink2 2.00a5.10 | disabled | the only versioned tag `pgscatalog/plink2` ships, and the build pgsc_calc uses | a newer tag from that publisher |
| Manta, Strelka2, duphold, Octopus, GRIDSS | disabled | legacy: archived or quiet upstream | bump by hand after a run on the fixture |
| TelomereHunter, haplogrep3 | disabled | personal-account images pinned by digest, with no versioned tags | they move to biocontainers, and the rule goes with them |
| fastp 1.3.6, Sniffles 2.8.0 | none needed | 1.3.7 and 2.8.1 have no biocontainer tag in the normal shape yet | Renovate proposes them once the tag exists |

## What to check by hand

- **PharmCAT**: diff the diplotype table on the HG002 fixture between the old and the new image, and run the CPIC parser tests on a `report.json` from the new version.
- **VEP major**: move `VEP_CACHE_RELEASE` to the new major in the same PR, download the new cache (about 26 GB) and rerun step 13. CI cannot run the offline cache.
- **PCGR/CPSR**: move `PCGR_DATA_BUNDLE` and `PCGR_VEP_CACHE_RELEASE` with the image, and rerun step 17 on a sample with an earlier result. CI does not run CPSR.
- **Nextflow**: run the stub and E2E jobs on the new version, then update the validated version in `nextflow.config` and the prose in [Nextflow Execution](nextflow.md) and [Lessons Learned](lessons-learned.md).
- **pypgx, when the hold ends**: check out the matching pypgx-bundle tag and compare the gene calls on the fixture.

## Checking a change to the Renovate config

The [Renovate dry run](https://github.com/GeiserX/Personal-Genome-Pipeline/blob/main/.github/workflows/renovate-dry-run.yml) workflow runs on every PR that changes `renovate.json`, `versions.env`, one of the files copied from it or the workflow itself, and by hand. It runs Renovate's lookups with no app and no write access (`--platform=local --dry-run=lookup`), validates `renovate.json`, and fails when:

- an `*_IMAGE` line of `versions.env` is not a detected dependency;
- a package name carries a `:` or `@`;
- `docs/versions.md` does not hold exactly the images of `versions.env`, or `conf/containers.config` or the smoke test matrix holds one it does not;
- a tool would get a different update in one file than in another;
- a lookup failed, or a dependency was skipped for any reason other than a rule that disables it.

Its job summary lists every dependency with the update Renovate would propose, or why it proposes none. The full debug log is the `renovate-dry-run-log` artifact.

To run the same lookups on your machine, with Node 24:

```bash
LOG_LEVEL=debug npx --yes renovate@44.132.2 --platform=local --dry-run=lookup
```
