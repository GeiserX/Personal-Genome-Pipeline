# End-to-end cases

Each `*.sh` file here (except `lib.sh`) is one case: it runs one pipeline step
on the HG002 fixture with the real containers and checks what the step wrote.
`scripts/ci/e2e-run.sh` runs cases in name order (C locale), keeps going
after a failure, and prints a result table. In CI the cases are split over
parallel jobs by `scripts/ci/e2e-shards.sh`, and each case runs at most once.

- Numbered files are the base cases. Later packages add
  `<package-key>-<what>.sh`; they run after the base cases and can read their
  outputs. Add a file; do not edit the runner, the workflow or another case.
- The base cases of [`scripts/ci/e2e-case-times.tsv`](../../scripts/ci/e2e-case-times.tsv)
  (stage `base`: the numbered files and `bash-step-08-hla`) run first, in
  their own job, and every shard starts from what they wrote. A case that reads what another, non-base case wrote must
  name that case in its `needs` column, so both run in the same shard, in
  name order. A new case with no row runs in some shard, counts as 300
  seconds and needs nothing. A case that runs well over 300 seconds also
  needs a row with its measured seconds, or the split balances it as a
  300-second case.
- Every case needs a row in
  [`scripts/ci/e2e-case-paths.tsv`](../../scripts/ci/e2e-case-paths.tsv):
  the globs of the files it runs (step scripts, modules, bin helpers,
  assets, a committed file under `tests/fixtures/`), or `-` for a base
  case. A fixture a case reads also goes in `SCOPE_FILES` of
  `scripts/ci/e2e-scope.sh` and in the push paths of `e2e.yml`. A pull request runs only the cases whose
  files, rows or paths it changed, the cases those need and two smoke cases;
  a change to shared wiring such as `main.nf` runs them all
  (`scripts/ci/e2e-scope.sh`). A case with no row turns `e2e` red.
- `lib.sh` holds the shared checks (`check`, `check_eq`, `check_ge`, `has`,
  `lacks`, `vcf_count`, ...). Every check is a count, a column or a content
  match, never an exit code alone.
- A case gets 3600 seconds (`CASE_TIMEOUT` of `scripts/ci/e2e-run.sh`). One
  that needs longer says so in a line `# case-timeout: SECONDS` of its file.
- `bin/docker` lowers `--cpus` to what the runner has; everything else passes
  through to the real Docker.

How the job runs and how to add a case: [docs/testing.md](../../docs/testing.md).
