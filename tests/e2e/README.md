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
  name order.
- Every case needs a row in `e2e-case-times.tsv`: its stage (`base` or
  `shard`), its measured seconds and its needs. `scripts/ci/e2e-shards.sh
  check` fails when a case has none, and `tests/test_e2e_shards.sh` runs it
  on the real table in Guard's unit tests. The split itself would still run
  a case with no row, as a 300-second shard case that needs nothing, so
  write the seconds the case took.
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
- `bin/docker` lowers `--cpus` to what the runner has; everything else passes
  through to the real Docker.

How the job runs and how to add a case: [docs/testing.md](../../docs/testing.md).
