# End-to-end cases

Each `*.sh` file here (except `lib.sh`) is one case: it runs one pipeline step
on the HG002 fixture with the real containers and checks what the step wrote.
`scripts/ci/e2e-run.sh` runs them all in name order (C locale), keeps going
after a failure, and prints a result table.

- Numbered files are the base cases. Later packages add
  `<package-key>-<what>.sh`; they run after the base cases and can read their
  outputs. Add a file; do not edit the runner, the workflow or another case.
- `lib.sh` holds the shared checks (`check`, `check_eq`, `check_ge`, `has`,
  `lacks`, `vcf_count`, ...). Every check is a count, a column or a content
  match, never an exit code alone.
- `bin/docker` lowers `--cpus` to what the runner has; everything else passes
  through to the real Docker.

How the job runs and how to add a case: [docs/testing.md](../../docs/testing.md).
