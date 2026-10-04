# PharmCAT report fixtures

Real PharmCAT output for `tests/test_cpic_parser.py` (one `report-<version>.json` per release in its `REAL_REPORTS`), trimmed with
`trim_report.py` to the fields `bin/pgx_parse.py` reads (every gene with all
its listed diplotypes and `relatedDrugs`; up to 4 matched annotations per drug
guideline).

| File | Source |
|---|---|
| `report-3.2.0.json` | PharmCAT 3.2.0 (`pgkb/pharmcat:3.2.0`) on the public GIAB HG002 slice of the e2e fixture (`fixture-v4`), written by step 7 in E2E run 37054440069 of pull request 72 (the `e2e-logs` artifact, `HG002.pharmcat-report.json`, 40 MB, sha256 `562bab215132c9da8a700c810b1a11f56b73aad5815d3e8c6609057b359ebcd2`). The slice covers few pharmacogenes: most are `Unknown/Unknown`, and CYP2C19 (528 possible diplotypes) and CYP2B6 (4) are ambiguous. |
| `report-3.4.0.json` | PharmCAT 3.4.0 (`pgkb/pharmcat:3.4.0`) on the same e2e fixture (`fixture-v4`), written by step 7 in E2E run 37156576994 of pull request 83 (`e2e-logs`, `HG002.pharmcat-report.json`, 40 MB, sha256 `48a2f35f8fa81d28ea06e09dd2ef2c9caf01d0ddb474af0fad975feee54a421d`). Every gene has the same possible diplotypes and phenotypes as in the 3.2.0 report: CYP2C19 528, CYP2B6 4, the rest `Unknown/Unknown`. |
| `pharmcat-docs-example.json` | `docs/examples/pharmcat.example.report.json` of the PharmCAT repository at tag `v3.2.0` (sha256 `23d58e7b2bf815159bb6136beec919d9c6eced0fe71cac7c38cebff6b4febf05`), PharmCAT's own example sample: 23 genes with one diplotype each, 10 of them with a non-normal phenotype. |

To refresh one or add a release: take the new `report.json`, then run this
with `VERSION` replaced by the PharmCAT release (for example `3.4.0`), and add
the release to `REAL_REPORTS` in `tests/test_cpic_parser.py`:

```bash
python3 tests/fixtures/pharmcat/trim_report.py report.json tests/fixtures/pharmcat/report-VERSION.json
```
