# Families fixtures

Invented names only, for `tests/test_family_report.py`.

- `two_families.tsv`: a families file with the family column, two families
  (a trio and a couple).
- `members.golden.tsv`: the `members.tsv` that `family_report.py prepare`
  wrote for `assets/pedigree.example.tsv` before the family column existed. A
  file without the column must still give these bytes.
- `stub_families.tsv`: two families for the stub run of
  `tests/e2e/families-file-1.sh`; family_b holds two samples of
  `assets/test_samplesheet.csv`, family_a none.
