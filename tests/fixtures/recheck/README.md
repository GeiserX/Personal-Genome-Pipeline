# Fixtures of bin/pgp_recheck.py

Invented data, read by `tests/test_pgp_recheck.py` and `tests/fake-docker/run-all-recheck.sh`. No real sample.

- `run_manifest.tsv`: the rows `bin/write_manifest.sh` writes that the check reads, for a made-up run.
- `offline.json`: upstream answers in the `--offline` format, with a newer ClinVar, VEP 117 and a newer PGS file.
- `offline-nothing-newer.json`: answers where nothing is newer than the run, for the replay case.
