# tests/fixtures/live/ops

Made-up inputs for `tests/test_pgp_live_ops.py`. No real sample: the run is `OPS-001`.

- `events.jsonl`: Nextflow weblog messages as `bin/pgp_live.py` stores them. A run starts at
  2026-01-02 10:00:00 UTC; step 06's task is killed with exit 137 and its second attempt ends at 7.6 of
  the 8 GB it asked for; step 11's task starts at 10:05 and is still running. `@WORK@` and `@OUT@` are
  replaced by temporary folders in the test.
- `run_status.tsv`: the `logs/run_status.tsv` that `run-all.sh` writes for the same run (steps 06 and 11 planned).
- `step_history.tsv`: a step history of a box called `opsbox` (`bin/run_eta.py`): step 11 took 60 and
  100 seconds in two runs, step 06 took 30.
