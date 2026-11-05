# tests/fixtures/live/demo

Invented inputs of the live-site demo (`tests/demo/build_live_demo.py`), added to the invented
sample DEMO-001 that `tests/demo/make_demo_sample.py` writes. No real sample, person or run.

- `DEMO-001/recheck/state.json`: a database re-check as `bin/owngenome_recheck.py` writes it: a newer
  ClinVar release and an Ensembl cache with nothing newer. The build sets `checked_utc` to one hour
  after the replayed run ends and each item's `steps` with `owngenome_recheck.affected_steps`, so the
  badges name the steps the real check would name.
- `DEMO-FAMILY/members.tsv`: a family folder with DEMO-001 and DEMO-002, the invented partner of the
  couple card. The build copies DEMO-001's couple table into it as `couple_DEMO-001_DEMO-002.tsv`, as
  step 39 names it.

The run itself is `../demo-001.weblog.jsonl`, written by `tests/demo/make_demo_events.py`.
