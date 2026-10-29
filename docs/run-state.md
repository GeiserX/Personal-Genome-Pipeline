# Where a Run Stands

`bin/run_state.py` tells you, during a run or after it, which steps are planned, which are running and which have finished. It writes one file, `<sample>/run_state.json`. It reads only the files that say how the run goes, never a result.

```bash
python3 bin/run_state.py --sample <sample> --sample-dir "$GENOME_DIR/<sample>"
```

Standard library only, so the plain `python3` on the host is enough. Run it as often as you like: each call reads the files again and replaces the JSON.

## What it reads

| File | Written by | What it says |
|---|---|---|
| `<sample>/logs/run_status.tsv` | `scripts/run-all.sh` | When the run started; a `skipped (<reason>)` row for each step it does not start; a `planned` row for each step it starts, written just before the pipeline starts; each step's result once the pipeline has ended; `meta finished_utc` when `run-all.sh` ends |
| `<sample>/nextflow/process_status.tsv` | `main.nf`, as the pipeline ends | Per process: tasks that succeeded, came from the cache, failed, were ignored or were left unfinished |
| `GENOME_DIR/pipeline_info/trace_*.txt` | Nextflow | One row per finished task. Only the trace files written since the run started are read (or the files you name with `--trace`), and only the rows tagged with this sample |
| an events file (`--events FILE`) | Nextflow's weblog messages, one JSON object per line | Tasks submitted, started and completed, while the run goes |
| `<sample>/summary.json` | the reports | Only each section's `source` path, to tell an earlier run's file from one this run wrote |

When `run-all.sh` ends, however it ends, a `planned` row gives way to the step's result, so `run_status.tsv` again holds one row per step. A step still `planned` at that point, for example a script step after a failed pipeline, becomes `not finished`.

## The steps catalogue

`assets/steps.tsv` has one row per step: its number, its label as `run-all.sh` prints it, its `--tools` name, whether `run-all.sh` records it in `run_status.tsv`, its Nextflow processes, the report sections it feeds, the steps whose output it reads, and its usual hours. The hours are the upper end of the step's time in the [runtime table](hardware-requirements.md#runtime-per-step); `hours_from` says whether that figure was measured or is the page's estimate, and is empty when the table has no time for the step.

`tests/test_run_state.py` keeps the catalogue true: every step `run-all.sh` plans has its row with the same label and tool, the recorded steps are the ones `bin/collect_summary.py` checks for stale results, every report section belongs to the row of its step, every process in `modules/local` belongs to a step, and the tools are the ones `main.nf` accepts.

## What a step's state means

| State | Meaning |
|---|---|
| `waiting` | Planned, no task of it has started |
| `running` | A task of it has started, or some have finished and others are still to come |
| `done` | Its tasks finished well, or `run-all.sh` recorded it `ok` |
| `failed` | A task of it failed and was not retried with success, or `run-all.sh` recorded it failed |
| `skipped` | `run-all.sh` did not start it; `reason` says why |
| `not finished` | The run ended before the step did |
| `not planned` | No row and no task: a `nextflow run` without `run-all.sh` that has not shown this step |

A result row of `run-all.sh` always wins over what the tasks say. Before that row exists, `process_status.tsv` decides, and before that the trace and the events.

A `waiting` or `running` step whose report file is older than the run start has `earlier_result`: that file comes from an earlier run and is being redone. It is not stale yet; the reports call it stale only if the step does not end `ok`.

While the run goes, a step whose tasks have all finished shows `done`. If a later process of the same step then starts, the step shows `running` again. Times are in UTC; a step's start and end come from its tasks.

## What it does not do yet

It does not serve a page and does not estimate the hours left. It is the data a live view of the run can be built on.
