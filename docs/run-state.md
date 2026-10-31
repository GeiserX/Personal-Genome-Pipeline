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

In `process_status.tsv`, a process with no task had not started, so its step is `not finished`. The one exception is `HAPLOCHECK` when step 20 did not run or failed: it never gets a task then. `run-all.sh` reads the file the same way.

A `waiting` or `running` step whose report file is older than the run start has `earlier_result`: that file comes from an earlier run and is being redone. It is not stale yet; the reports call it stale only if the step does not end `ok`.

Only tasks tagged with the sample count. A task shared by every sample has no tag, so it never moves its step. Step 28 (MultiQC) has only such a task, so it shows `waiting` until its result row. The same holds for the shared preparation processes (`T1K_BUILD`, `KIR_BUILD`, `Y_POSITIONS`, `LPA_REGIONS`, `BLOODTYPE_REGIONS`, `CHIP_PREPARE`) and for `SOMALIER_RELATE`: they never count toward their steps. `MINIMAP2_INDEX` and `PRS_PREPARE` run once per run under the reference or score-set name, so they never count either, even when that name equals the sample id.

While the run goes, a step whose tasks have all finished shows `done`. If a later process of the same step then starts, the step shows `running` again. Times are in UTC; a step's start and end come from its tasks.

## How long is left

`bin/run_eta.py` says how long each step usually takes on this box and how many hours a run has left. Both are ranges, and each step's range says where it comes from.

```bash
python3 bin/run_eta.py "$GENOME_DIR"                    # every step's usual time on this box
python3 bin/run_eta.py "$GENOME_DIR" --sample <sample>  # and the hours left of that sample's run
```

The first argument is the run's `--outdir`, which is `GENOME_DIR` with `run-all.sh`. It prints JSON; `--out FILE` writes it to a file instead. Times are in seconds. Standard library only.

**The history.** Each call first reads `GENOME_DIR/pipeline_info/trace_*.txt` into `GENOME_DIR/pipeline_info/step_history/<box>.tsv`. It keeps one row per task that ended `COMPLETED` with exit 0: the trace file's name, the task hash and attempt, its step and process, whether it runs once per run, the CPUs and memory it asked for, its `peak_rss`, its wall time, how long its whole step took for that sample, when it ended, and the box's cores. A cached task, a failed attempt and a process of no step never enter it. It holds no sample id and no path. Reading the same trace again adds nothing, and the rows of a trace you later delete stay. The box is `PGP_BOX` if you set it, else the short host name, reduced to letters, digits, `_` and `-` (`--box` overrides both). Run it on the box that ran the pipeline: the history is named after the box that reads the trace.

Only Nextflow writes a trace, so the history fills from `run-all.sh` and from a bare `nextflow run`, both through `nextflow.config`. The single step scripts write no trace and get no history of their own. A trace from before `nextflow.config` gained the extra columns is still read: its task start and end come from `submit`, `duration` and `realtime`.

**A step's usual time.** On a box with history, it is the shortest and the longest time the step took in the last 5 runs that ran it on this box, with `n` the number of runs (`source: box`). A step's time is its first task start to its last task end for one sample. A run does not count for a step when part of the step came from the cache, or when a task of the step failed or was stopped in that run (a failed attempt that a later attempt replaced is fine). With no history, it is the step's hours from `assets/steps.tsv`, the upper end of the [runtime table](hardware-requirements.md#runtime-per-step) measured on 16 threads (`source: docs`): the low end is that figure and the high end is that figure times 16 / cores on a box with fewer than 16 cores. A box with more cores is not scaled below it, so on 16 cores or more the docs give one figure for both ends. A step the table does not measure gets no number (`source: none`): that is a step with no time in the table, or with only the step page's estimate.

**The hours left.** It counts the steps the state file shows `waiting` or `running`. A running step counts only what is left of it after the time it has run.

- The low end is the longest chain of those steps, following `upstream` in `assets/steps.tsv`, at their low ends.
- The high end is the longer of the same chain at the high ends and a simple CPU-work bound: each step's high end times the CPUs its tasks asked for (one CPU for a step with no history), summed, divided by the box's cores. The bound does not model which tasks run side by side; it only stops the estimate from assuming more cores than the box has.
- A running step already past its high end makes the high end `null`, and `note` says `step N is taking longer than usual`.
- A step with no number is listed in `unknown_steps` and adds nothing, so the range leaves it out instead of counting it as 0.

With `--sample`, the trace files of that run, while it is still going, are left out of the step times, so a half-done step does not count as a short one.

## Memory and disk headroom

`bin/run_headroom.py` says how close each task came to the memory it asked for, and how full the disks of the work folder and the results folder are. It works the same on any Linux, macOS or WSL2 machine and needs no host path set in advance.

```bash
python3 bin/run_headroom.py --work-dir "$GENOME_DIR/<sample>/nextflow/work" --outdir "$GENOME_DIR" \
  --trace "$GENOME_DIR"/pipeline_info/trace_<time>.txt
python3 bin/run_headroom.py --events events.jsonl   # folders and max_memory from the weblog metadata
```

It prints plain lines and exits 0. Standard library only.

**Memory per task.** The tasks come from trace files (`--trace`, repeatable) or from Nextflow's weblog messages (`--events`, the `process_*` messages; the latest message of a task wins). For each task, from its last attempt, the ratio is its `peak_rss` over the memory it asked for. At or above the warn percent the line carries a hint. It is a hint, not a fault: `peak_rss` adds up the resident memory of the task's processes, shared and file-mapped pages included, so it can be above what the container was charged ([measured runtimes](hardware-requirements.md#runtime-per-step) show CPSR above its limit without being killed). A task with no request or no `peak_rss` gets no ratio. A task still running shows its request only. An earlier attempt that ended with one of the exits `conf/base.config` retries (143, 137, 104, 134, 139, 247, which usually mean killed for memory or time) counts as one retry.

**Memory of the machine.** Its RAM, as the system running the command sees it. When the pipeline's `max_memory` (`--max-memory`, else the weblog metadata) is above that RAM, it warns: a task can then ask for more memory than there is. Inside WSL2 the RAM is the VM's share. On macOS with Docker Desktop the containers get only the memory set in Docker Desktop, which can be less than the RAM it reports.

**Disks.** The work folder (`--work-dir`, else `workDir` from the weblog metadata) and the results folder (`--outdir`, else `outdir`). For each it gives the used share, the free and total bytes, and the inode use when the filesystem reports inodes. The used share is used over used plus available, as `df` shows it. A folder not created yet is measured on its nearest existing parent and says so. Two folders on one filesystem give one line, `work and results`. The labels are `work` and `results`, never a path.

**Warn percent.** 90 for both by default. `--disk-warn-percent` or `PGP_DISK_WARN_PERCENT` sets the disk one, `--mem-warn-percent` or `PGP_MEM_WARN_PERCENT` the memory one, from 1 to 100.

**As a library.** `snapshot(events, work_dir, outdir, disk_warn_percent, mem_warn_percent, max_memory, trace_rows)` returns `{disk, tasks, host, warn_used_ratio, mem_warn_ratio}`; `task_memory`, `host_memory` and `disk` give each part alone. It holds no path, and it formats nothing: a page or a metrics endpoint formats the numbers itself. It measures when called; it does not poll and does not project when a disk will fill.

## What it does not do yet

It does not serve a page. It is the data a live view of the run can be built on.
