# Run Alerts

The live site (`bin/pgp_live.py`, its guide is `docs/live-site.md`) can warn you before a run goes wrong, and say so to Prometheus. The plugin `bin/pgp_live_ops.py` adds a **Run health** card to each run's page and a `GET /metrics` endpoint. Both read the same numbers:

- **Hours left**, as a range, from `bin/run_eta.py` ([How long is left](run-state.md#how-long-is-left)).
- **Memory and disk headroom**, from `bin/run_headroom.py` ([Memory and disk headroom](run-state.md#memory-and-disk-headroom)).
- **When the run last showed activity**: its newest message to the site, or the newest change to its `logs/run_status.tsv` or its trace files.

It works the same on any Linux, macOS or WSL2 machine. It needs nothing beyond the live site: no Prometheus is required to see the card.

## The five alerts

| Alert | When | Level |
|---|---|---|
| failed | a step of the run is in state `failed`, or the run ended with a non-zero exit (`pgp_run_finished`), or Nextflow reported an error | critical |
| stalled | the run is going, nothing happened for `--stall-minutes` (60), and no running step is still inside its usual time | warning |
| late | a running step is already past the high end of its usual time on this box | warning |
| memory | a task was retried after a memory or time kill, or the run's `max_memory` is above this machine's memory (warning); a finished task reached `--mem-warn-percent` (90) of the memory it asked for (info, a hint only: `peak_rss` counts shared pages) | warning or info |
| disk | the work or the results disk, or its inodes, is at or above `--disk-warn-percent` (90) | warning, critical from 98% |

A long task posts nothing between its start and its end, so a quiet run is not stalled while its running step is inside its usual time. A step with no usual time (no history on this box and no measured time in the [runtime table](hardware-requirements.md#runtime-per-step)) cannot cover a quiet run: after `--stall-minutes` with no activity, that run shows as stalled. Raise `--stall-minutes` if your box runs such a step for longer than that.

Late and stalled need the step history of this box. The site reads `GENOME_DIR/pipeline_info/step_history/<box>.tsv` and never writes it, so `GENOME_DIR` can stay read-only. To add a finished run to it, run `python3 bin/run_eta.py "$GENOME_DIR"` once. With no history, each step uses its measured time from the runtime table.

## Flags

The plugin adds these flags to `bin/pgp_live.py`. Each also works as a `KEY=VALUE` line in its `--config` file (`disk_warn_percent=85`).

| Flag | Default | What it does |
|---|---|---|
| `--disk-warn-percent N` | `PGP_DISK_WARN_PERCENT` or 90 | A disk at or above N% used, or N% of its inodes, warns. One knob for the card and the rules |
| `--mem-warn-percent N` | `PGP_MEM_WARN_PERCENT` or 90 | A task at or above N% of its memory request is a hint |
| `--stall-minutes N` | 60 | Minutes with no activity before a running run is stalled |
| `--metrics-token-file FILE` | none | A file holding the token Prometheus sends as `Authorization: Bearer <token>` |
| `--eta-box NAME` | `PGP_BOX` or the short host name | The step history to read |
| `--eta-cores N` | this machine's cores | The cores of the box that runs the pipeline |

The paths it measures are the run's folder, its `nextflow/work` folder, and for a run started with a bare `nextflow run` the absolute `workDir` and `outdir` from its weblog messages. They are measured as this machine sees them, so run the site on the machine that runs the pipeline. A relative `outdir` is not measured.

## Scrape it with Prometheus

Without a password the site listens on loopback only, and `/metrics` is open like every page. With a password, `/metrics` answers a logged-in browser or the bearer token, nothing else. Make a token and keep it out of the repository:

```bash
python3 -c "import secrets; print(secrets.token_urlsafe(32))" > ~/.config/owngenome/metrics.token
chmod 600 ~/.config/owngenome/metrics.token
python3 bin/pgp_live.py --root "${GENOME_DIR:?Set GENOME_DIR to your data directory}" \
  --password-file ~/.config/owngenome/live.hash --metrics-token-file ~/.config/owngenome/metrics.token
```

Then in `prometheus.yml`:

```yaml
scrape_configs:
  - job_name: owngenome
    scrape_interval: 1m
    authorization:
      credentials_file: /etc/prometheus/owngenome.token
    static_configs:
      - targets: ["127.0.0.1:8787"]
rule_files:
  - /etc/prometheus/rules/pgp-alerts.yml
```

Copy `assets/monitoring/pgp-alerts.yml` to the rules folder. It is an example: the thresholds the page also shows come from the site's flags (`pgp_disk_warn_ratio` carries the disk one, and `pgp_run_stalled` already applies `--stall-minutes`), so the rules and the card agree. The `for:` times and the 98% line of `OwnGenomeDiskFull` are yours to tune. `assets/monitoring/pgp-alerts.test.yml` tests every rule; run `promtool test rules assets/monitoring/pgp-alerts.test.yml` after a change. Alertmanager then sends the alerts wherever you route them.

## The metrics

The labels are `run` (the run's folder name, the sample id) and `fs` (`work`, `results`, or `work and results` when both are on one filesystem). Nothing else: no path, no task, no step, so the number of series grows with the number of runs only. A value that is not known is left out, never written as 0. The sample id does reach Prometheus as a label: use a name that says nothing about the person, as for the page itself.

| Metric | Labels | Meaning |
|---|---|---|
| `pgp_runs` | | Runs the site knows |
| `pgp_disk_warn_ratio` | | `--disk-warn-percent` / 100 |
| `pgp_host_memory_bytes` | | RAM of the machine the site runs on |
| `pgp_metrics_errors` | | Runs whose numbers could not be read in this scrape |
| `pgp_run_running` | run | 1 while the run is going |
| `pgp_run_failed` | run | 1 when the run failed (see the table above) |
| `pgp_run_failed_steps` | run | Steps in state `failed` |
| `pgp_run_stalled` | run | 1 when the run is stalled |
| `pgp_run_seconds_since_activity` | run | Seconds since the last activity |
| `pgp_run_late_steps` | run | Running steps past their usual high end |
| `pgp_run_unknown_steps` | run | Steps left with no usual time |
| `pgp_run_seconds_left_low` | run | Low end of the time left |
| `pgp_run_seconds_left_high` | run | High end of the time left; absent while a step is late |
| `pgp_run_tasks_near_memory_request` | run | Finished tasks at or above `--mem-warn-percent` of their request |
| `pgp_run_memory_retries` | run | Task attempts retried after a memory or time kill |
| `pgp_run_max_memory_bytes` | run | The run's `max_memory` |
| `pgp_run_max_memory_above_ram` | run | 1 when `max_memory` is above this machine's RAM |
| `pgp_disk_used_ratio` | run, fs | Used share of the filesystem |
| `pgp_disk_inodes_used_ratio` | run, fs | Used share of its inodes |
| `pgp_disk_free_bytes` | run, fs | Bytes free |
| `pgp_disk_size_bytes` | run, fs | Size of the filesystem |

The hours left, the late steps and the unknown steps are written only while a run is going. A run started with a bare `nextflow run` and the weblog config of the live site guide gives the same series as one started by `run-all.sh`.

## How it is tested

`tests/test_pgp_live_ops.py` runs the plugin on the made-up run of `tests/fixtures/live/ops/` with a fake disk and a fake host, then against planted mistakes that must fail it. The `promtool` job in `.github/workflows/lint.yml` runs `pgp-alerts.test.yml`, then changes one threshold and requires the test to fail.
