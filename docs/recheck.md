# Newer Databases Since a Run

## What This Does

A finished sample keeps the results of the databases it ran with. ClinVar publishes a new file every month, and Ensembl, PCGR, pypgx, PharmCAT, IPD-IMGT/HLA and the PGS Catalog release less often. `bin/owngenome_recheck.py` says which of them has a newer release than the one your run used, which steps that would change, and what to do about it. With `--rerun` it then runs the same `run-all.sh` command again, so Nextflow redoes only those steps and takes the rest from its cache.

```bash
export GENOME_DIR=/path/to/your/data
python3 bin/owngenome_recheck.py --sample your_sample
```

Standard library only, so the plain `python3` on the host is enough. It needs the sample's `run_manifest.tsv`, which `run-all.sh` writes (`bin/write_manifest.sh`).

## What it checks

| Item | What the run used (`run_manifest.tsv`) | Where the newer release comes from | Steps it feeds |
|---|---|---|---|
| ClinVar | `clinvar_file_date`, the `##fileDate` of the ClinVar file | The `Last-Modified` date of NCBI's `clinvar.vcf.gz` | 06, 13 (VEP adds its significance), 20 (the chrM table only), then 23, 27 (the contraindication check), 30, 31, 38, 39 |
| Ensembl / VEP cache | `VEP_CACHE_RELEASE` | Ensembl's REST API | 13, then 23, 30, 31, 38 |
| PCGR / CPSR data | `PCGR_DATA_BUNDLE` | PCGR's newest GitHub release | 17 |
| pypgx and its bundle | `PYPGX_BUNDLE_VERSION` | PyPI, the pypgx biocontainer and the pypgx-bundle tags | 32, then 27, 36 |
| PharmCAT | `PHARMCAT_IMAGE` | PharmCAT's newest GitHub release | 07, 27, then 36 |
| IPD-IMGT/HLA | `hla_database` | The `release_version.txt` of the IPD-IMGT/HLA repository | 08 |
| Each PGS score | `pgs:<ID>`, its `HmPOS_date` | The PGS Catalog's REST API, then the `Last-Modified` date of the harmonised GRCh38 file | 25, then 26 |

The pins of steps 50 to 53 (`MEGANE_IMAGE`, `MOCHA_PANEL_VERSION`, `MOCHA_IMAGE`, `VNTYPER_IMAGE`, `LOCITYPER_PANEL_VERSION`, `LOCITYPER_IMAGE`) have no upstream lookup: each is compared only with this checkout's `versions.env`. An image is compared by its full reference, so a new `@sha256` digest at the same tag counts as a change, for these and for `PHARMCAT_IMAGE`. A pin read only by steps this run did not run (an opt-in step not in `TOOLS`, or one skipped for missing data) is not listed: a rerun would not run that step either.

The HTML report (24) and the text report are written again after any of them. The ClinVar release inside the PCGR bundle (`cpsr_clinvar_release`) moves with the PCGR bundle, never with step 06's ClinVar file.

The lookups are the ones the monthly CI check uses (`scripts/ci/freshness.py`). It also reads the newest release of this pipeline (`--owngenome-repo`, default `GeiserX/OwnGenome`) and that release's `versions.env`, to say whether a newer pin is already in a release.

### What is sent

Only release names: the pinned versions above and the PGS IDs your run scored. Nothing about the sample, no genotype, no result. It prints every host the first time it asks it: `ftp.ncbi.nlm.nih.gov`, `rest.ensembl.org` and Ensembl's FTP host, `api.github.com`, `raw.githubusercontent.com`, the Docker Hub, ghcr.io and quay.io registries, `pypi.org`, `www.pgscatalog.org` and `ftp.ebi.ac.uk`. `GITHUB_TOKEN`, when set, goes to `api.github.com` only. On a machine with no network, write the answers into a file and pass `--offline FILE` (the format is in `tests/fixtures/recheck/offline.json`).

## What it says

Each item gets one class:

| Class | Meaning | What to do |
|---|---|---|
| `rerun-ready` | What is installed here is newer than what the run used: a newer ClinVar file or PGS file on disk, or a pin in this checkout (an image re-pinned to a new digest at the same tag too) | Rerun (for a pin, run `./scripts/setup.sh` first so its data is there) |
| `refresh-available` | ClinVar only: NCBI's file is more than 2 days newer than the one installed | `./scripts/setup.sh --refresh clinvar $GENOME_DIR`, then rerun; `--rerun --refresh-clinvar` does both, and refuses while any sample of `$GENOME_DIR` is running, since they share the ClinVar files ([the ClinVar lock](#the-clinvar-lock)) |
| `in-pipeline-release` | A pipeline release pins a newer version | Update the pipeline to that release, run `setup.sh`, then rerun |
| `upstream-only` | Upstream has a newer release that no pipeline release pins | Nothing yet for a pinned item; for a PGS score, download it again ([step 25](25-prs.md)) |
| `no-newer-release` | The lookup answered and nothing is newer | Nothing |
| `unknown` | A pin of steps 50 to 53 that the run's manifest does not name (a run of an older version): nothing was compared | Nothing; the next run records it. It does not make the command exit 1 |
| `error` | The lookup failed | Check again later. A failed lookup is never read as "up to date", and the command exits 1 |

The 2-day margin for ClinVar is there because NCBI dates the file a day or two after the release in its `##fileDate`. A newer PGS file date means the catalog wrote the harmonised file again; it does not prove the weights changed.

## Which steps a change redoes

`assets/steps.tsv` has a `reads` column: the manifest keys of the databases each step reads itself. A step whose key moved is marked, then every step downstream of it through the `upstream` column, then the reports. A key with `@self` (step 20's ClinVar) changes only a table of that step's own, so nothing below it is marked.

The Nextflow workflow checks this on every pull request: it runs the all-tools stub again with `-resume` and a ClinVar file of another date, and the steps with a task that ran again must be exactly the steps `bin/owngenome_recheck.py --affected clinvar_file_date` names.

## Rerunning

```bash
python3 bin/owngenome_recheck.py --sample your_sample --rerun
python3 bin/owngenome_recheck.py --sample your_sample --rerun --refresh-clinvar   # refresh ClinVar first
```

`run-all.sh` writes its own command into `logs/run_status.tsv` (`meta command`: a `cd` to the directory it was run from, the switches it reads that were set, then its arguments, each shell-quoted; `meta switches`: the names of those switches). `--rerun` runs that command again with those switches only, so a `TOOLS` set in your shell today does not change the run. `NXF_OPTS` and `NXF_JVM_ARGS` are not recorded, because they can hold a proxy password; the replay takes them from your shell. Nextflow's `-resume` then redoes only the tasks whose inputs changed.

It refuses while a run of the sample is going on (`run_status.tsv` has no `finished_utc`), and while another `--rerun` of the sample holds `recheck/rerun.lock`. A run that was killed without its exit trap leaves no `finished_utc` either: run `run-all.sh` by hand once, and the next check can rerun again.

When the run ends, it compares the steps that ran a task again with the steps the check named, and records the difference. Then it checks again against the run manifest the rerun wrote, with the answers it already has, so `state.json` no longer calls ready what the rerun has just used. The script steps (04b, 14, 29, 03a to 03d) skip when their output exists, so a named one is listed as not redone: delete its output and run its script. A sample that started from FASTQ redoes more once: the first `run-all.sh` call after a run that ended with every task done starts from the BAM and VCF that run published ([the delta report](delta-report.md)), so every step after the variant calls runs again and is listed as redone but not named. Later reruns redo only what the check names.

A sample last run by an older `run-all.sh` has no `meta command` row: run the same command once by hand.

### The ClinVar lock

The ClinVar files are shared by every sample of `$GENOME_DIR`. `--refresh-clinvar` holds `$GENOME_DIR/clinvar/.refresh.lock` exclusive while `setup.sh --refresh clinvar` replaces them, and lets go before the rerun starts. `run-all.sh` holds the same lock shared from its start to its end, and `scripts/06-clinvar-screen.sh` while it runs: a run that starts during a refresh waits for it ("Waiting for the ClinVar refresh ..."), and a refresh refuses (exit 3) while a run holds the lock. It also refuses while any sample's `logs/run_status.tsv` says a run is going on, which covers a run of an older checkout. The shell side needs `flock` (part of util-linux on Linux; on macOS `brew install flock`); without it nothing waits. A bare `nextflow run`, and `setup.sh --refresh clinvar` run by hand, take no lock.

### A bare Nextflow run

A bare `nextflow run main.nf` writes no run manifest, so there is nothing to check. To get the command back, `--launch-dir` prints the last command of that launch directory (its `.nextflow/history`) with `-resume` added. It never runs it.

```bash
python3 bin/owngenome_recheck.py --launch-dir /path/to/launch/dir
```

## The files it writes

- `<sample>/recheck/state.json`, replaced whole on each check: `checked_utc`, the hosts asked, the pipeline release it read, one entry per item (`item`, `label`, `ran_with`, `installed`, `upstream`, `class`, `steps`, `action`, and `error` or `in_owngenome_release` when they apply; an image item also has `ran_with_ref` and `installed_ref`, the full references with any digest), `affected_steps` (the steps a rerun would redo for the `rerun-ready` and `refresh-available` items), and after a rerun a `rerun` block (`started_utc`, `finished_utc`, `rc`, `command`, `predicted_vs_redone`). The last rerun stays until another one replaces it.
- `<sample>/recheck/log.tsv`, only ever appended: one row per item and check, and a row as each rerun starts and ends.
- `<sample>/recheck/rerun.lock`, an empty file each `--rerun` locks while it runs.
- `$GENOME_DIR/clinvar/.refresh.lock`, an empty file: [the ClinVar lock](#the-clinvar-lock). `run-all.sh` and step 06 make it too when it is missing.

## In the reports and on the live site

The next reports read `recheck/state.json`: both show one "Database re-check" line in the run facts. A check older than the sample's `run_manifest.tsv` describes an earlier run, so the reports leave it out until the next check. When the sample has a previous run in `history/` ([the delta report](delta-report.md)), a section whose results changed since it carries a badge such as "New since 2026-10-01T09:15:00Z", in both reports, and the delta report heads its changes the same way. The live site's re-check plugin (`bin/owngenome_live_recheck.py`, `assets/live/plugins/recheck.js`) shows the same line as a banner and marks each named step on the step map with a dot whose tooltip says "newer ClinVar available" or "redone 2026-10-05".

The Nextflow `HTML_REPORT` process sees no `history/` folder, so its report has no "New since" badge.

## Checking on a schedule

Nothing here installs itself. To check once a month, copy one of these. Each runs the check alone; uncomment `--rerun` only if you want the rerun to start by itself, which can take hours and the whole machine.

crontab (`crontab -e`), at 06:00 on the 5th of each month:

```text
0 6 5 * * GENOME_DIR=/path/to/your/data /usr/bin/python3 /path/to/OwnGenome/bin/owngenome_recheck.py --sample your_sample >> /path/to/your/data/your_sample/recheck-cron.log 2>&1
# 0 6 5 * * GENOME_DIR=/path/to/your/data /usr/bin/python3 /path/to/OwnGenome/bin/owngenome_recheck.py --sample your_sample --rerun >> /path/to/your/data/your_sample/recheck-cron.log 2>&1
```

systemd (Linux, as your user): `~/.config/systemd/user/owngenome-recheck.service`

```ini
[Unit]
Description=Check a genome against newer databases

[Service]
Type=oneshot
Environment=GENOME_DIR=/path/to/your/data
ExecStart=/usr/bin/python3 /path/to/OwnGenome/bin/owngenome_recheck.py --sample your_sample
# ExecStart=/usr/bin/python3 /path/to/OwnGenome/bin/owngenome_recheck.py --sample your_sample --rerun
```

and `~/.config/systemd/user/owngenome-recheck.timer`, then `systemctl --user enable --now owngenome-recheck.timer`:

```ini
[Unit]
Description=Monthly database re-check

[Timer]
OnCalendar=monthly
Persistent=true

[Install]
WantedBy=timers.target
```

`Persistent=true` runs a check missed while the machine was off at the next start.

launchd (macOS): `~/Library/LaunchAgents/org.owngenome.recheck.plist`, then `launchctl load ~/Library/LaunchAgents/org.owngenome.recheck.plist`:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>org.owngenome.recheck</string>
  <key>ProgramArguments</key>
  <array>
    <string>/usr/bin/python3</string>
    <string>/path/to/OwnGenome/bin/owngenome_recheck.py</string>
    <string>--sample</string><string>your_sample</string>
    <!-- <string>--rerun</string> -->
  </array>
  <key>EnvironmentVariables</key>
  <dict><key>GENOME_DIR</key><string>/path/to/your/data</string></dict>
  <key>StartCalendarInterval</key>
  <dict><key>Day</key><integer>5</integer><key>Hour</key><integer>6</integer></dict>
  <key>StandardOutPath</key><string>/path/to/your/data/your_sample/recheck-launchd.log</string>
  <key>StandardErrorPath</key><string>/path/to/your/data/your_sample/recheck-launchd.log</string>
</dict>
</plist>
```

launchd runs an interval missed during sleep when the Mac wakes; one missed while it was off is skipped.

## Tests

`tests/test_owngenome_recheck.py` checks the classes, the steps each item names, a failed lookup, the files it writes, the refusal during a run, the comparison after a rerun and the live-site plugin, each against a planted mistake too. `tests/fake-docker/run-all-recheck.sh` checks that `--rerun` replays the recorded command exactly, with a control where one dropped argument makes it differ.
