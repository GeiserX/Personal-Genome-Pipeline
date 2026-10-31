# Watch a Run Live

A full run takes a day or more. The live site shows it while it goes, in a browser on any device that can reach the machine: which step is running, which finished, which failed or was skipped, and one line on what each finished step found. The page repaints as messages arrive; nobody has to reload it.

`bin/pgp_live.py` is the server. It uses Python's standard library only, so the `python3` on the host is enough, and it has its own password, access log and limits, so it needs no web server in front of it.

## Start it

```bash
export GENOME_DIR=/path/to/your/data
./scripts/pgp-live.sh hash-password                         # once: asks a password, stores only its hash
./scripts/pgp-live.sh                                       # the site, on http://127.0.0.1:8787
PGP_LIVE_URL=http://127.0.0.1:8788 ./scripts/run-all.sh your_name male   # in another terminal
```

`./scripts/pgp-live.sh your_name` does the same and prints the address of that run's page. `hash-password` writes `GENOME_DIR/.pgp-live/password` with mode 600. The file holds a PBKDF2-SHA256 hash, never the password. Open http://127.0.0.1:8787, log in, and pick the run.

With `PGP_LIVE_URL` set, `run-all.sh` writes `<sample>/nextflow/live.config` and passes it to Nextflow with `-c`. That config loads the nf-weblog plugin (pinned in `run-all.sh`; Nextflow downloads it on the first live run) and points it at `PGP_LIVE_URL/ingest/<sample>`, so every task's progress is posted. The script steps that run after the pipeline (GRIDSS, the extra callers, the benchmark) post when they start and how they end, and `run-all.sh` posts once more as it ends, with its exit code. A site that is down or slow costs one warning; the run never waits for it and never fails for it. `PGP_LIVE_URL` must be a plain `http://host:port` URL.

A bare `nextflow run main.nf` can report too. Write a config like this one, with the sample id at the end of the URL (the page counts only the tasks tagged with it), and pass it with `-c`:

```groovy
plugins {
    id 'nf-weblog@1.2.0'
}
weblog {
    enabled = true
    url     = 'http://127.0.0.1:8788/ingest/your_name'
}
```

## What the page shows

- **The step map.** One box per step the run planned, drawn from `assets/steps.tsv`: a step sits right of the steps whose output it reads. The colour is the step's state from `bin/run_state.py` ([Where a Run Stands](run-state.md)): waiting, running, done, failed, skipped, not finished. Skipped steps stay on the map, greyed.
- **A step's detail.** Click a box: its state and reason, start and end, the usual hours, each process with its tasks waiting, running, done, from the cache or failed, and the headline of its sections with a link into the report.
- **Found so far.** One line per finished section, worded by the same code that writes the report, caveat included. A line marked for attention means the report has something to read there; no line ever says a section is fine. The ClinVar line, with no hit, is information, never an all-clear.
- **The report.** Once written, `report.html` opens the run's `<sample>_report.html`, and the files it links to beside it.

`GET /` lists the runs: each folder under `GENOME_DIR` holding `logs/run_status.tsv`, `nextflow/` or `run_manifest.tsv`, and each run the site has messages for.

## Privacy and what is stored

- The site serves only its own pages, the run's state, the report and the files the report links to by a relative link. Any other file is 404, and so is a link that leaves the run folder or a symlink out of `GENOME_DIR`. It never lists a folder.
- Each message is cut to what a progress view needs before it is stored: the task, its process and status, times, CPU and memory figures. Nextflow's parameters, command line and task scripts are dropped. Three paths stay on the server for later views (the output folder, the CPU and memory caps, the work directory); none of them reaches a page.
- Messages go to `<state dir>/runs/<run>/events.jsonl` (state dir: `GENOME_DIR/.pgp-live`, or `PGP_LIVE_STATE_DIR`). The site writes nowhere else, so `GENOME_DIR` can be mounted read-only.
- The access log, `<state dir>/access.log`, has one tab-separated row per page or file read and per login: time, client, the first 8 characters of a hash of the session, action, run, file, status. It never holds a query string, a request body, a result or the password. It moves to `access.log.1` past `--access-log-max-mb` (10).
- Every page carries a Content-Security-Policy that lets it talk to this server only, and the page loads no script from anywhere else. The run's own HTML (the report, CPSR, PharmCAT, MultiQC) is sent sandboxed: no script runs and nothing outside is fetched.

## The password and the limits

| Flag | Default | What it does |
|---|---|---|
| `--password-file FILE` | `PGP_LIVE_PASSWORD_HASH` | The hash to check logins against. Without one, the site listens on loopback only and asks for nothing; any other `--bind` exits 2 |
| `--bind`, `--port` | `127.0.0.1`, `8787` | The browser's address |
| `--ingest-bind`, `--ingest-port` | `127.0.0.1`, `8788` | Where Nextflow and `run-all.sh` post |
| `--session-hours` | 12 | A login lasts at most this long |
| `--session-idle-minutes` | 120 | A login unused this long ends |
| `--login-max-failures` | 5 | Failed logins per client in 15 minutes; then 429 until the window passes |
| `--max-streams` | 32 | Open live pages at once |
| `--max-body-kb` | 1024 | Largest message taken (413 above it) |
| `--allowed-host NAME` | none | A host name the browser may use besides loopback and the bind address (repeatable); any other Host header is 403 |
| `--trust-proxy CIDR` | none | Read `X-Forwarded-For` only from this proxy (repeatable) |
| `--tls-cert`, `--tls-key` | none | Serve HTTPS; the cookie is then `Secure` |
| `--cookie-secure` | off | Mark the cookie `Secure` when TLS ends at a proxy |
| `--config FILE` | none | `key=value` lines setting any flag above (`session_hours=8`), read before the flags |

The session cookie is random (32 bytes), `HttpOnly` and `SameSite=Strict`. A POST from another site is refused.

The ingest port has no password: it takes messages from the processes of the run's own machine. That is why it listens on loopback by default: a local user can post fake progress there, but the port serves nothing, so nobody can read a run through it. A browser's post to it (one with an `Origin` header) is refused.

## Deploy it on a server

The same command runs on any Linux, macOS or WSL2 machine. Pick the reach you need.

**Loopback and ssh.** Leave the defaults and forward the port from your own machine:

```bash
ssh -L 8787:127.0.0.1:8787 user@your-server     # then open http://127.0.0.1:8787
```

**On the local network.** Set the password, then bind the machine's address:

```bash
PGP_LIVE_BIND=0.0.0.0 PGP_LIVE_ALLOWED_HOSTS="genome.example.org" ./scripts/pgp-live.sh
```

`PGP_LIVE_ALLOWED_HOSTS` lists every name or address you will type in the browser: with `0.0.0.0` the site knows none of the machine's addresses, so `http://192.0.2.10:8787` is 403 until `192.0.2.10` is listed. Without TLS the password crosses the network in clear; use this on a network you trust, or add `--tls-cert` and `--tls-key`.

**In Docker.** The image is `PYTHON_IMAGE` from `versions.env`, the data read-only, the state folder writable outside it:

```bash
. ./versions.env
docker run --rm --user "$(id -u):$(id -g)" \
  -p 127.0.0.1:8787:8787 -p 127.0.0.1:8788:8788 \
  -v "$PWD":/pgp:ro -v "$GENOME_DIR":/genome:ro -v "$HOME/.pgp-live":/state \
  "$PYTHON_IMAGE" python3 /pgp/bin/pgp_live.py --root /genome --state-dir /state \
  --bind 0.0.0.0 --ingest-bind 0.0.0.0 --password-file /state/password
```

Inside the container both listeners bind every address, and `-p 127.0.0.1:...` keeps them on the host's loopback. Write the hash first with `python3 bin/pgp_live.py hash-password > ~/.pgp-live/password` and `chmod 600` it.

**Behind a reverse proxy.** Point the proxy at `127.0.0.1:8787`, pass the Host header, and turn off response buffering for `/runs/<run>/events` (the server sends `X-Accel-Buffering: no`, which nginx honours). Then:

```bash
./scripts/pgp-live.sh --allowed-host live.example.org --trust-proxy 127.0.0.1 --cookie-secure
```

**As a service.** A systemd unit for the loopback setup:

```ini
[Unit]
Description=Personal Genome Pipeline live site
After=network.target

[Service]
User=genome
Environment=GENOME_DIR=/path/to/your/data
ExecStart=/path/to/Personal-Genome-Pipeline/scripts/pgp-live.sh
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

## Plugins

The core files (`bin/pgp_live.py`, `assets/live/index.html`, `live.js`, `live.css`, `login.html`) are meant to stay small. A new view adds two files instead of editing them:

- `bin/pgp_live_<name>.py` with `register(app)`. At start the server imports every such file, sorted, and calls it. `app.route(method, pattern, fn, auth=True)` adds a page or an endpoint (behind the login unless `auth=False`), `app.run_file(glob)` lets the browser fetch more files of a run folder, `app.report_filter(fn)` rewrites the served report, `app.state_hook(name, fn)` adds `plugins[name]` to `state.json`, `app.run_kind(fn)` lists more runs, `app.add_argument(...)` adds a flag, `app.log(...)` writes an access-log row, and `app.config` holds the parsed flags.
- `assets/live/plugins/<name>.js`. The run page loads every script in that folder after `live.js`; it calls `window.PGP.register({name, onState, onEvent, panel})`. `onState(state)` gets each new `state.json`, `onEvent(message)` each message from the stream, and `panel(element)` gets a card of its own in the side column.

The page has a replay mode for a demo with no server: `<html data-mode="replay">` reads `frames.json` (`{"interval_ms": 1000, "frames": [{"t": ..., "state": <state.json>}]}`) and plays it with a play button and a slider, through the same drawing code.

`tests/test_pgp_live.py` checks the server, the login, the stored messages, the plugin interface and the page in a headless Chrome; `tests/test_live_post.sh` and `tests/fake-docker/run-all-live.sh` check what `run-all.sh` posts.
