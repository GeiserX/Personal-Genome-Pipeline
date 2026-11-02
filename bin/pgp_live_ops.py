#!/usr/bin/env python3
"""pgp_live_ops.py: a plugin of bin/pgp_live.py. How a run is doing, on the page and in Prometheus.

Per run it gives the hours left (bin/run_eta.py), the memory and disk headroom
(bin/run_headroom.py), when the run last showed activity, and five alerts:
failed, stalled, late, memory and disk. Two outputs, one source:

  state.json   plugins.ops, drawn by assets/live/plugins/ops.js in a card
  GET /metrics the same numbers in the Prometheus text format, for the
               example rules in assets/monitoring/pgp-alerts.yml

Flags it adds to pgp_live.py (also as KEY=VALUE lines in --config):

  --disk-warn-percent N   a disk at or above N% used warns (PGP_DISK_WARN_PERCENT, 90)
  --mem-warn-percent N    a task at or above N% of its memory request is a hint (PGP_MEM_WARN_PERCENT, 90)
  --stall-minutes N       a running run with no activity for N minutes, and no running
                          step still inside its usual time, is stalled (60)
  --metrics-token-file F  a file holding the token Prometheus sends as
                          'Authorization: Bearer <token>'. With a password set and no
                          token, /metrics answers only a logged-in browser
  --eta-box NAME          the step history to read (default PGP_BOX or the short host name)
  --eta-cores N           the cores of the box that runs the pipeline (default this machine's)

It writes nothing: the step history is read, never updated (run bin/run_eta.py
GENOME_DIR after a run to add it), so GENOME_DIR can stay read-only. The paths
it measures are the run folder, the run's work folder (nextflow/work, or an
absolute workDir from the weblog) and an absolute outdir from the weblog, as
this machine sees them: run the site where the run writes. No path, task name
or parameter reaches a page or a metric; the labels are run (the run's folder
name) and fs ('work', 'results' or 'work and results'). Standard library only.
"""
import hmac
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_eta  # noqa: E402
import run_headroom  # noqa: E402
import run_state  # noqa: E402

DEFAULT_STALL_MINUTES = 60
CONTENT_TYPE = "text/plain; version=0.0.4; charset=utf-8"
ALERTS = ("failed", "stalled", "late", "memory", "disk")

# name: (type, help). The one list of what /metrics exports; docs/alerts.md lists the same names.
METRICS = {
    "pgp_runs": ("gauge", "Runs the site knows."),
    "pgp_disk_warn_ratio": ("gauge", "The disk warn share (--disk-warn-percent / 100)."),
    "pgp_host_memory_bytes": ("gauge", "RAM of the machine the site runs on."),
    "pgp_metrics_errors": ("gauge", "Runs whose numbers could not be read in this scrape."),
    "pgp_run_running": ("gauge", "1 while the run is going."),
    "pgp_run_failed": ("gauge", "1 when the run ended badly: a failed step, a non-zero exit or a Nextflow error."),
    "pgp_run_failed_steps": ("gauge", "Steps of the run in state failed."),
    "pgp_run_stalled": ("gauge", "1 when the run is going but has shown no activity for --stall-minutes."),
    "pgp_run_seconds_since_activity": ("gauge", "Seconds since the run's last message or file change."),
    "pgp_run_late_steps": ("gauge", "Running steps already past the high end of their usual time."),
    "pgp_run_unknown_steps": ("gauge", "Steps left with no usual time."),
    "pgp_run_seconds_left_low": ("gauge", "Low end of the time left."),
    "pgp_run_seconds_left_high": ("gauge", "High end of the time left (absent when a step is late)."),
    "pgp_run_tasks_near_memory_request": ("gauge", "Finished tasks whose peak_rss reached the memory warn share of their request."),
    "pgp_run_memory_retries": ("gauge", "Task attempts retried after a memory or time kill."),
    "pgp_run_max_memory_bytes": ("gauge", "The run's max_memory."),
    "pgp_run_max_memory_above_ram": ("gauge", "1 when the run's max_memory is above this machine's RAM."),
    "pgp_disk_used_ratio": ("gauge", "Used share of the filesystem."),
    "pgp_disk_inodes_used_ratio": ("gauge", "Used share of the filesystem's inodes."),
    "pgp_disk_free_bytes": ("gauge", "Bytes free on the filesystem."),
    "pgp_disk_size_bytes": ("gauge", "Size of the filesystem."),
}


# --- the numbers of one run -----------------------------------------------------------------

def current_run(events):
    """The messages of the newest run, as bin/pgp_live.py's current_run: after the last
    pgp_run_finished that is not the newest message, and from the last weblog 'started'."""
    start = 0
    for i, ev in enumerate(events[:-1]):
        if ev.get("event") == "pgp_run_finished":
            start = i + 1
    for i in range(len(events) - 1, start - 1, -1):
        if events[i].get("event") == "started":
            start = i
            break
    return events[start:]


def _abs_dir(p):
    return p if isinstance(p, str) and os.path.isabs(p) else None


def _mtime(p):
    try:
        return os.path.getmtime(p)
    except OSError:
        return None


def last_activity(run_dir, events, started=None):
    """Epoch of the newest sign of life: a message's utcTime, the run's run_status.tsv, or a
    trace file of the run."""
    times = [run_eta.epoch_of(ev.get("utcTime")) for ev in events]
    if run_dir and os.path.isdir(run_dir):
        times.append(_mtime(os.path.join(run_dir, "logs", "run_status.tsv")))
        times += [_mtime(p) for p in run_state.trace_files(run_dir, started)]
    times = [t for t in times if t is not None]
    return max(times) if times else None


def folders(run_dir, events):
    """(work, results) to measure: the run folder's nextflow/work, else an absolute workDir from
    the weblog; the run folder, else an absolute outdir from the weblog."""
    m_work, m_out, _ = run_headroom.metadata(events)
    have = run_dir and os.path.isdir(run_dir)
    work = os.path.join(run_dir, "nextflow", "work") if have and os.path.isdir(
        os.path.join(run_dir, "nextflow", "work")) else _abs_dir(m_work)
    return work, (run_dir if have else _abs_dir(m_out))


def trace_rows(run_dir, sample, started):
    """This run's own tasks from its trace files (a run with no weblog messages has only these)."""
    if not run_dir or not os.path.isdir(run_dir):
        return []
    return [r for p in run_state.trace_files(run_dir, started) for r in run_eta.read_trace(p)
            if r.get("tag") is not None and run_state.is_sample(r["tag"], sample)]


def run_ops(run_dir, state, events, cfg, now=None, steps=None, physical=None):
    """{run, hours_left, memory, disk, activity, alerts} for one run: the plugin's state hook,
    and the source of every metric. No path in it."""
    now = now if now is not None else time.time()
    steps = steps if steps is not None else run_state.read_steps(getattr(cfg, "steps", run_state.STEPS_TSV))
    root = getattr(cfg, "root", None) or (os.path.dirname(run_dir) if run_dir else None)
    cur = current_run(list(events))
    sample = state.get("sample") or (os.path.basename(run_dir) if run_dir else "")
    run = state.get("run") or {}
    # run-all.sh posts pgp_run_finished as its last message; bin/pgp_live.py marks the run finished
    # on it too, before calling this hook.
    ended = bool(cur) and cur[-1].get("event") == "pgp_run_finished"
    run_st = "finished" if ended and run.get("state") == "running" else run.get("state")
    running = run_st == "running"
    started = run_eta.epoch_of(run.get("started_utc"))

    # Hours left: run_eta's ranges from the history as it is, without adding to it.
    cores = getattr(cfg, "eta_cores", None) or os.cpu_count() or 1
    history = run_eta.read_history(run_eta.history_path(root, getattr(cfg, "eta_box", None))) if root else []
    exclude = (run_eta.live_traces(root) if root else set()) | \
        (run_eta.current_traces(run_dir, state) if run_dir and os.path.isdir(run_dir) else set())
    ranges = run_eta.step_ranges(history, steps, cores, exclude)
    rem = run_eta.remaining(state, ranges, now, cores) if running else None

    # Headroom: run_headroom's snapshot on the folders found here. The weblog metadata is taken
    # out of the events first, so snapshot() cannot fall back to a relative outdir, which would be
    # read against this server's own folder.
    work, results = folders(run_dir, cur)
    _, _, cap = run_headroom.metadata(cur)
    plain = [{k: v for k, v in ev.items() if k != "metadata"} for ev in cur]
    snap = run_headroom.snapshot(plain, work or "", results or "", getattr(cfg, "disk_warn_percent", None),
                                 getattr(cfg, "mem_warn_percent", None), cap, trace_rows(run_dir, sample, started),
                                 steps, physical)
    snap["disk"] = [d for d in snap["disk"] if d["used_ratio"] is not None or d["free_bytes"] is not None]

    # What went wrong or is going wrong.
    failed_steps = [s["step"] for s in state.get("steps", []) if s["state"] == "failed"]
    exits = [ev.get("exit") for ev in cur if ev.get("event") == "pgp_run_finished"]
    exit_code = exits[-1] if exits else None
    error = any(ev.get("event") == "error" for ev in cur)
    finished = run_st == "finished"
    failed = bool(failed_steps) or (finished and ((exit_code not in (None, 0)) or error))
    act = last_activity(run_dir, cur, started)
    quiet = None if act is None else max(0.0, now - act)
    stall_s = float(getattr(cfg, "stall_minutes", None) or DEFAULT_STALL_MINUTES) * 60
    covered = []  # running steps still inside their usual high end: a long task posts nothing meanwhile
    for s in state.get("steps", []):
        r = ranges.get(s["step"]) or {}
        st = run_eta.epoch_of(s.get("started_utc"))
        if s["state"] == "running" and r.get("high_s") is not None and st is not None and now - st < r["high_s"]:
            covered.append(s["step"])
    stalled = running and quiet is not None and quiet >= stall_s and not covered
    late = list(rem["late_steps"]) if rem else []
    mem_hints = [t for t in snap["tasks"] if t["warn"]]
    retries = sum(t["oom_retry"] for t in snap["tasks"])
    disk_warn = [d for d in snap["disk"] if d["warn"]]

    alerts = []
    if failed:
        why = (f"step{'s' if len(failed_steps) > 1 else ''} {', '.join(failed_steps)} failed" if failed_steps
               else f"the run ended with exit {exit_code}" if exit_code not in (None, 0) else "Nextflow reported an error")
        alerts.append({"kind": "failed", "level": "critical", "text": why})
    if stalled:
        alerts.append({"kind": "stalled", "level": "warning",
                       "text": f"no activity for {quiet / 60:.0f} minutes, and no running step is inside its usual time"})
    if late:
        alerts.append({"kind": "late", "level": "warning",
                       "text": "; ".join(f"step {k} is taking longer than usual" for k in late)})
    if retries or snap["host"]["warn"]:
        parts = []
        if retries:
            parts.append(f"{retries} task attempt{'s' if retries > 1 else ''} retried after a memory or time kill")
        if snap["host"]["warn"]:
            parts.append("max_memory is above this machine's memory")
        alerts.append({"kind": "memory", "level": "warning", "text": "; ".join(parts)})
    elif mem_hints:
        alerts.append({"kind": "memory", "level": "info",
                       "text": f"{len(mem_hints)} task{'s' if len(mem_hints) > 1 else ''} reached "
                               f"{snap['mem_warn_ratio'] * 100:.0f}% of the memory asked for (peak_rss counts "
                               "shared pages, so this may be fine)"})
    for d in disk_warn:
        alerts.append({"kind": "disk", "level": "critical" if (d["used_ratio"] or 0) >= 0.98 else "warning",
                       "text": f"disk {d['fs']} is at or above {snap['warn_used_ratio'] * 100:.0f}% used"})

    return {
        "run": {"state": run_st, "failed": failed, "failed_steps": failed_steps, "exit": exit_code,
                "stalled": stalled, "seconds_since_activity": quiet, "stall_seconds": stall_s},
        "hours_left": None if rem is None else {
            "low_s": rem["low_s"], "high_s": rem["high_s"], "late_steps": late,
            "unknown_steps": rem["unknown_steps"], "basis": rem["basis"], "cores": cores},
        "memory": {"host_bytes": snap["host"]["memory_bytes"], "max_memory_bytes": snap["host"]["max_memory_bytes"],
                   "max_memory_above_ram": snap["host"]["warn"], "warn_ratio": snap["mem_warn_ratio"],
                   "near_request": len(mem_hints), "retries": retries,
                   "tasks": [{k: t[k] for k in ("step", "process", "state", "request_bytes", "peak_rss_bytes",
                                                 "ratio", "warn", "oom_retry")} for t in mem_hints]},
        "disk": [{k: d[k] for k in ("fs", "used_ratio", "free_bytes", "size_bytes", "inodes_used_ratio",
                                     "missing", "warn")} for d in snap["disk"]],
        "disk_warn_ratio": snap["warn_used_ratio"],
        "alerts": alerts,
    }


# --- the Prometheus text --------------------------------------------------------------------

def _label(v):
    return str(v).replace("\\", "\\\\").replace("\n", "\\n").replace('"', '\\"')


def _num(v):
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    return repr(float(v))


def metrics_text(runs, disk_warn_ratio, host_bytes=None, errors=0):
    """The Prometheus text format of RUNS, a list of (run name, run_ops() result). A value that is
    not known is left out, never written as 0."""
    samples = {k: [] for k in METRICS}

    def put(name, value, **labels):
        if value is None:
            return
        lab = ",".join(f'{k}="{_label(v)}"' for k, v in labels.items())
        samples[name].append(f"{name}{{{lab}}} {_num(value)}" if lab else f"{name} {_num(value)}")

    put("pgp_runs", len(runs))
    put("pgp_disk_warn_ratio", disk_warn_ratio)
    put("pgp_host_memory_bytes", host_bytes)
    put("pgp_metrics_errors", errors)
    for name, o in runs:
        r, h, m = o["run"], o["hours_left"], o["memory"]
        put("pgp_run_running", r["state"] == "running", run=name)
        put("pgp_run_failed", r["failed"], run=name)
        put("pgp_run_failed_steps", len(r["failed_steps"]), run=name)
        put("pgp_run_stalled", r["stalled"], run=name)
        put("pgp_run_seconds_since_activity", r["seconds_since_activity"], run=name)
        if h is not None:
            put("pgp_run_late_steps", len(h["late_steps"]), run=name)
            put("pgp_run_unknown_steps", len(h["unknown_steps"]), run=name)
            put("pgp_run_seconds_left_low", h["low_s"], run=name)
            put("pgp_run_seconds_left_high", h["high_s"], run=name)
        put("pgp_run_tasks_near_memory_request", m["near_request"], run=name)
        put("pgp_run_memory_retries", m["retries"], run=name)
        put("pgp_run_max_memory_bytes", m["max_memory_bytes"], run=name)
        if m["max_memory_bytes"] is not None:
            put("pgp_run_max_memory_above_ram", m["max_memory_above_ram"], run=name)
        for d in o["disk"]:
            put("pgp_disk_used_ratio", d["used_ratio"], run=name, fs=d["fs"])
            put("pgp_disk_inodes_used_ratio", d["inodes_used_ratio"], run=name, fs=d["fs"])
            put("pgp_disk_free_bytes", d["free_bytes"], run=name, fs=d["fs"])
            put("pgp_disk_size_bytes", d["size_bytes"], run=name, fs=d["fs"])
    out = []
    for name, (typ, help_) in METRICS.items():
        if samples[name]:
            out += [f"# HELP {name} {help_}", f"# TYPE {name} {typ}"] + samples[name]
    return "\n".join(out) + "\n"


# --- the plugin -------------------------------------------------------------------------------

def read_token(path):
    if not path:
        return None
    with open(path) as f:
        tok = f.read().strip()
    return tok or None


def metrics_allowed(req, token):
    """Without a password the site is loopback only and /metrics is open, like every page. With
    one: a logged-in browser, or the bearer token of --metrics-token-file."""
    live = getattr(req, "live", None)
    if not getattr(live, "password_hash", None):
        return True
    if getattr(req, "sid", None):
        return True
    auth = (req.headers.get("Authorization") or "") if getattr(req, "headers", None) is not None else ""
    scheme, _, given = auth.partition(" ")
    return bool(token) and scheme.lower() == "bearer" and hmac.compare_digest(given.strip().encode(), token.encode())


def register(app):
    app.add_argument("--disk-warn-percent", help="a disk at or above this use warns (default PGP_DISK_WARN_PERCENT or 90)")
    app.add_argument("--mem-warn-percent", help="a task at or above this share of its request is a hint (default PGP_MEM_WARN_PERCENT or 90)")
    app.add_argument("--stall-minutes", type=float, default=DEFAULT_STALL_MINUTES,
                     help="a running run quiet this long, with no step inside its usual time, is stalled")
    app.add_argument("--metrics-token-file", help="the bearer token Prometheus sends to /metrics")
    app.add_argument("--eta-box", help="the step history to read (default PGP_BOX or the short host name)")
    app.add_argument("--eta-cores", type=int, help="the cores of the box that runs the pipeline")
    token = {}

    def hook(run_dir, state, events):
        try:
            return run_ops(run_dir, state, events, app.config)
        except ValueError as e:  # a warn percent outside 1 to 100: say so on the page, keep the page
            return {"error": str(e)[:200]}

    def metrics(req, _m):
        if "t" not in token:
            token["t"] = read_token(getattr(app.config, "metrics_token_file", None))
        if not metrics_allowed(req, token["t"]):
            app.log("metrics", status=401)
            return req.send_body(401, "text/plain; charset=utf-8", "login or bearer token required\n")
        live = req.live
        runs, errors = [], 0
        for name in sorted(live.runs()):
            try:
                st = live.state(live.run(name))
                o = st["plugins"].get("ops")
                if o is None or "error" in o:
                    raise ValueError(o and o["error"])
                runs.append((name, o))
            except Exception as e:  # one unreadable run must not take the scrape down
                errors += 1
                print(f"pgp_live_ops: {name}: {type(e).__name__}", file=sys.stderr)
        try:
            warn = run_headroom.warn_percent(getattr(app.config, "disk_warn_percent", None),
                                             "PGP_DISK_WARN_PERCENT") / 100.0
        except ValueError:
            warn = None
        app.log("metrics", status=200)
        req.send_body(200, CONTENT_TYPE, metrics_text(runs, warn, run_headroom.physical_memory(), errors))

    app.state_hook("ops", hook)
    app.route("GET", r"/metrics", metrics, auth=False)
