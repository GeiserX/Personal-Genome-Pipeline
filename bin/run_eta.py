#!/usr/bin/env python3
"""run_eta.py: how long each step usually takes on this box, and how many hours a run has left.

  run_eta.py OUTDIR [--sample S] [--cores N] [--box NAME] [--out FILE]

OUTDIR is the --outdir of the Nextflow run (GENOME_DIR with run-all.sh). Each
call first reads the trace files OUTDIR/pipeline_info/trace_*.txt into this
box's history, OUTDIR/pipeline_info/step_history/<box>.tsv, then prints the
estimate as JSON (or writes it to --out).

The history keeps one row per task that ended COMPLETED with exit 0. A task
taken from the cache, a failed task and a task of an unknown process never
enter it. It holds no sample id and no path: the trace file's name, the task
hash, its step and process, what it asked for, what it used and how long its
step took on this box. The box is PGP_BOX, or the short host name, reduced to
letters, digits, '_' and '-'.

Per step, the range is the shortest and the longest time the step took in the
last 5 runs on this box that hold it, with n the number of runs. A step with no
history on this box takes the hours of assets/steps.tsv (the upper end of the
runtime table in docs/hardware-requirements.md, measured on 16 threads): the low
end is that figure, the high end is that figure times 16 / cores when the box has
fewer than 16 cores. A step that table does not measure (no time, or only the step
page's estimate) gets no number. A run in which a task of the step failed or was
stopped, or took part of the step from the cache, does not count for that step.

With --sample, the hours left of that sample's run (bin/run_state.py gives the
state): the low end is the longest chain of waiting and running steps at their
low ends; the high end is the longer of that chain at the high ends and a simple
CPU-work bound (each step's high end times the CPUs it asked for, summed, divided
by the cores). A running step already past its high end makes the high end
unknown. A step with no number is listed under unknown_steps and adds nothing.
Times are seconds. Standard library only.
"""
import argparse
import fcntl
import glob
import json
import os
import re
import socket
import sys
import tempfile
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_state  # noqa: E402

# Nextflow's default trace columns, in its order, then the ones the estimate needs.
# nextflow.config's trace.fields is the two lists joined; tests/test_run_eta.py checks it.
DEFAULT_FIELDS = ["task_id", "hash", "native_id", "name", "status", "exit", "submit", "duration",
                  "realtime", "%cpu", "peak_rss", "peak_vmem", "rchar", "wchar"]
ADDED_FIELDS = ["process", "tag", "attempt", "cpus", "memory", "start", "complete"]
HISTORY_COLUMNS = ["trace_file", "hash", "attempt", "step", "process", "run_level", "cpus", "memory_bytes",
                   "peak_rss_bytes", "realtime_s", "step_span_s", "finished_utc", "box_cpus"]
LAST_RUNS = 5
DOCS_THREADS = 16
MEM_UNITS = {"B": 1, "KB": 1024, "MB": 1024 ** 2, "GB": 1024 ** 3, "TB": 1024 ** 4, "PB": 1024 ** 5}
NUMBER = re.compile(r"^\d+(?:\.\d+)?$")


# --- units: human ('1h 2m 3s', '18.3 GB') and raw (milliseconds, bytes) -------------------

def blank(v):
    return v is None or str(v).strip() in ("", "-")


def parse_duration(v):
    """Seconds from a trace duration: '1h 2m 3s', '850ms', or raw milliseconds."""
    if blank(v):
        return None
    v = str(v).strip()
    if NUMBER.match(v):
        return float(v) / 1000.0
    return run_state.parse_duration(v)


def parse_memory(v):
    """Bytes from a trace memory value: '18.3 GB', '40 MB', or raw bytes."""
    if blank(v):
        return None
    v = str(v).strip()
    if NUMBER.match(v):
        return int(float(v))
    m = re.match(r"^([\d.]+)\s*([KMGTP]?B)$", v, re.I)
    if not m:
        return None
    return int(float(m.group(1)) * MEM_UNITS[m.group(2).upper()])


def parse_time(v):
    """Epoch seconds from a trace time: 'YYYY-MM-DD HH:MM:SS.mmm' (the host's local time) or raw ms."""
    if blank(v):
        return None
    v = str(v).strip()
    if NUMBER.match(v):
        return float(v) / 1000.0 if float(v) > 0 else None
    return run_state.parse_trace_time(v)


def iso(epoch):
    return run_state.iso(epoch)


def epoch_of(v):
    try:
        return datetime.strptime(v, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()
    except (TypeError, ValueError):
        return None


# --- the trace -------------------------------------------------------------------------------

def read_trace(path):
    """Every row of a trace file, old (default columns) or new (with the added ones), human or
    raw units, as task dicts. Old traces have no start or complete: the task ended at submit +
    duration and ran for its realtime before that."""
    tasks = []
    try:
        f = open(path, errors="replace")
    except OSError:
        return tasks
    with f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            proc, tag = run_state.task_name(r.get("name"))
            if not blank(r.get("process")):
                proc = r["process"].rsplit(":", 1)[-1]
            if "tag" in r:
                tag = None if blank(r["tag"]) else r["tag"]
            submit = parse_time(r.get("submit"))
            dur = parse_duration(r.get("duration"))
            real = parse_duration(r.get("realtime"))
            complete = parse_time(r.get("complete"))
            if complete is None and submit is not None and dur is not None:
                complete = submit + dur
            start = parse_time(r.get("start"))
            if start is None and complete is not None:
                start = complete - real if real is not None else submit
            try:
                cpus = int(r["cpus"]) if not blank(r.get("cpus")) else None
            except ValueError:
                cpus = None
            tasks.append({
                "trace_file": os.path.basename(path), "hash": (r.get("hash") or "").replace("/", ""),
                "attempt": r.get("attempt") if not blank(r.get("attempt")) else "1",
                "process": proc, "tag": tag, "status": (r.get("status") or "").upper(),
                "exit": (r.get("exit") or "").strip(), "cpus": cpus,
                "memory_bytes": parse_memory(r.get("memory")), "peak_rss_bytes": parse_memory(r.get("peak_rss")),
                "realtime_s": real, "start": start, "complete": complete})
    return tasks


def process_steps(steps):
    """{process: step} from the catalogue; a process listed by two steps belongs to the first."""
    out = {}
    for s in steps:
        for p in s["processes"]:
            out.setdefault(p, s["step"])
    return out


def run_level(t):
    return t["process"] in run_state.SHARED_PROCESSES or t["tag"] is None


def finished(ts):
    """True when every task of a step's group ended well: per process and tag, the attempt that
    ended last is COMPLETED with exit 0. A failed attempt that a later one replaced is fine."""
    last = {}
    for i, t in enumerate(ts):
        k = (t["process"], t["tag"])
        order = (t["complete"] if t["complete"] is not None else float("-inf"), i)
        if k not in last or order > last[k][0]:
            last[k] = (order, t)
    return all(t["status"] == "COMPLETED" and t["exit"] == "0" for _, t in last.values())


def history_rows(tasks, proc_step, box_cpus):
    """The history rows of one trace. A step's span is its first start to its last end among the
    tasks of one sample (or of the run-level tasks); a span with a cached task in it, or with a
    task whose last attempt failed or was stopped, is not the step's full time and is left empty."""
    groups = {}
    for t in tasks:
        step = proc_step.get(t["process"])
        if step is None:
            continue
        key = (step, "" if run_level(t) else t["tag"].split(":", 1)[0])
        groups.setdefault(key, []).append(t)
    span_of = {}
    for key, ts in groups.items():
        if any(t["status"] == "CACHED" for t in ts) or not finished(ts):
            span_of[key] = None
            continue
        starts = [t["start"] for t in ts if t["start"] is not None]
        ends = [t["complete"] for t in ts if t["complete"] is not None]
        span_of[key] = max(ends) - min(starts) if starts and ends else None
    rows, seen = [], set()
    for (step, g), ts in groups.items():
        for t in ts:
            if t["status"] != "COMPLETED" or t["exit"] != "0" or not t["hash"]:
                continue
            k = (t["hash"], t["attempt"])
            if k in seen:
                continue  # the trace can list a finished task twice
            seen.add(k)
            rows.append({
                "trace_file": t["trace_file"], "hash": t["hash"], "attempt": t["attempt"], "step": step,
                "process": t["process"], "run_level": "yes" if g == "" else "no", "cpus": t["cpus"],
                "memory_bytes": t["memory_bytes"], "peak_rss_bytes": t["peak_rss_bytes"],
                "realtime_s": t["realtime_s"], "step_span_s": span_of[(step, g)],
                "finished_utc": iso(t["complete"]), "box_cpus": box_cpus})
    return rows


# --- the history -----------------------------------------------------------------------------

def box_name(box=None):
    raw = box or os.environ.get("PGP_BOX") or socket.gethostname().split(".", 1)[0]
    return re.sub(r"[^A-Za-z0-9_-]", "_", raw) or "box"


def history_path(outdir, box=None):
    return os.path.join(outdir, "pipeline_info", "step_history", box_name(box) + ".tsv")


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.3f}".rstrip("0").rstrip(".")
    return str(v)


def read_history(path):
    rows = []
    try:
        f = open(path, errors="replace")
    except OSError:
        return rows
    with f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def write_atomic(path, text):
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), prefix=".run_eta.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        mask = os.umask(0)
        os.umask(mask)
        os.chmod(tmp, 0o666 & ~mask)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


def known_tasks(by_file, name):
    """(hash, attempt) of every task the history holds from a trace other than NAME."""
    return {(r["hash"], r["attempt"]) for f, rs in by_file.items() if f != name for r in rs}


def ingest(outdir, box=None, steps=None, box_cpus=None):
    """Read OUTDIR/pipeline_info/trace_*.txt into this box's history and return the history.
    A trace already read is read again and its rows replaced, so a trace still being written
    ends with its full rows; a task already in the history from another trace is not added
    twice; rows of a trace that was deleted stay. Two calls in a row leave the file as it was."""
    steps = steps if steps is not None else run_state.read_steps()
    box_cpus = box_cpus if box_cpus is not None else os.cpu_count()
    path = history_path(outdir, box)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    proc_step = process_steps(steps)
    with open(path + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        old = read_history(path)
        by_file = {}
        for r in old:
            by_file.setdefault(r["trace_file"], []).append({k: r.get(k, "") for k in HISTORY_COLUMNS})
        for tp in sorted(glob.glob(os.path.join(outdir, "pipeline_info", "trace_*.txt"))):
            name = os.path.basename(tp)
            others = known_tasks(by_file, name)
            new = [{k: fmt(v) for k, v in r.items()} for r in history_rows(read_trace(tp), proc_step, box_cpus)
                   if (r["hash"], r["attempt"]) not in others]
            if new or name in by_file:
                by_file[name] = new
        rows = [r for f in sorted(by_file) for r in by_file[f]]
        text = "\t".join(HISTORY_COLUMNS) + "\n" + "".join(
            "\t".join(r[k] for k in HISTORY_COLUMNS) + "\n" for r in rows)
        try:
            with open(path) as f:
                same = f.read() == text
        except OSError:
            same = False
        if not same:
            write_atomic(path, text)
    return rows


# --- the estimate ----------------------------------------------------------------------------

def num(v):
    try:
        return float(v) if v not in (None, "") else None
    except ValueError:
        return None


def docs_scale(cores):
    """The catalogue's hours were measured on 16 threads: a box with fewer cores scales them up,
    a box with more is never scaled below them."""
    return max(1.0, DOCS_THREADS / float(cores or 1))


def step_ranges(history, steps, cores, exclude=()):
    """{step: {low_s, high_s, n, source, label, cpus}}. source 'box': the last LAST_RUNS runs on
    this box; 'docs': the catalogue's measured hours scaled to the cores; 'none': no number
    (the catalogue has no time for the step, or only an estimate)."""
    runs = {}
    for r in history:
        if r["trace_file"] in exclude:
            continue
        span = num(r.get("step_span_s"))
        run = runs.setdefault(r["step"], {}).setdefault(r["trace_file"], {"spans": [], "end": "", "cpus": 0})
        if span is not None:
            run["spans"].append(span)
        run["end"] = max(run["end"], r.get("finished_utc") or "")
        run["cpus"] = max(run["cpus"], int(num(r.get("cpus")) or 0))
    scale = docs_scale(cores)
    out = {}
    for s in steps:
        mine = sorted((v for v in runs.get(s["step"], {}).values() if v["spans"]), key=lambda v: v["end"])[-LAST_RUNS:]
        if mine:
            spans = [max(v["spans"]) for v in mine]
            out[s["step"]] = {"low_s": min(spans), "high_s": max(spans), "n": len(mine), "source": "box",
                              "label": f"on this box, {len(mine)} run{'s' if len(mine) > 1 else ''}",
                              "cpus": max(v["cpus"] for v in mine) or 1}
        elif s["hours"] is not None and s["hours_from"] == "measured":
            h = s["hours"] * 3600.0
            out[s["step"]] = {"low_s": h, "high_s": h * scale, "n": 0, "source": "docs",
                              "label": "from the docs, measured on 16 threads"
                                       + (f", scaled to {cores} cores" if scale > 1 else ""),
                              "cpus": 1}
        else:
            out[s["step"]] = {"low_s": None, "high_s": None, "n": 0, "source": "none",
                              "label": "not measured", "cpus": None}
    return out


def left_of(lo, hi, gone):
    """What is left of a step that has run GONE seconds: (low, high, late). Past its high end,
    late is True and the step adds nothing to the high end, which is then unknown."""
    lo, hi = max(0.0, lo - gone), hi - gone
    return (lo, 0.0, True) if hi < 0 else (lo, hi, False)


def cpu_work(dur, ranges, left, unknown):
    """CPU-seconds still to come: each step's high end times the CPUs its tasks asked for."""
    return sum(dur[k][1] * ((ranges.get(k) or {}).get("cpus") or 1) for k in left if k not in unknown)


def remaining(state, ranges, now, cores):
    """Hours left of a run: {low_s, high_s (None when a running step is past its high end),
    unknown_steps, late_steps, note, basis}."""
    steps = {s["step"]: s for s in state.get("steps", [])}
    finished = (state.get("run") or {}).get("state") == "finished"
    left = [] if finished else [k for k, s in steps.items() if s["state"] in ("waiting", "running")]
    dur, unknown, late = {}, [], []
    for k in left:
        r = ranges.get(k) or {}
        if r.get("low_s") is None:
            unknown.append(k)
            dur[k] = (0.0, 0.0)
            continue
        lo, hi = r["low_s"], r["high_s"]
        started = epoch_of(steps[k].get("started_utc"))
        if steps[k]["state"] == "running" and started is not None:
            lo, hi, is_late = left_of(lo, hi, max(0.0, now - started))
            if is_late:
                late.append(k)
        dur[k] = (lo, hi)

    def chain(k, i, memo, path=()):
        if k in memo:
            return memo[k]
        if k in path:
            return 0.0
        up = [chain(u, i, memo, path + (k,)) for u in steps[k].get("upstream", []) if u in dur]
        memo[k] = dur[k][i] + max(up, default=0.0)
        return memo[k]
    low = max((chain(k, 0, {}) for k in left), default=0.0)
    high_path = max((chain(k, 1, {}) for k in left), default=0.0)
    work = cpu_work(dur, ranges, left, unknown)
    bound = work / float(cores or 1)
    high = None if late else max(high_path, bound, low)
    note = "; ".join(f"step {k} is taking longer than usual" for k in late) or None
    return {"low_s": low, "high_s": high, "unknown_steps": unknown, "late_steps": late, "note": note,
            "basis": ("low: the longest chain of waiting and running steps at their low ends; high: the longer of "
                      f"that chain at the high ends and their CPU work ({work / 3600:.1f} CPU-hours) over {cores} cores"
                      + ("; steps with no number add nothing" if unknown else ""))}


def current_traces(sample_dir, state):
    """The trace files of a run still going: its own partial spans are not history yet."""
    if (state.get("run") or {}).get("state") != "running":
        return set()
    started = epoch_of((state.get("run") or {}).get("started_utc"))
    return {os.path.basename(p) for p in run_state.trace_files(sample_dir, started)}


def estimate(outdir, sample=None, state=None, now=None, cores=None, box=None, steps=None):
    """{box, cores, steps: {step: range}, remaining: {...} or None}."""
    steps = steps if steps is not None else run_state.read_steps()
    cores = cores or os.cpu_count() or 1
    now = now if now is not None else datetime.now(tz=timezone.utc).timestamp()
    history = ingest(outdir, box, steps)
    sample_dir = os.path.join(outdir, sample) if sample else None
    if sample and state is None:
        state = run_state.collect(sample, sample_dir)
    exclude = current_traces(sample_dir, state) if state and sample_dir else set()
    ranges = step_ranges(history, steps, cores, exclude)
    return {"box": box_name(box), "cores": cores, "generated_utc": iso(now),
            "steps": ranges, "remaining": remaining(state, ranges, now, cores) if state else None}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("outdir", help="the run's --outdir (GENOME_DIR with run-all.sh)")
    ap.add_argument("--sample", help="also the hours left of this sample's run")
    ap.add_argument("--cores", type=int, help="the box's cores (default: all of this machine's)")
    ap.add_argument("--box", help="the history's name for this box (default: PGP_BOX or the short host name)")
    ap.add_argument("--out", help="write the JSON here (default: print it)")
    a = ap.parse_args(argv)
    est = estimate(a.outdir, a.sample, cores=a.cores, box=a.box)
    text = json.dumps(est, indent=1) + "\n"
    if a.out:
        write_atomic(a.out, text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
