#!/usr/bin/env python3
"""run_headroom.py: how close each task came to its memory request, and how full the work and results disks are.

  run_headroom.py [--events FILE] [--trace FILE ...] [--work-dir DIR] [--outdir DIR]
                  [--max-memory SIZE] [--disk-warn-percent N] [--mem-warn-percent N]

The one place that measures memory and disk headroom; anything that shows them
(a page, metrics) calls snapshot() and formats its result. Standard library
only, on any Linux, macOS or WSL2 box.

Tasks come from Nextflow weblog messages (--events, one JSON object per line,
the process_* messages) or trace files (--trace). Per task, the ratio is its
peak_rss over the memory it asked for. A ratio at or above the warn percent is
a hint, not a fault: peak_rss counts shared and file-mapped pages, so it can be
above what the container was charged. A task with no request, or no peak_rss,
has no ratio; a task still running shows its request only. An attempt that
ended with an exit conf/base.config retries (killed for memory or time) and
that a later attempt replaced counts as a retry.

The disks are the work folder and the results folder: --work-dir and --outdir,
else workDir and outdir from the weblog metadata. There is no default path. A
folder not created yet is measured on its nearest existing parent. Two folders
on the same filesystem are reported once. Labels are 'work' and 'results',
never a path.

The warn percent is 90 for both, or PGP_DISK_WARN_PERCENT and
PGP_MEM_WARN_PERCENT, or the options. The command prints plain lines and
exits 0; a warn percent outside 1 to 100 exits 2.
"""
import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_eta  # noqa: E402
import run_state  # noqa: E402

DEFAULT_WARN_PERCENT = 90
# The exits conf/base.config's errorStrategy retries: usually killed for memory or time.
# tests/test_run_headroom.py checks the two lists agree.
MEMORY_EXITS = {"143", "137", "104", "134", "139", "247"}
DONE = {"COMPLETED", "FAILED", "ABORTED", "CACHED"}


def warn_percent(value=None, env=None):
    """The warn percent: VALUE, else the environment variable ENV, else DEFAULT_WARN_PERCENT. 1 to 100."""
    raw = value if value not in (None, "") else (os.environ.get(env) if env else None)
    if raw in (None, ""):
        return float(DEFAULT_WARN_PERCENT)
    try:
        p = float(raw)
    except (TypeError, ValueError):
        raise ValueError(f"warn percent must be a number from 1 to 100, got {raw!r}")
    if not 1 <= p <= 100:
        raise ValueError(f"warn percent must be a number from 1 to 100, got {raw!r}")
    return p


# --- memory per task -----------------------------------------------------------------------

def _row(r, i):
    """A weblog trace object or a run_eta.read_trace row, as one shape."""
    if "memory_bytes" in r or "peak_rss_bytes" in r:  # run_eta.read_trace
        mem, peak, proc, tag = r.get("memory_bytes"), r.get("peak_rss_bytes"), r.get("process"), r.get("tag")
        key = (r.get("hash"), r.get("attempt"))
    else:  # weblog trace object
        mem, peak = run_eta.parse_memory(r.get("memory")), run_eta.parse_memory(r.get("peak_rss"))
        proc, tag = run_state.task_name(r.get("name"))
        if not run_eta.blank(r.get("process")):
            proc = str(r["process"]).rsplit(":", 1)[-1]
        if not run_eta.blank(r.get("tag")):
            tag = str(r["tag"])
        key = r.get("task_id") if r.get("task_id") is not None else (r.get("hash"), r.get("attempt"))
    try:
        attempt = int(r.get("attempt") or 1)
    except (TypeError, ValueError):
        attempt = 1
    ex = r.get("exit")
    return {"key": key, "order": i, "process": proc, "tag": tag, "attempt": attempt,
            "status": str(r.get("status") or "").upper(), "exit": "" if run_eta.blank(ex) else str(ex).strip(),
            "memory": mem, "peak": peak}


def task_memory(rows, warn=None, steps=None):
    """[{step, process, attempt, state, request_bytes, peak_rss_bytes, ratio, warn, oom_retry}], one per
    task (process and tag), from its last attempt. ROWS are weblog trace objects (the latest message
    of a task wins) or run_eta.read_trace rows."""
    warn = warn_percent(warn)
    proc_step = run_eta.process_steps(steps if steps is not None else run_state.read_steps())
    latest = {}
    for i, r in enumerate(rows):
        t = _row(r, i)
        if t["process"]:
            latest[t["key"]] = t
    groups = {}
    for t in latest.values():
        groups.setdefault((t["process"], t["tag"]), []).append(t)
    out = []
    for (proc, _tag), ts in sorted(groups.items(), key=lambda kv: min(t["order"] for t in kv[1])):
        ts.sort(key=lambda t: (t["attempt"], t["order"]))
        last = ts[-1]
        retries = len({t["attempt"] for t in ts if t["attempt"] < last["attempt"] and t["exit"] in MEMORY_EXITS})
        running = last["status"] not in DONE
        ratio = None
        if not running and last["memory"] and last["peak"] is not None:
            ratio = last["peak"] / float(last["memory"])
        out.append({"step": proc_step.get(proc), "process": proc, "attempt": last["attempt"],
                    "state": "running" if running else last["status"].lower(),
                    "request_bytes": last["memory"], "peak_rss_bytes": None if running else last["peak"],
                    "ratio": ratio, "warn": ratio is not None and ratio * 100 >= warn, "oom_retry": retries})
    return out


# --- memory of the host --------------------------------------------------------------------

def physical_memory():
    """This machine's RAM in bytes (inside WSL2 or a VM: what that system sees), or None."""
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def host_memory(max_memory=None, physical=None):
    """{memory_bytes, max_memory_bytes, warn}: warn when the pipeline's max_memory is above the RAM."""
    phys = physical if physical is not None else physical_memory()
    cap = run_eta.parse_memory(max_memory) if max_memory is not None else None
    return {"memory_bytes": phys, "max_memory_bytes": cap,
            "warn": phys is not None and cap is not None and cap > phys}


# --- disks ---------------------------------------------------------------------------------

def _statvfs(path):
    return os.statvfs(path)


def _dev(path):
    return os.stat(path).st_dev


def nearest_existing(path):
    """PATH, or its nearest parent that exists."""
    p = os.path.abspath(path)
    while not os.path.exists(p):
        up = os.path.dirname(p)
        if up == p:
            break
        p = up
    return p


def _measure(path):
    """(used_ratio, free_bytes, size_bytes, inodes_used_ratio) of the filesystem holding PATH."""
    try:
        s = _statvfs(path)
        frsize = s.f_frsize or s.f_bsize
        used = (s.f_blocks - s.f_bfree) * frsize
        free = s.f_bavail * frsize
        size = s.f_blocks * frsize
        inodes = (s.f_files - s.f_ffree) / float(s.f_files) if s.f_files > 0 else None
    except (AttributeError, OSError):
        u = shutil.disk_usage(path)
        used, free, size, inodes = u.used, u.free, u.total, None
    return (used / float(used + free) if used + free > 0 else None), free, size, inodes


def disk(paths, warn=None):
    """[{fs, used_ratio, free_bytes, size_bytes, inodes_used_ratio, missing, warn}] for PATHS, a list of
    (label, path) pairs; a pair with no path is left out. Folders on one filesystem give one entry,
    whose fs names them all ('work and results')."""
    warn = warn_percent(warn)
    out, by_dev = [], {}
    for label, path in paths:
        if not path:
            continue
        at = nearest_existing(path)
        try:
            dev = _dev(at)
        except OSError:
            continue
        if dev in by_dev:
            e = by_dev[dev]
            e["fs"] += " and " + label
            e["missing"] = e["missing"] or at != os.path.abspath(path)
            continue
        try:
            used, free, size, inodes = _measure(at)
        except OSError:
            continue
        e = {"fs": label, "used_ratio": used, "free_bytes": free, "size_bytes": size,
             "inodes_used_ratio": inodes, "missing": at != os.path.abspath(path),
             "warn": any(r is not None and r * 100 >= warn for r in (used, inodes))}
        by_dev[dev] = e
        out.append(e)
    return out


# --- one snapshot ----------------------------------------------------------------------------

def metadata(events):
    """(workDir, outdir, max_memory) from the last weblog message that carries metadata."""
    work = outdir = cap = None
    for ev in events:
        md = ev.get("metadata") if isinstance(ev, dict) else None
        if not isinstance(md, dict):
            continue
        wf, params = md.get("workflow") or {}, md.get("parameters") or {}
        work = wf.get("workDir") or work
        outdir = params.get("outdir") or outdir
        cap = params.get("max_memory") or cap
    return work, outdir, cap


def event_rows(events):
    return [ev["trace"] for ev in events
            if isinstance(ev, dict) and str(ev.get("event", "")).startswith("process_") and isinstance(ev.get("trace"), dict)]


def snapshot(events=(), work_dir=None, outdir=None, disk_warn_percent=None, mem_warn_percent=None,
             max_memory=None, trace_rows=(), steps=None, physical=None):
    """{disk, tasks, host, warn_used_ratio, mem_warn_ratio}. Folders and max_memory not given come
    from the weblog metadata in EVENTS; tasks come from TRACE_ROWS then EVENTS (the later wins)."""
    dw = warn_percent(disk_warn_percent, "PGP_DISK_WARN_PERCENT")
    mw = warn_percent(mem_warn_percent, "PGP_MEM_WARN_PERCENT")
    events = list(events)
    m_work, m_out, m_cap = metadata(events)
    return {"disk": disk([("work", work_dir or m_work), ("results", outdir or m_out)], dw),
            "tasks": task_memory(list(trace_rows) + event_rows(events), mw, steps),
            "host": host_memory(max_memory or m_cap, physical),
            "warn_used_ratio": dw / 100.0, "mem_warn_ratio": mw / 100.0}


# --- the command -----------------------------------------------------------------------------

def size(b):
    if b is None:
        return "unknown"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if b < 1024 or unit == "TB":
            return f"{b:.0f} {unit}" if unit == "B" else f"{b:.1f} {unit}"
        b /= 1024.0


def pct(r):
    return f"{r * 100:.0f}%"


def lines(snap):
    out = []
    for d in snap["disk"]:
        s = (f"disk {d['fs']}: {pct(d['used_ratio']) if d['used_ratio'] is not None else 'unknown'} used, "
             f"{size(d['free_bytes'])} free of {size(d['size_bytes'])}")
        if d["inodes_used_ratio"] is not None:
            s += f", inodes {pct(d['inodes_used_ratio'])} used"
        if d["missing"]:
            s += " (folder not created yet: measured on the disk that will hold it)"
        if d["warn"]:
            s += f"; WARN: at or above {snap['warn_used_ratio'] * 100:.0f}%"
        out.append(s)
    if not snap["disk"]:
        out.append("disk: no work or results folder given (--work-dir, --outdir, or weblog metadata)")
    h = snap["host"]
    s = f"memory host: {size(h['memory_bytes'])}"
    if h["max_memory_bytes"] is not None:
        s += f", max_memory {size(h['max_memory_bytes'])}"
    if h["warn"]:
        s += ("; WARN: max_memory is above this machine's memory, so a task can ask for more than there is "
              "(lower --max_memory; under WSL2 this machine is the VM, and Docker Desktop's containers get only "
              "the memory set in its settings)")
    out.append(s)
    for t in snap["tasks"]:
        name = f"task {t['process']}" + (f" (step {t['step']})" if t["step"] else "")
        if t["state"] == "running":
            s = f"{name}: running, asked for {size(t['request_bytes'])}"
        elif t["ratio"] is None:
            s = f"{name}: {t['state']}, no ratio (request or peak_rss unknown)"
        else:
            s = (f"{name}: peak_rss {size(t['peak_rss_bytes'])} of {size(t['request_bytes'])} asked "
                 f"({pct(t['ratio'])})")
            if t["warn"]:
                s += "; hint: close to its request (peak_rss counts shared pages, so this may be fine)"
        if t["oom_retry"]:
            s += f"; retried {t['oom_retry']} time{'s' if t['oom_retry'] > 1 else ''} after a memory or time kill"
        out.append(s)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--events", help="Nextflow weblog messages, one JSON object per line")
    ap.add_argument("--trace", action="append", default=[], help="a Nextflow trace file (repeatable)")
    ap.add_argument("--work-dir", help="the run's work folder (default: workDir from the weblog metadata)")
    ap.add_argument("--outdir", help="the run's results folder (default: outdir from the weblog metadata)")
    ap.add_argument("--max-memory", help="the pipeline's max_memory, e.g. 64.GB (default: from the weblog metadata)")
    ap.add_argument("--disk-warn-percent", help="warn at this disk use (default: PGP_DISK_WARN_PERCENT or 90)")
    ap.add_argument("--mem-warn-percent", help="hint at this peak_rss share of a request (default: PGP_MEM_WARN_PERCENT or 90)")
    a = ap.parse_args(argv)
    try:
        snap = snapshot(run_state.read_events(a.events), a.work_dir, a.outdir, a.disk_warn_percent,
                        a.mem_warn_percent, a.max_memory, [r for p in a.trace for r in run_eta.read_trace(p)])
    except ValueError as e:
        ap.error(str(e))
    sys.stdout.write("".join(s + "\n" for s in lines(snap)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
