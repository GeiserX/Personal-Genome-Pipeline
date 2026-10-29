#!/usr/bin/env python3
"""run_state.py: where one sample's run stands, step by step, as one JSON file.

  run_state.py --sample S --sample-dir DIR [--trace FILE ...] [--events FILE]
               [--steps assets/steps.tsv] [--out DIR/run_state.json]

It reads only files that say how the run goes, never a result:
  DIR/logs/run_status.tsv       scripts/run-all.sh: when the run started, the
                                steps it planned ('planned'), skipped (with the
                                reason) and how each ended; 'meta finished_utc'
                                once run-all.sh has ended
  DIR/nextflow/process_status.tsv  main.nf, as the pipeline ends: the tasks of
                                each process that succeeded, came from the
                                cache, failed, were ignored or were left
                                unfinished
  the Nextflow trace            the trace_*.txt files of DIR/../pipeline_info
                                written since the run started (or --trace):
                                one row per finished task; only rows tagged
                                with this sample are read
  an events file (--events)     Nextflow weblog messages, one JSON object per
                                line: tasks submitted, started and completed
                                while the run goes
  DIR/summary.json              only each section's 'source' path, to tell a
                                result an earlier run left from one this run
                                wrote (its file time against the run start)

assets/steps.tsv says which processes and which summary sections belong to
each step. Standard library only. build() takes what the files hold and
returns the state, so a live server and the tests share it.

States of a step:
  waiting       planned, no task of it has started
  running       a task of it has started, or some have finished and others wait
  done          its tasks finished well (or run-all.sh recorded it ok)
  failed        a task of it failed for good (or run-all.sh recorded it failed)
  skipped       run-all.sh did not plan it; 'reason' says why
  not finished  the run ended before the step did
  not planned   no row and no task: a run without run-all.sh that does not
                show this step
A waiting or running step whose section file is older than the run start has
'earlier_result': that file is from an earlier run and is being redone.
"""
import argparse
import glob
import json
import os
import re
import sys
from datetime import datetime, timezone

SCHEMA_VERSION = 1
HERE = os.path.dirname(os.path.abspath(__file__))
STEPS_TSV = os.path.join(os.path.dirname(HERE), "assets", "steps.tsv")

# Trace and weblog task statuses
DONE_OK = {"COMPLETED", "CACHED"}
DONE_BAD = {"FAILED", "ABORTED"}


def iso(epoch):
    if epoch is None:
        return None
    return datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_row(v):
    """run-all.sh's 'YYYY-MM-DD HH:MM:SS UTC' as ISO 8601; anything else as it is."""
    try:
        return iso(datetime.strptime(v, "%Y-%m-%d %H:%M:%S UTC").replace(tzinfo=timezone.utc).timestamp())
    except (TypeError, ValueError):
        return v


# --- reading the files --------------------------------------------------------------------

def read_steps(path=STEPS_TSV):
    """assets/steps.tsv as a list of dicts, in file order."""
    rows, header = [], None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            r = dict(zip(header, cols + [""] * (len(header) - len(cols))))
            for k in ("processes", "sections", "upstream"):
                r[k] = [x for x in r[k].split(",") if x]
            r["hours"] = float(r["hours"]) if r["hours"] else None
            rows.append(r)
    return rows


def read_run_status(path):
    """logs/run_status.tsv: {'meta': {...}, 'steps': {step: last row}}; None without it."""
    if not os.path.isfile(path):
        return None
    st = {"meta": {}, "steps": {}}
    with open(path, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3 or c[0].startswith("#"):
                continue
            if c[0] == "meta":
                st["meta"][c[1]] = c[2]
            elif c[0] == "step":
                st["steps"][c[1]] = c[2]  # the last row wins: 'planned' gives way to the result
    return st


def read_process_status(path):
    """process_status.tsv: {process name without the workflow prefix: counts}."""
    if not os.path.isfile(path):
        return None
    out = {}
    with open(path, errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) != len(header):
                continue
            r = dict(zip(header, c))
            try:
                counts = {k: int(r[k]) for k in ("succeeded", "cached", "failed", "ignored", "unfinished")}
            except (KeyError, ValueError):
                continue
            out[r["process"].rsplit(":", 1)[-1]] = counts
    return out


def task_name(name):
    """'SV:MANTA (s1)' -> ('MANTA', 's1'); a task without a tag has tag None."""
    m = re.match(r"^\s*(\S+?)\s*(?:\((.*)\))?\s*$", name or "")
    if not m:
        return None, None
    return m.group(1).rsplit(":", 1)[-1], m.group(2)


def is_sample(tag, sample):
    """A task is the sample's when its tag is the sample id, or the id and ':<part>'."""
    return tag is not None and (tag == sample or tag.startswith(sample + ":"))


def parse_trace_time(v):
    """The trace's submit time, 'YYYY-MM-DD HH:MM:SS.mmm' in the host's local time."""
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(v, fmt).timestamp()
        except (TypeError, ValueError):
            pass
    return None


def parse_duration(v):
    """A trace duration ('1h 2m 3s', '450ms', '2.5s', '-') in seconds, or None."""
    if not v or v == "-":
        return None
    total, found = 0.0, False
    for num, unit in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", v):
        total += float(num) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
        found = True
    return total if found else None


def read_trace(paths, sample):
    """Rows of the trace files tagged with the sample, as task dicts."""
    tasks = []
    for p in paths:
        try:
            f = open(p, errors="replace")
        except OSError:
            continue
        with f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                r = dict(zip(header, line.rstrip("\n").split("\t")))
                proc, tag = task_name(r.get("name"))
                if not proc or not is_sample(tag, sample):
                    continue
                start = parse_trace_time(r.get("submit"))
                dur = parse_duration(r.get("duration"))
                tasks.append({"id": r.get("task_id") or r.get("hash") or r.get("name"), "name": r.get("name"),
                              "process": proc, "status": (r.get("status") or "").upper(),
                              "exit": r.get("exit"), "start": start,
                              "end": start + dur if start is not None and dur is not None else None})
    return tasks


def read_events(path):
    """An events file: one Nextflow weblog JSON message per line; a bad line is skipped."""
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open(path, errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except ValueError:
                continue
            if isinstance(ev, dict):
                out.append(ev)
    return out


def event_tasks(events, sample):
    """Weblog process_* messages of the sample's tasks, as task dicts, in file order."""
    status_of = {"process_submitted": "SUBMITTED", "process_started": "RUNNING"}
    tasks = []
    for ev in events:
        kind, tr = ev.get("event"), ev.get("trace")
        if not isinstance(tr, dict) or not str(kind).startswith("process_"):
            continue
        name = tr.get("name") or (f"{tr.get('process')} ({tr['tag']})" if tr.get("tag") else tr.get("process"))
        proc, tag = task_name(name)
        if not proc or not is_sample(tag, sample):
            continue
        st = status_of.get(kind) or str(tr.get("status") or "").upper()
        if kind == "process_completed" and st not in DONE_OK | DONE_BAD:
            st = "COMPLETED" if str(tr.get("exit")) == "0" else "FAILED"

        def ms(k):
            v = tr.get(k)
            return v / 1000.0 if isinstance(v, (int, float)) and v > 0 else None
        tasks.append({"id": str(tr.get("task_id") or tr.get("hash") or name), "name": name, "process": proc,
                      "status": st, "exit": tr.get("exit"), "start": ms("start") or ms("submit"),
                      "end": ms("complete")})
    return tasks


def run_finished_by_events(events):
    return any(ev.get("event") in ("completed", "error") for ev in events)


# --- the state ------------------------------------------------------------------------------

def latest_tasks(trace_tasks, ev_tasks):
    """One entry per task id (the trace first, then the events in order: the later wins), then
    one per task name: the attempt with the highest id is the task's outcome (a retry replaces
    a failed attempt)."""
    by_id = {}
    for t in trace_tasks + ev_tasks:
        old = by_id.get(t["id"])
        if old is not None:
            if old["status"] in DONE_OK | DONE_BAD and t["status"] not in DONE_OK | DONE_BAD:
                continue  # a late 'started' message does not undo a finished task
            t = dict(old, **{k: v for k, v in t.items() if v is not None})
        by_id[t["id"]] = t

    def order(t):
        try:
            return (0, int(t["id"]))
        except (TypeError, ValueError):
            return (1, str(t["id"]))
    by_name = {}
    for t in sorted(by_id.values(), key=order):
        by_name[t["name"]] = t
    return list(by_name.values())


def counts_of(tasks):
    c = {"waiting": 0, "running": 0, "done": 0, "cached": 0, "failed": 0}
    for t in tasks:
        s = t["status"]
        c["cached" if s == "CACHED" else "done" if s == "COMPLETED" else "failed" if s in DONE_BAD
          else "running" if s == "RUNNING" else "waiting"] += 1
    return c


def state_from_tasks(tasks):
    if not tasks:
        return None
    c = counts_of(tasks)
    if c["failed"]:
        return "failed"
    if c["running"] or (c["waiting"] and (c["done"] or c["cached"])):
        return "running"
    if c["waiting"]:
        return "waiting"
    return "done"


def state_from_process_status(rows):
    """The step from main.nf's end-of-run table (rows of the step's processes that have one)."""
    if not rows:
        return None
    if any(r["failed"] + r["ignored"] for r in rows):
        return "failed"
    if any(r["unfinished"] for r in rows):
        return "not finished"
    if any(r["succeeded"] + r["cached"] for r in rows):
        return "done"
    return "not finished"


def from_row(row):
    """run-all.sh's last row of a step -> (state, reason)."""
    if row is None:
        return None, None
    if row == "planned":
        return "planned", None
    if row == "ok":
        return "done", None
    if row.startswith("failed"):
        return "failed", None
    if row == "not finished":
        return "not finished", None
    m = re.match(r"^skipped(?: \((.*)\))?$", row)
    if m:
        return "skipped", m.group(1)
    return None, None


def build(sample, steps, status=None, process_status=None, trace_tasks=(), ev_tasks=(), events=(),
          sources=None, now=None):
    """The run state. steps: read_steps(); status: read_run_status(); process_status:
    read_process_status(); trace_tasks / ev_tasks: read_trace() / event_tasks() rows; events:
    read_events() (for the run's own start and end messages); sources: {section key: (path,
    mtime)} of the files the last summary read."""
    now = now if now is not None else datetime.now(tz=timezone.utc).timestamp()
    meta = (status or {}).get("meta", {})
    rows = (status or {}).get("steps", {})
    try:
        started = float(meta["started_epoch"]) if "started_epoch" in meta else None
    except ValueError:
        started = None
    finished = "finished_utc" in meta or (process_status is not None and status is None) or run_finished_by_events(events)
    tasks = latest_tasks(list(trace_tasks), list(ev_tasks))
    out_steps = []
    for s in steps:
        mine = [t for t in tasks if t["process"] in s["processes"]]
        state, reason = from_row(rows.get(s["step"]))
        note = None
        if state in (None, "planned"):
            ps = [process_status[p] for p in s["processes"] if process_status and p in process_status]
            live = state_from_process_status(ps) if ps else state_from_tasks(mine)
            if live:
                state = live
            elif state == "planned":
                state = "not finished" if finished else "waiting"
            else:
                state = "not planned"
        earlier = None
        if state in ("waiting", "running") and started is not None:
            for k in s["sections"]:
                src = (sources or {}).get(k)
                if src and src[1] is not None and src[1] < started:
                    earlier = {"section": k, "source": src[0], "file_date": iso(src[1])}
                    note = f"an earlier result ({iso(src[1])}) is being redone"
                    break
        starts = [t["start"] for t in mine if t.get("start") is not None]
        ends = [t["end"] for t in mine if t.get("end") is not None]
        procs = {p: counts_of([t for t in mine if t["process"] == p]) for p in s["processes"]
                 if any(t["process"] == p for t in mine)}
        out_steps.append({
            "step": s["step"], "label": s["label"], "tool": s["tool"] or None, "state": state,
            "reason": reason, "note": note, "earlier_result": earlier,
            "started_utc": iso(min(starts)) if starts else None,
            "ended_utc": iso(max(ends)) if ends and state in ("done", "failed") else None,
            "processes": procs, "sections": s["sections"], "upstream": s["upstream"],
            "hours": s["hours"], "hours_from": s["hours_from"] or None,
        })
    if finished:
        run_state = "finished"
    elif status is not None or tasks:
        run_state = "running"
    else:
        run_state = "not started"
    return {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "generated_utc": iso(now),
        "run": {"state": run_state, "started_utc": iso(started), "finished_utc": utc_row(meta.get("finished_utc")),
                "declared_sex": meta.get("declared_sex"), "status_file": status is not None},
        "steps": out_steps,
    }


# --- the command ----------------------------------------------------------------------------

def summary_sources(sample_dir):
    """{section: (source path, its file time now)} from DIR/summary.json; only 'source' is read."""
    p = os.path.join(sample_dir, "summary.json")
    try:
        with open(p) as f:
            secs = json.load(f).get("sections") or {}
    except (OSError, ValueError, AttributeError):
        return {}
    out = {}
    for k, v in secs.items():
        src = v.get("source") if isinstance(v, dict) else None
        if not src:
            continue
        path = os.path.join(sample_dir, src)
        try:
            out[k] = (src, os.path.getmtime(path))
        except OSError:
            continue
    return out


def trace_files(sample_dir, started):
    """The trace files of the run: those under ../pipeline_info written since it started
    (the newest one when the start is not known)."""
    found = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(sample_dir)), "pipeline_info",
                                          "trace_*.txt")), key=os.path.getmtime)
    if started is None:
        return found[-1:]
    return [p for p in found if os.path.getmtime(p) >= started]


def collect(sample, sample_dir, trace=None, events_path=None, steps_path=STEPS_TSV):
    status = read_run_status(os.path.join(sample_dir, "logs", "run_status.tsv"))
    try:
        started = float(status["meta"]["started_epoch"]) if status else None
    except (KeyError, ValueError):
        started = None
    pt_path = os.path.join(sample_dir, "nextflow", "process_status.tsv")
    pt = read_process_status(pt_path)
    if pt is not None and started is not None and os.path.getmtime(pt_path) < started:
        pt = None  # an earlier run's table
    events = read_events(events_path)
    return build(sample, read_steps(steps_path), status, pt,
                 read_trace(trace if trace else trace_files(sample_dir, started), sample),
                 event_tasks(events, sample), events, summary_sources(sample_dir))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sample-dir", required=True)
    ap.add_argument("--trace", action="append", help="a Nextflow trace file (default: the run's under ../pipeline_info)")
    ap.add_argument("--events", help="Nextflow weblog messages, one JSON object per line")
    ap.add_argument("--steps", default=STEPS_TSV, help="the steps catalogue (default: assets/steps.tsv)")
    ap.add_argument("--out", help="default: SAMPLE_DIR/run_state.json")
    a = ap.parse_args(argv)
    state = collect(a.sample, a.sample_dir, a.trace, a.events, a.steps)
    out = a.out or os.path.join(a.sample_dir, "run_state.json")
    tmp = out + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f, indent=1)
        f.write("\n")
    os.replace(tmp, out)
    print(f"Run state written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
