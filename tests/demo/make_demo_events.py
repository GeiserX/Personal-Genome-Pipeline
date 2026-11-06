#!/usr/bin/env python3
"""make_demo_events.py: write the invented weblog of a DEMO-001 run, for the live-site demo.

  make_demo_events.py [--out tests/fixtures/live/demo-001.weblog.jsonl]
  make_demo_events.py --check [--out FILE]

The messages are the ones Nextflow's weblog and scripts/run-all.sh post to
bin/owngenome_live.py while a run goes: 'started', then process_submitted,
process_started and process_completed for one task of each process of every
step in assets/steps.tsv that has processes (step 38, the family workflow,
runs on its own and is left out), then 'completed' and pgp_run_finished.

Nothing is recorded from a real run. A step starts when the steps it reads
from (the upstream column) have ended, and step 24, the report, when every
other step has. It lasts the middle of its measured range in assets/steps.tsv
(hours_low to hours), or 10 minutes when the table has no measured time; its
tasks are all submitted when it starts and run one after the other. A process
two steps share (PRS: steps 25 and 26) runs once, in the first of them. Two
made-up events show the run-health card: the first attempt of CLINVAR_SCREEN
is killed with exit 137 and a second attempt finishes, and HLA_TYPING ends at
92% of the memory it asked for.

--check fails when FILE differs from what this script writes now, so a change
to assets/steps.tsv that the fixture misses is caught (tests/test_live_demo.py
runs it). Standard library only.
"""
import argparse
import json
import os
import sys
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "bin"))
import run_state  # noqa: E402

SAMPLE = "DEMO-001"
RUN_NAME = "demo_run"
RUN_ID = "00000000-0000-0000-0000-00000000d001"
START = 1767344400  # 2026-01-02 09:00:00 UTC
OUT = os.path.join(REPO, "tests", "fixtures", "live", "demo-001.weblog.jsonl")
LEFT_OUT = {"38"}  # the family workflow: its own run, not a sample's
LAST = {"24"}      # the report reads every section, so it runs after every other step
DEFAULT_S = 600
GIB = 1 << 30
RETRIED = "CLINVAR_SCREEN"   # attempt 1 killed with exit 137, attempt 2 finishes
NEAR_LIMIT = "HLA_TYPING"    # peak_rss at 92% of its request


def utc(epoch):
    return datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def duration(s):
    """Seconds a step lasts in the demo: the middle of its measured range, else DEFAULT_S."""
    if s["hours"] is not None and s.get("hours_from") == "measured":
        low = float(s["hours_low"]) if s.get("hours_low") else s["hours"]
        return int(round((low + s["hours"]) / 2 * 3600))
    return DEFAULT_S


def schedule(steps):
    """[(step, start, end)] in start order: each step starts when its upstream steps have ended."""
    by = {s["step"]: s for s in steps}
    end = {}

    def finish(step, seen=()):
        if step in end:
            return end[step]
        s = by[step]
        ups = [u for u in (s["upstream"] if step not in LAST else [x for x in by if x not in LAST])
               if u in by and u not in seen]
        start = max([finish(u, seen + (step,)) for u in ups] or [START])
        end[step] = start + duration(s)
        return end[step]
    for s in steps:
        finish(s["step"])
    return sorted(((s["step"], end[s["step"]] - duration(s), end[s["step"]]) for s in steps),
                  key=lambda x: (x[1], x[0]))


def events(steps_path=run_state.STEPS_TSV):
    steps = [s for s in run_state.read_steps(steps_path) if s["processes"] and s["step"] not in LEFT_OUT]
    by = {s["step"]: s for s in steps}
    timed = []  # (epoch, order, message)
    task = [0]

    def msg(kind, t, trace):
        timed.append((t, len(timed), {"runName": RUN_NAME, "runId": RUN_ID, "event": kind, "utcTime": utc(t),
                                      "trace": trace}))

    def one_task(proc, submit, start, end, attempt=1, fail=False, peak_share=0.4):
        task[0] += 1
        cpus = 4 if end - start > 3600 else 2
        mem = 8 * GIB if end - start > 3600 else 2 * GIB
        base = {"task_id": task[0], "hash": f"de/{task[0]:06d}", "name": f"OWNGENOME:{proc} ({SAMPLE})",
                "process": f"OWNGENOME:{proc}", "tag": SAMPLE, "attempt": attempt, "cpus": cpus, "memory": mem,
                "submit": submit * 1000}
        msg("process_submitted", submit, dict(base, status="SUBMITTED"))
        msg("process_started", start + 1, dict(base, status="RUNNING", start=(start + 1) * 1000))
        done = dict(base, status="FAILED" if fail else "COMPLETED", exit=137 if fail else 0,
                    start=(start + 1) * 1000, complete=end * 1000, realtime=(end - start - 1) * 1000,
                    peak_rss=int(mem * (0.99 if fail else peak_share)))
        msg("process_completed", end, done)

    emitted = set()
    for step, start, end in schedule(steps):
        # A process two steps share (PRS: steps 25 and 26) runs once, in the first of them.
        procs = [p for p in by[step]["processes"] if p not in emitted]
        emitted.update(procs)
        if not procs:
            continue
        span = max(len(procs), (end - start) // len(procs))
        for i, proc in enumerate(procs):
            a = start + i * span
            b = end if i == len(procs) - 1 else a + span
            # Every task of the step is submitted when the step starts, so the step stays running
            # until its last process ends.
            if proc == RETRIED:
                mid = a + (b - a) // 3
                one_task(proc, start, a, mid, attempt=1, fail=True)
                one_task(proc, mid, mid + 5, b, attempt=2)  # Nextflow resubmits at once
            else:
                one_task(proc, start, a, b, peak_share=0.92 if proc == NEAR_LIMIT else 0.4)
    last = max(t for t, _, _ in timed)
    out = [{"runName": RUN_NAME, "runId": RUN_ID, "event": "started", "utcTime": utc(START),
            "metadata": {"parameters": {"outdir": "/genome", "max_cpus": 16, "max_memory": "12.GB"},
                         "workflow": {"commandLine": "nextflow run main.nf -resume"}}}]
    out += [m for _, _, m in sorted(timed, key=lambda x: (x[0], x[1]))]
    out.append({"runName": RUN_NAME, "runId": RUN_ID, "event": "completed", "utcTime": utc(last + 30),
                "metadata": {"workflow": {"success": True}}})
    out.append({"event": "pgp_run_finished", "exit": 0, "utcTime": utc(last + 60)})
    return out


def text(evs):
    return "".join(json.dumps(e, separators=(",", ":"), sort_keys=True) + "\n" for e in evs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--check", action="store_true", help="fail when --out differs from what would be written")
    a = ap.parse_args(argv)
    want = text(events())
    if a.check:
        try:
            with open(a.out) as f:
                have = f.read()
        except OSError:
            have = None
        if have != want:
            print(f"ERROR: {a.out} is not what tests/demo/make_demo_events.py writes now: run it again", file=sys.stderr)
            return 1
        print(f"{a.out}: up to date")
        return 0
    with open(a.out, "w") as f:
        f.write(want)
    print(f"Written: {a.out} ({want.count(chr(10))} messages)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
