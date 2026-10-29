#!/usr/bin/env python3
"""bin/run_state.py and assets/steps.tsv on synthetic run files.

Checks (sample demo1; the trace and weblog files are tests/fixtures/run_state):
  1. a planned step with no task is waiting;
  2. a weblog process_started makes its step running;
  3. process_completed with exit 0 for every task of a step makes it done;
     a task that failed and was retried with success is done;
  4. a failed task makes its step failed;
  5. a skipped row keeps its reason;
  6. a planned step whose section file predates the run is waiting with
     'earlier_result' (being redone), not stale; a file this run wrote is not;
  7. trace rows and weblog messages of another sample, and untagged shared
     tasks, are not counted;
  8. once run-all.sh has ended (finished_utc), a planned step with no result
     is not finished; main.nf's process_status.tsv decides a step without a
     result row; run-all.sh's result row wins over the tasks;
  9. the command writes run_state.json from the files and opens no result
     file (a planted VCF and the summary's results are not read);
 10. assets/steps.tsv against the code: every step run-all.sh plans or
     records has a row with the same label and --tools name, the recorded
     steps are collect_summary.RUN_ALL_STEPS, every SECTIONS key belongs to
     the row of its step, every process in modules/local belongs to a row,
     and the tools are main.nf's known tools.
Each check is then run against a planted mistake (a mutant of the builder or
of the catalogue) and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_run_state.py
"""
import builtins
import copy
import importlib.util
import json
import os
import re
import shutil
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "run_state")
SAMPLE = "demo1"
STARTED = 1767347940  # 2026-01-02 09:59:00 UTC, before every task of the fixtures
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def status(extra_rows="", finished=False):
    rows = {"meta": {"started_epoch": str(STARTED), "declared_sex": "female"}, "steps": {}}
    for line in ("step\t16\tplanned\nstep\t06\tplanned\nstep\t07\tplanned\nstep\t04\tplanned\n"
                 "step\t10\tplanned\nstep\t11\tplanned\nstep\t27\tplanned\n"
                 "step\t13\tskipped (data not installed: vep_cache/homo_sapiens/115_GRCh38/info.txt)\n"
                 + extra_rows).splitlines():
        c = line.split("\t")
        rows["steps"][c[1]] = c[2]
    if finished:
        rows["meta"]["finished_utc"] = "2026-01-02 11:00:00 UTC"
    return rows


def build(rs, steps, st, pt=None, sources=None):
    events = rs.read_events(os.path.join(FIX, "events.jsonl"))
    return rs.build(SAMPLE, steps, st, pt, rs.read_trace([os.path.join(FIX, "trace.txt")], SAMPLE),
                    rs.event_tasks(events, SAMPLE), events, sources, now=STARTED + 3600)


def by_step(state):
    return {s["step"]: s for s in state["steps"]}


# --- the checks: each returns {name: (ok, detail)} ------------------------------------------

def state_checks(rs, steps):
    out = {}
    sources = {"cpic": ("cpic/demo1_cpic.tsv", STARTED - 86400), "roh": ("roh/demo1_roh.txt", STARTED + 400)}
    s = by_step(build(rs, steps, status(), sources=sources))
    out["planned, no task: waiting"] = (s["27"]["state"] == "waiting", s["27"])
    out["process_started: running"] = (s["11"]["state"] == "running", s["11"])
    out["every task completed with exit 0: done"] = (s["07"]["state"] == "done" and s["06"]["state"] == "done"
                                                     and s["16"]["state"] == "done", [s[k]["state"] for k in ("07", "06", "16")])
    out["a failed task retried with success: done"] = (s["04"]["state"] == "done", s["04"])
    out["a failed task: failed"] = (s["10"]["state"] == "failed", s["10"])
    out["a skipped row keeps its reason"] = (
        s["13"]["state"] == "skipped" and s["13"]["reason"] == "data not installed: vep_cache/homo_sapiens/115_GRCh38/info.txt",
        s["13"])
    out["an older file of a planned step: being redone, not stale"] = (
        s["27"]["state"] == "waiting" and (s["27"]["earlier_result"] or {}).get("section") == "cpic"
        and "redone" in (s["27"]["note"] or ""), s["27"])
    out["a file this run wrote is no earlier result"] = (s["11"]["earlier_result"] is None, s["11"])
    out["another sample and untagged tasks are not counted"] = (
        s["11"]["processes"] == {"ROH": {"waiting": 0, "running": 1, "done": 0, "cached": 0, "failed": 0}}
        and s["25"]["state"] == "not planned", (s["11"]["processes"], s["25"]["state"]))
    out["a step's times come from its tasks"] = (
        s["07"]["started_utc"] == "2026-01-02T10:05:04Z" and s["07"]["ended_utc"] == "2026-01-02T10:05:55Z", s["07"])

    fin = build(rs, steps, status(finished=True))
    f = by_step(fin)
    out["run ended: a planned step without a result is not finished"] = (
        f["27"]["state"] == "not finished" and fin["run"]["state"] == "finished"
        and fin["run"]["finished_utc"] == "2026-01-02T11:00:00Z", (f["27"]["state"], fin["run"]))
    pt = {"ROH": {"succeeded": 0, "cached": 0, "failed": 0, "ignored": 0, "unfinished": 1},
          "CPIC_LOOKUP": {"succeeded": 1, "cached": 0, "failed": 0, "ignored": 0, "unfinished": 0}}
    p = by_step(build(rs, steps, status(finished=True), pt))
    out["process_status.tsv decides a step without a result"] = (
        p["11"]["state"] == "not finished" and p["27"]["state"] == "done", (p["11"]["state"], p["27"]["state"]))
    zero = {"succeeded": 0, "cached": 0, "failed": 0, "ignored": 0, "unfinished": 0}
    pt0 = {"SURVIVOR_PREP": dict(zero, succeeded=3), "SURVIVOR_MERGE": dict(zero),
           "MITO_EXTRACT_CHRM": dict(zero, succeeded=1), "MITO_HAPLOGROUP": dict(zero, succeeded=1),
           "HAPLOCHECK": dict(zero)}
    z = by_step(build(rs, steps, status(finished=True), pt0))
    z20 = by_step(build(rs, steps, status(finished=True), dict(pt0, MITO_VARIANTS=dict(zero, succeeded=1))))
    out["a process with no task: not finished, except haplocheck without step 20"] = (
        z["22"]["state"] == "not finished" and z["12"]["state"] == "done" and z20["12"]["state"] == "not finished",
        (z["22"]["state"], z["12"]["state"], z20["12"]["state"]))
    r = by_step(build(rs, steps, status("step\t10\tok\nstep\t11\tfailed\n", finished=True)))
    out["run-all.sh's result row wins over the tasks"] = (r["10"]["state"] == "done" and r["11"]["state"] == "failed",
                                                         (r["10"]["state"], r["11"]["state"]))
    return out


def command_checks(rs):
    """The command on a sample folder; records every file it opens."""
    out = {}
    w = tempfile.mkdtemp(prefix="run_state_")
    opened = []
    real_open = builtins.open
    try:
        d = os.path.join(w, SAMPLE)
        for sub in ("logs", "nextflow", "vcf", "cpic"):
            os.makedirs(os.path.join(d, sub))
        os.makedirs(os.path.join(w, "pipeline_info"))
        with open(os.path.join(d, "logs", "run_status.tsv"), "w") as f:
            f.write(f"# run-all.sh\nmeta\tstarted_epoch\t{int(time.time()) - 60}\nmeta\tdeclared_sex\tfemale\n"
                    "step\t27\tplanned\nstep\t06\tplanned\n")
        with open(os.path.join(d, "vcf", f"{SAMPLE}.vcf.gz"), "w") as f:
            f.write("not read\n")
        with open(os.path.join(d, "cpic", f"{SAMPLE}_cpic.tsv"), "w") as f:
            f.write("gene\tdiplotype\n")
        os.utime(os.path.join(d, "cpic", f"{SAMPLE}_cpic.tsv"), (time.time() - 86400,) * 2)
        with open(os.path.join(d, "summary.json"), "w") as f:
            json.dump({"sections": {"cpic": {"source": f"cpic/{SAMPLE}_cpic.tsv", "data": {"genes": ["not read"]}}}}, f)
        shutil.copy(os.path.join(FIX, "trace.txt"), os.path.join(w, "pipeline_info", "trace_2026-01-02_10-00-00.txt"))
        ev = os.path.join(w, "events.jsonl")
        shutil.copy(os.path.join(FIX, "events.jsonl"), ev)

        def spy(path, *a, **k):
            opened.append(os.path.abspath(str(path)))
            return real_open(path, *a, **k)
        builtins.open = spy
        try:
            rc = rs.main(["--sample", SAMPLE, "--sample-dir", d, "--events", ev])
        finally:
            builtins.open = real_open
        with open(os.path.join(d, "run_state.json")) as f:
            got = json.load(f)
        s = by_step(got)
        out["the command writes run_state.json"] = (
            rc == 0 and got["sample"] == SAMPLE and s["06"]["state"] == "done" and s["27"]["state"] == "waiting"
            and (s["27"]["earlier_result"] or {}).get("section") == "cpic", (rc, s.get("06"), s.get("27")))
        read = sorted({os.path.relpath(p, w) for p in opened if p.startswith(w) and not p.endswith(".tmp")})
        allowed = {f"{SAMPLE}/logs/run_status.tsv", f"{SAMPLE}/summary.json", "events.jsonl",
                   "pipeline_info/trace_2026-01-02_10-00-00.txt"}
        out["the command reads only status, trace, event and summary files"] = (set(read) <= allowed and
                                                                              f"{SAMPLE}/logs/run_status.tsv" in read, read)
    finally:
        builtins.open = real_open
        shutil.rmtree(w)
    return out


def run_all_steps():
    """(step, label, tool) of every plan call in scripts/run-all.sh, and the steps it records otherwise."""
    src = open(os.path.join(REPO, "scripts", "run-all.sh")).read()
    planned = set(re.findall(r'"(\d+b?) ([^":]+):(\w+)"', src))
    planned |= set(re.findall(r'plan "(\d+b?) ([^"]+)" (\w+)', src))
    other = set(re.findall(r"printf 'step\\t(\d+b?)\\tskipped", src)) - {s for s, _, _ in planned}
    scripts = set(re.findall(r"\b\w+:(\d+b?)-[\w-]+\b", src))  # GRIDSS:04b-gridss ...
    m = re.search(r'compgen -G "03\[([a-d-]+)\]', src)
    callers = {f"03{c}" for c in "abcd"} if m and m.group(1) == "a-d" else set()
    return planned, other | scripts | callers


def module_processes():
    names = set()
    for root, _, files in os.walk(os.path.join(REPO, "modules", "local")):
        for fn in files:
            if fn.endswith(".nf"):
                names |= set(re.findall(r"^process\s+([A-Z0-9_]+)", open(os.path.join(root, fn)).read(), re.M))
    return names


def main_nf_tools():
    src = open(os.path.join(REPO, "main.nf")).read()
    block = re.search(r"def known_tools = \[(.*?)\]", src, re.S).group(1)
    return set(re.findall(r"'(\w+)'", block))


def catalogue_checks(steps, sections, run_all_set):
    out = {}
    row = {s["step"]: s for s in steps}
    planned, other = run_all_steps()
    bad = [(n, l, t, row.get(n, {}).get("label"), row.get(n, {}).get("tool")) for n, l, t in sorted(planned)
           if n not in row or row[n]["label"] != l or row[n]["tool"] != t]
    out["every run-all.sh plan step has its row, label and tool"] = (len(planned) >= 40 and not bad, bad or len(planned))
    recorded = {s["step"] for s in steps if s["recorded"] == "yes"}
    want = {n for n, _, _ in planned} | other
    out["the recorded steps are what run-all.sh records"] = (recorded == want, sorted(recorded ^ want))
    out["RUN_ALL_STEPS is the recorded steps"] = (set(run_all_set) == recorded, sorted(set(run_all_set) ^ recorded))
    owner = {k: s["step"] for s in steps for k in s["sections"]}
    keys = {k: step for k, _, step, _ in sections}
    wrong = [(k, step, owner.get(k)) for k, step in keys.items() if owner.get(k) != step]
    out["every SECTIONS key belongs to the row of its step"] = (not wrong and set(owner) <= set(keys),
                                                                (wrong, sorted(set(owner) - set(keys))))
    procs = {}
    for s in steps:
        for p in s["processes"]:
            procs.setdefault(p, []).append(s["step"])
    mods = module_processes()
    lost = sorted(mods - set(procs))
    twice = {p: v for p, v in procs.items() if len(v) > 1 and p != "PRS"}  # step 26 runs inside PRS
    out["every module process belongs to a step"] = (len(mods) > 100 and not lost and not twice
                                                    and set(procs) <= mods, (lost, twice, sorted(set(procs) - mods)))
    tools = {s["tool"] for s in steps if s["tool"]}
    out["the tools are main.nf's known tools"] = (tools == main_nf_tools(), sorted(tools ^ main_nf_tools()))
    nums = [s["step"] for s in steps]
    out["one row per step, every upstream step has a row"] = (
        len(nums) == len(set(nums)) and all(u in row for s in steps for u in s["upstream"]), nums)
    return out


# --- the planted mistakes ---------------------------------------------------------------------

def mutant(**repl):
    """A fresh copy of bin/run_state.py with the named functions or values replaced."""
    spec = importlib.util.spec_from_file_location("run_state_mutant", os.path.join(REPO, "bin", "run_state.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for k, v in repl.items():
        setattr(mod, k, v)
    return mod


def mutants(steps):
    def no_running(tasks):
        st = run_state.state_from_tasks(tasks)
        return "waiting" if st == "running" else st

    def no_reason(row):
        st, _ = run_state.from_row(row)
        return st, None

    def first_attempt(trace, ev):
        return list({t["name"]: t for t in reversed(trace + ev)}.values())
    m = {
        "planned, no task: waiting": mutant(from_row=lambda row: (None, None) if row == "planned" else run_state.from_row(row)),
        "process_started: running": mutant(state_from_tasks=no_running),
        "every task completed with exit 0: done": mutant(event_tasks=lambda ev, sample: run_state.event_tasks(
            [e for e in ev if e.get("event") != "process_completed"], sample)),
        "a failed task retried with success: done": mutant(latest_tasks=first_attempt),
        "a failed task: failed": mutant(DONE_BAD=set()),
        "a skipped row keeps its reason": mutant(from_row=no_reason),
        "another sample and untagged tasks are not counted": mutant(is_sample=lambda tag, sample: True),
        "run ended: a planned step without a result is not finished": mutant(run_finished_by_events=lambda ev: False,
                                                                             utc_row=lambda v: None),
        "process_status.tsv decides a step without a result": mutant(state_from_process_status=lambda rows, ps=None: None),
        "a process with no task: not finished, except haplocheck without step 20": mutant(
            no_haplocheck_task_expected=lambda ps: True,
            state_from_process_status=lambda rows, ps=None: run_state.state_from_process_status(
                {p: r for p, r in rows.items() if any(r.values())} or rows, ps)),
        "an older file of a planned step: being redone, not stale": mutant(
            build=lambda sample, steps, st, pt, tt, et, ev, sources=None, now=None:
            run_state.build(sample, steps, st, pt, tt, et, ev, None, now)),
        "a file this run wrote is no earlier result": mutant(
            build=lambda sample, steps, st, pt, tt, et, ev, sources=None, now=None:
            run_state.build(sample, steps, st, pt, tt, et, ev, {k: (v[0], 0) for k, v in (sources or {}).items()}, now)),
        "a step's times come from its tasks": mutant(parse_trace_time=lambda v: None, event_tasks=lambda ev, sample: [
            dict(t, start=None, end=None) for t in run_state.event_tasks(ev, sample)]),
        "run-all.sh's result row wins over the tasks": mutant(
            from_row=lambda row: ("planned", None) if row in ("ok", "failed") else run_state.from_row(row)),
    }
    return m


def main():
    steps = run_state.read_steps()
    # 1-8 on the builder, 9 on the command, 10 on the catalogue
    results = {}
    results.update(state_checks(run_state, steps))
    results.update(command_checks(run_state))
    results.update(catalogue_checks(steps, collect_summary.SECTIONS, collect_summary.RUN_ALL_STEPS))
    for name, (ok, detail) in results.items():
        check(name, ok, detail)

    print("\nNegative controls: each check against a planted mistake must fail")
    for name, mod in mutants(steps).items():
        got = state_checks(mod, steps)[name]
        check(f"fails on a planted mistake: {name}", not got[0], "the check passed on the mutant")
    check("fails on a planted mistake: the command reports an earlier result",
          not command_checks(mutant(summary_sources=lambda d: {}))["the command writes run_state.json"][0],
          "passed on the mutant")

    def opens_vcf(d):
        with open(os.path.join(d, "vcf", f"{SAMPLE}.vcf.gz")) as f:
            f.read()
        return run_state.summary_sources(d)
    got = command_checks(mutant(summary_sources=opens_vcf))["the command reads only status, trace, event and summary files"]
    check("fails on a planted mistake: the command opens a VCF", not got[0], "passed on the mutant")
    same = state_checks(mutant(), steps)
    check("an unchanged copy passes (the controls test the mistake, not the copy)",
          all(ok for ok, _ in same.values()), [k for k, (ok, _) in same.items() if not ok])

    cat = [
        ("every run-all.sh plan step has its row, label and tool", [s for s in steps if s["step"] != "49"], None),
        ("every run-all.sh plan step has its row, label and tool",
         [dict(s, label="Renamed") if s["step"] == "07" else s for s in steps], None),
        ("the recorded steps are what run-all.sh records", [dict(s, recorded="no") if s["step"] == "04b" else s for s in steps], None),
        ("RUN_ALL_STEPS is the recorded steps", steps, set(collect_summary.RUN_ALL_STEPS) - {"33"}),
        ("every SECTIONS key belongs to the row of its step",
         [dict(s, sections=[k for k in s["sections"] if k != "mosaic_sex_loss"]) for s in steps], None),
        ("every module process belongs to a step",
         [dict(s, processes=[p for p in s["processes"] if p != "Y_STR"]) for s in steps], None),
        ("the tools are main.nf's known tools", [dict(s, tool="") if s["step"] == "34" else s for s in steps], None),
        ("one row per step, every upstream step has a row", steps + [copy.deepcopy(steps[0])], None),
    ]
    for name, st, ras in cat:
        got = catalogue_checks(st, collect_summary.SECTIONS, collect_summary.RUN_ALL_STEPS if ras is None else ras)[name]
        check(f"fails on a planted catalogue mistake: {name}", not got[0], "the check passed on the mutant")

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
