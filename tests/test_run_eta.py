#!/usr/bin/env python3
"""bin/run_eta.py on synthetic traces (tests/fixtures/run_eta) and synthetic run states.

Checks:
  1. nextflow.config's trace.fields is Nextflow's 14 default columns in order, then
     process, tag, attempt, cpus, memory, start, complete;
  2. an old trace (default columns) and a new one (with the added columns) give the
     same step spans;
  3. '1h 2m 3s', '850ms', '18.3 GB' and a trace time parse equal to their raw values;
  4. cached and failed tasks, and a process of no step, never enter the history;
  5. ingesting twice, or a copy of a trace under another name, adds nothing;
  6. two boxes write two history files;
  7. the history holds no sample id and no '/';
  8. three runs give the min-max range with n=3, and only the last 5 runs count;
  9. with no history, the low end is the catalogue's measured lower end (hours_low)
     and the high end its upper end (hours), doubled on an 8-core box and not scaled
     below it on a 32-core box; on 16 and on 32 cores every measured step has
     low < high, never one figure;
 10. the hours left are the hand-computed longest chain, and the CPU-work bound when
     it is larger;
 11. a running step past its high end makes the high end unknown and says so;
 12. a step the catalogue does not time, or only estimates, lands in unknown_steps, never 0;
 13. a step with a task whose last attempt failed or was stopped gives no span;
 14. the command writes the estimate and the history from an OUTDIR;
 15. without --sample, the trace of a run that run-all.sh has not finished is left
     out of the step times, and counts once run_status.tsv says it finished; a
     run_status.tsv that cannot be read leaves out the newest trace, never a crash.
Each check is then run against a planted mistake (a mutant of run_eta.py or a
broken config) and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_run_eta.py
"""
import contextlib
import glob
import importlib.util
import io
import json
import os
import re
import shutil
import sys
import tempfile
import time

os.environ["TZ"] = "UTC"  # trace times are the host's local time
time.tzset()
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import run_eta  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "run_eta")
STEPS = run_state.read_steps()
H = 3600.0
FAILS = 0
TMP = None


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def mutant(**repl):
    """A fresh copy of bin/run_eta.py with the named functions or values replaced."""
    spec = importlib.util.spec_from_file_location("run_eta_mutant", os.path.join(REPO, "bin", "run_eta.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for k, v in repl.items():
        setattr(mod, k, v(mod) if callable(v) and getattr(v, "wants_mod", False) else v)
    return mod


def with_mod(f):
    f.wants_mod = True
    return f


def outdir(*traces):
    """A fresh OUTDIR whose pipeline_info holds the named fixture traces (name -> fixture)."""
    d = tempfile.mkdtemp(dir=TMP)
    os.makedirs(os.path.join(d, "pipeline_info"))
    for name, fx in traces:
        shutil.copy(os.path.join(FIX, fx), os.path.join(d, "pipeline_info", name))
    return d


def spans(history):
    out = {}
    for r in history:
        if r["step_span_s"]:
            out.setdefault(r["step"], set()).add(float(r["step_span_s"]))
    return out


def config_fields(text):
    block = re.search(r"^trace \{\n(.*?)^\}", text, re.S | re.M)
    m = re.search(r"^\s*fields\s*=\s*'([^']*)'", block.group(1), re.M) if block else None
    return m.group(1).split(",") if m else None


# --- the checks: each returns {name: (ok, detail)} ------------------------------------------

def config_checks(text):
    got = config_fields(text)
    return {"trace.fields: the 14 defaults in order, then the 7 added": (
        got == run_eta.DEFAULT_FIELDS + run_eta.ADDED_FIELDS and len(run_eta.DEFAULT_FIELDS) == 14
        and run_eta.ADDED_FIELDS == ["process", "tag", "attempt", "cpus", "memory", "start", "complete"], got)}


def trace_checks(m):
    out = {}
    old = m.ingest(outdir(("trace_2026-01-02_09-59-00.txt", "trace_old.txt")), "t", STEPS, 16)
    new = m.ingest(outdir(("trace_2026-01-02_09-59-00.txt", "trace_new.txt")), "t", STEPS, 16)
    want = {"01b": {20 * 60.0}, "02": {2 * 60.0, 4 * H}, "03": {11 * H + 10}, "04": {30 * 60.0}}
    out["an old and a new trace give the same spans"] = (spans(old) == spans(new) == want, (spans(old), spans(new)))

    out["human and raw units parse equal"] = (
        m.parse_duration("1h 2m 3s") == m.parse_duration("3723000") == 3723.0
        and m.parse_duration("850ms") == m.parse_duration("850") == 0.85
        and m.parse_memory("18.3 GB") == m.parse_memory(str(int(18.3 * 1024 ** 3)))
        and m.parse_memory("40 MB") == m.parse_memory(str(40 * 1024 ** 2))
        and m.parse_time("2026-01-02 10:00:00.000") == m.parse_time("1767348000000") == 1767348000.0,
        [m.parse_duration("1h 2m 3s"), m.parse_duration("850ms"), m.parse_memory("18.3 GB"),
         m.parse_time("2026-01-02 10:00:00.000")])

    hashes = [r["hash"] for r in new]
    out["cached, failed and unknown-process tasks never enter the history"] = (
        sorted(hashes) == ["aa000001", "aa000002", "aa000003", "aa000006", "aa000007"], hashes)

    d = outdir(("trace_2026-01-02_09-59-00.txt", "trace_new.txt"))
    m.ingest(d, "t", STEPS, 16)
    path = m.history_path(d, "t")
    first = open(path).read()
    m.ingest(d, "t", STEPS, 16)
    shutil.copy(os.path.join(FIX, "trace_new.txt"), os.path.join(d, "pipeline_info", "trace_2026-01-03_08-00-00.txt"))
    m.ingest(d, "t", STEPS, 16)
    out["ingesting twice, or a copied trace, adds nothing"] = (open(path).read() == first, open(path).read())

    m.ingest(d, "box-b", STEPS, 16)
    files = sorted(os.listdir(os.path.join(d, "pipeline_info", "step_history")))
    out["two boxes write two files"] = ([f for f in files if f.endswith(".tsv")] == ["box-b.tsv", "t.tsv"], files)

    body = open(path).read()
    out["the history holds no sample id and no '/'"] = (
        "demo1" not in body and "/" not in body and body.split("\n", 1)[0].split("\t") == m.HISTORY_COLUMNS, body)
    out["the box name keeps only letters, digits, '_' and '-'"] = (m.box_name("my box.local") == "my_box_local",
                                                                  m.box_name("my box.local"))
    return out


def hist(step, runs):
    """History rows: one run per (trace name, span hours)."""
    return [{"trace_file": f"trace_{i}.txt", "hash": f"h{i}", "attempt": "1", "step": step, "process": "X",
             "run_level": "no", "cpus": "8", "memory_bytes": "", "peak_rss_bytes": "", "realtime_s": "",
             "step_span_s": str(hours * H), "finished_utc": f"2026-01-{i + 1:02d}T00:00:00Z", "box_cpus": "16"}
            for i, hours in enumerate(runs)]


def range_checks(m):
    out = {}
    r = m.step_ranges(hist("02", [4, 3, 5]), STEPS, 16)["02"]
    out["three runs give min-max with n=3"] = (
        (r["low_s"], r["high_s"], r["n"], r["source"]) == (3 * H, 5 * H, 3, "box"), r)
    r = m.step_ranges(hist("02", [9, 3, 4, 5, 6, 7]), STEPS, 16)["02"]
    out["only the last 5 runs count"] = ((r["low_s"], r["high_s"], r["n"]) == (3 * H, 7 * H, 5), r)
    row = [s for s in STEPS if s["step"] == "03"][0]
    dv, dv_low = row["hours"] * H, float(row["hours_low"]) * H
    r8 = m.step_ranges([], STEPS, 8)["03"]
    r32 = m.step_ranges([], STEPS, 32)["03"]
    out["an 8-core box: the measured low end, the high end doubled"] = (
        r8["source"] == "docs" and r8["low_s"] == dv_low and r8["high_s"] == 2 * dv and "16 threads" in r8["label"], r8)
    out["a 32-core box is not scaled below the docs hours"] = (
        r32["low_s"] == dv_low and r32["high_s"] == dv, r32)
    measured = [s["step"] for s in STEPS if s["hours_from"] == "measured"]
    one = {c: [k for k in measured if not (m.step_ranges([], STEPS, c)[k]["source"] == "docs"
                                           and m.step_ranges([], STEPS, c)[k]["low_s"]
                                           < m.step_ranges([], STEPS, c)[k]["high_s"])] for c in (16, 32)}
    out["on 16 and 32 cores every measured step has low < high"] = (
        len(measured) >= 25 and not one[16] and not one[32], (len(measured), one))

    # Steps 14 and 04b have hours in the catalogue, but only the step page's estimate.
    est = m.step_ranges([], STEPS, 16)
    unk = m.remaining(state([("03", "done", None), ("14", "waiting", None)]), est, NOW, 16)
    out["a step the table only estimates gets no number"] = (
        est["14"]["source"] == est["04b"]["source"] == "none" and est["14"]["low_s"] is None
        and unk["unknown_steps"] == ["14"] and unk["low_s"] == 0, (est["14"], est["04b"], unk))

    # Step 07: the preprocessor ended, then PharmCAT was stopped (a killed run). Step 03: a
    # failed first attempt that a second one replaced still counts.
    def task(h, proc, status, ex, start, end, attempt="1"):
        return {"trace_file": "trace_x.txt", "hash": h, "attempt": attempt, "process": proc, "tag": "s1",
                "status": status, "exit": ex, "cpus": 1, "memory_bytes": None, "peak_rss_bytes": None,
                "realtime_s": end - start, "start": start, "complete": end}
    rows = m.history_rows([task("p1", "PHARMCAT_PREPROCESS", "COMPLETED", "0", 0.0, 40.0),
                           task("p2", "PHARMCAT", "ABORTED", "-", 41.0, 50.0),
                           task("d1", "DEEPVARIANT", "FAILED", "137", 0.0, 100.0),
                           task("d2", "DEEPVARIANT", "COMPLETED", "0", 110.0, 300.0, "2")],
                          m.process_steps(STEPS), 16)
    got = {(r["step"], r["process"]): r["step_span_s"] for r in rows}
    out["a step with a task that failed or was stopped gives no span"] = (
        got == {("07", "PHARMCAT_PREPROCESS"): None, ("03", "DEEPVARIANT"): 300.0}, got)
    return out


def state(rows, run="running"):
    """A run state with these steps: (step, state, started_utc); upstream from the catalogue."""
    up = {s["step"]: s["upstream"] for s in STEPS}
    return {"run": {"state": run}, "steps": [{"step": k, "state": st, "started_utc": t, "upstream": up[k]}
                                             for k, st, t in rows]}


NOW = 1767348000.0  # 2026-01-02T10:00:00Z
RANGES = {"02": {"low_s": 4 * H, "high_s": 5 * H, "cpus": 8},
          "03": {"low_s": 10 * H, "high_s": 12 * H, "cpus": 8},
          "13": {"low_s": 4 * H, "high_s": 5 * H, "cpus": 8},
          "23": {"low_s": 0.1 * H, "high_s": 0.2 * H, "cpus": 1},
          "04": {"low_s": 0.5 * H, "high_s": 1 * H, "cpus": 8},
          "33": {"low_s": None, "high_s": None, "cpus": None}}


def close(a, b):
    return a is not None and abs(a - b) < 1e-6


def remaining_checks(m):
    out = {}
    st = state([("02", "done", None), ("03", "waiting", None), ("13", "waiting", None),
                ("23", "waiting", None), ("04", "waiting", None)])
    # 03 -> 13 -> 23 is the longest chain: low 10 + 4 + 0.1, high 12 + 5 + 0.2.
    big = m.remaining(st, RANGES, NOW, 64)
    out["the hours left are the longest chain"] = (
        close(big["low_s"], 14.1 * H) and close(big["high_s"], 17.2 * H), big)
    # CPU work: 12*8 + 5*8 + 0.2*1 + 1*8 = 144.2 CPU-hours; over 4 cores 36.05 h.
    small = m.remaining(st, RANGES, NOW, 4)
    out["the CPU-work bound when it is larger"] = (close(small["low_s"], 14.1 * H) and close(small["high_s"], 36.05 * H),
                                                   small)

    late = m.remaining(state([("02", "done", None), ("03", "running", "2026-01-01T21:00:00Z"),
                              ("13", "waiting", None)]), RANGES, NOW, 64)
    out["a running step past its high end: high unknown, and said"] = (
        late["high_s"] is None and late["late_steps"] == ["03"]
        and late["note"] == "step 03 is taking longer than usual" and close(late["low_s"], 4 * H), late)
    going = m.remaining(state([("02", "done", None), ("03", "running", "2026-01-02T06:00:00Z"),
                               ("13", "waiting", None)]), RANGES, NOW, 64)
    out["a running step counts only what is left of it"] = (
        close(going["low_s"], 10 * H) and close(going["high_s"], 13 * H), going)

    ranges = m.step_ranges([], STEPS, 16)
    unk = m.remaining(state([("02", "done", None), ("33", "waiting", None)]), ranges, NOW, 16)
    out["a step with no time lands in unknown_steps, never 0"] = (
        unk["unknown_steps"] == ["33"] and ranges["33"]["low_s"] is None and ranges["33"]["source"] == "none", unk)
    done = m.remaining(state([("03", "waiting", None)], run="finished"), RANGES, NOW, 16)
    out["a finished run has nothing left"] = (done["low_s"] == 0 and done["high_s"] == 0, done)
    return out


def command_checks(m):
    d = outdir(("trace_2026-01-02_09-59-00.txt", "trace_new.txt"))
    o = os.path.join(d, "eta.json")
    with contextlib.redirect_stdout(io.StringIO()):
        rc = m.main([d, "--box", "ci", "--cores", "8", "--out", o])
    est = json.load(open(o)) if os.path.isfile(o) else {}
    ok = (rc == 0 and est.get("box") == "ci" and est["steps"]["02"]["source"] == "box"
          and est["steps"]["02"]["high_s"] == 4 * H and est["steps"]["13"]["source"] == "docs"
          and est["remaining"] is None and os.path.isfile(m.history_path(d, "ci")))
    out = {"the command writes the estimate and the history": (ok, est.get("steps", {}).get("02"))}

    # A run-all.sh run still going: its run_status.tsv has a start and no 'meta finished_utc'.
    d = outdir(("trace_2026-01-02_09-59-00.txt", "trace_new.txt"))
    logs = os.path.join(d, "demo1", "logs")
    os.makedirs(logs)
    status = os.path.join(logs, "run_status.tsv")
    with open(status, "w") as f:
        f.write(f"meta\tstarted_epoch\t{int(time.time()) - 60}\nstep\t02\tplanned\n")
    going = m.estimate(d, cores=16, box="ci")["steps"]["02"]["source"]
    with open(status, "a") as f:
        f.write("meta\tfinished_utc\t2026-01-02 12:00:00 UTC\n")
    ended = m.estimate(d, cores=16, box="ci")["steps"]["02"]["source"]
    out["without --sample, a run not finished yet is left out"] = ((going, ended) == ("docs", "box"), (going, ended))

    def unreadable(path):
        raise PermissionError(13, "Permission denied", path)
    keep = run_state.read_run_status
    run_state.read_run_status = unreadable
    try:
        got = m.estimate(d, cores=16, box="ci")["steps"]["02"]["source"]
    except OSError as e:
        got = repr(e)
    finally:
        run_state.read_run_status = keep
    out["an unreadable run_status.tsv leaves out the newest trace"] = (got == "docs", got)
    return out


def all_checks(m):
    r = {}
    r.update(trace_checks(m))
    r.update(range_checks(m))
    r.update(remaining_checks(m))
    r.update(command_checks(m))
    return r


def mutants():
    def forced_ok(mod):
        orig = mod.read_trace
        return lambda p: [dict(t, status="COMPLETED", exit="0") for t in orig(p)]

    def old_start_is_submit(mod):
        def f(p):
            ts = run_eta.read_trace(p)
            if "process" not in open(p).readline():
                for t, line in zip(ts, open(p).read().splitlines()[1:]):
                    t["start"] = run_eta.parse_time(line.split("\t")[6])
            return ts
        return f

    def slash(mod):
        orig = mod.history_rows
        return lambda *a: [dict(r, trace_file="pipeline_info/" + r["trace_file"]) for r in orig(*a)]

    def sample_in(mod):
        orig = mod.history_rows
        return lambda *a: [dict(r, process=r["process"] + "_demo1") for r in orig(*a)]

    def none_is_zero(mod):
        orig = mod.step_ranges
        return lambda *a, **k: {s: (dict(r, low_s=0.0, high_s=0.0) if r["source"] == "none" else r)
                                for s, r in orig(*a, **k).items()}

    def estimates_count(mod):
        orig = mod.step_ranges
        return lambda history, steps, *a, **k: orig(history, [dict(s, hours_from="measured") for s in steps], *a, **k)

    def no_low(mod):
        orig = mod.step_ranges
        return lambda history, steps, *a, **k: orig(history, [dict(s, hours_low="") for s in steps], *a, **k)

    def no_guard(mod):
        def f(outdir):
            out = set()
            for p in glob.glob(os.path.join(outdir, "*", "logs", "run_status.tsv")):
                meta = (run_state.read_run_status(p) or {}).get("meta", {})
                if "finished_utc" not in meta and "started_epoch" in meta:
                    d = os.path.dirname(os.path.dirname(p))
                    out |= {os.path.basename(t) for t in run_state.trace_files(d, float(meta["started_epoch"]))}
            return out
        return f

    def finished_ignored(mod):
        orig = mod.remaining
        return lambda st, *a: orig(dict(st, run={"state": "running"}), *a)

    def out_ignored(mod):
        orig = mod.main
        def f(argv=None):
            i = argv.index("--out")
            return orig(argv[:i] + argv[i + 2:])
        return f

    for f in (forced_ok, old_start_is_submit, slash, sample_in, none_is_zero, estimates_count, no_low, no_guard,
              finished_ignored, out_ignored):
        f.wants_mod = True
    return {
        "an old and a new trace give the same spans": mutant(read_trace=old_start_is_submit),
        "human and raw units parse equal": mutant(MEM_UNITS=dict(run_eta.MEM_UNITS, GB=1000 ** 3)),
        "cached, failed and unknown-process tasks never enter the history": mutant(read_trace=forced_ok),
        "ingesting twice, or a copied trace, adds nothing": mutant(known_tasks=lambda by_file, name: set()),
        "two boxes write two files": mutant(box_name=lambda box=None: "same"),
        "the history holds no sample id and no '/'": mutant(history_rows=slash),
        "three runs give min-max with n=3": mutant(LAST_RUNS=2),
        "only the last 5 runs count": mutant(LAST_RUNS=6),
        "an 8-core box: the measured low end, the high end doubled": mutant(DOCS_THREADS=8),
        "a 32-core box is not scaled below the docs hours": mutant(docs_scale=lambda c: 16.0 / c),
        "on 16 and 32 cores every measured step has low < high": mutant(step_ranges=no_low),
        "without --sample, a run not finished yet is left out": mutant(live_traces=lambda outdir: set()),
        "an unreadable run_status.tsv leaves out the newest trace": mutant(live_traces=no_guard),
        "the hours left are the longest chain": mutant(remaining=None),  # replaced below
        "the CPU-work bound when it is larger": mutant(cpu_work=lambda *a: 0.0),
        "a running step past its high end: high unknown, and said": mutant(
            left_of=lambda lo, hi, gone: (max(0.0, lo - gone), max(0.0, hi - gone), False)),
        "a step with no time lands in unknown_steps, never 0": mutant(step_ranges=none_is_zero),
        "a step the table only estimates gets no number": mutant(step_ranges=estimates_count),
        "a step with a task that failed or was stopped gives no span": mutant(finished=lambda ts: True),
        "a running step counts only what is left of it": mutant(left_of=lambda lo, hi, gone: (lo, hi, False)),
        "a finished run has nothing left": mutant(remaining=finished_ignored),
        "the command writes the estimate and the history": mutant(main=out_ignored),
        "the box name keeps only letters, digits, '_' and '-'": mutant(box_name=lambda box=None: box or "box"),
        "the history mutant for samples": mutant(history_rows=sample_in),
    }


def main():
    global TMP
    TMP = tempfile.mkdtemp(prefix="test_run_eta.")
    try:
        cfg = open(os.path.join(REPO, "nextflow.config")).read()
        results = config_checks(cfg)
        results.update(all_checks(run_eta))
        for name, (ok, detail) in results.items():
            check(name, ok, detail)

        print("\nNegative controls: each check against a planted mistake must fail")
        broken = cfg.replace(",memory,start", ",start")
        check("fails on a planted mistake: trace.fields without memory",
              not config_checks(broken)["trace.fields: the 14 defaults in order, then the 7 added"][0], "passed")
        ms = mutants()
        # The longest-chain check against a builder that ignores upstream steps.
        m = mutant()
        orig_remaining = m.remaining

        def no_upstream(st, *a):
            return orig_remaining(dict(st, steps=[dict(s, upstream=[]) for s in st["steps"]]), *a)
        m.remaining = no_upstream
        ms["the hours left are the longest chain"] = m
        sample_mutant = ms.pop("the history mutant for samples")
        for name, mod in ms.items():
            check(f"fails on a planted mistake: {name}", not all_checks(mod)[name][0], "the check passed on the mutant")
        check("fails on a planted mistake: a sample id in the history",
              not all_checks(sample_mutant)["the history holds no sample id and no '/'"][0], "passed")
        same = all_checks(mutant())
        check("an unchanged copy passes (the controls test the mistake, not the copy)",
              all(ok for ok, _ in same.values()), [k for k, (ok, _) in same.items() if not ok])
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
