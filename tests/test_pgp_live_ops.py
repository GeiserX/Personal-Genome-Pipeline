#!/usr/bin/env python3
"""bin/pgp_live_ops.py on the made-up run of tests/fixtures/live/ops, with a fake disk and a fake host.

Checks:
  1. half a minute into step 11 (usual time 60 to 100 s on this box), the run is not late and not
     stalled, and the hours left is a range with both ends;
  2. two hours later, with no message since, step 11 is late and the run is stalled;
  3. a run quiet past --stall-minutes whose running step is still inside its usual time is not stalled;
  4. --stall-minutes 300 keeps the quiet run of check 2 from being stalled;
  5. a failed task makes the run failed; a pgp_run_finished with exit 1 does too, exit 0 does not;
  6. the retried kill counts one retry, and 7.6 of 8 GB asked is a memory hint;
  7. max_memory 30 GB on a fake 16 GB host is a memory alert, on a fake 64 GB host it is not;
  8. a disk 95% used is a disk alert at the default 90, not with --disk-warn-percent 96;
  9. a bare 'nextflow run' (no run folder, folders from the weblog) gives the same metric series as
     the run-all.sh run;
 10. the result holds none of the run's paths, and the metrics carry only the labels run and fs;
 11. /metrics is open without a password, needs a login or the bearer token with one, and reads
     the token file only when a token is needed;
 12. every metric /metrics can write is in docs/alerts.md, and every metric the example rules use
     is one /metrics writes;
 13. the step history is read, never written;
 14. with bin/pgp_live.py present, the server loads the plugin: state.json has plugins.ops and
     GET /metrics answers in the Prometheus text format.
Checks 1 to 9 then run against planted mistakes (mutants of pgp_live_ops.py) and must fail
there: a check that cannot fail is no check.

Run: python3 tests/test_pgp_live_ops.py
"""
import contextlib
import importlib.util
import io
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import types
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
import run_eta  # noqa: E402
import run_headroom  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "live", "ops")
SAMPLE = "OPS-001"
T0 = 1767348000  # 2026-01-02T10:00:00Z, the fixture's start
G = 1024 ** 3
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1
    return ok


def load(path, name="pgp_live_ops"):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class FakeFs:
    def __init__(self, used_ratio):
        self.f_frsize = self.f_bsize = 4096
        self.f_blocks = 1000000
        self.f_bfree = self.f_bavail = int(self.f_blocks * (1 - used_ratio))
        self.f_files, self.f_ffree = 1000, 900


@contextlib.contextmanager
def disk_at(ratio):
    keep = run_headroom._statvfs
    run_headroom._statvfs = lambda p: FakeFs(ratio)
    try:
        yield
    finally:
        run_headroom._statvfs = keep


def make_tree(tmp, run_all=True, extra_events=()):
    """GENOME_DIR in TMP: the history of box 'opsbox', and for run-all.sh the run folder. Returns
    (root, run_dir or None, events, work, out)."""
    root = os.path.join(tmp, "genome")
    work, out = os.path.join(tmp, "nfwork"), os.path.join(tmp, "nfout")
    for d in (root, work, out, os.path.join(root, "pipeline_info", "step_history")):
        os.makedirs(d, exist_ok=True)
    shutil.copy(os.path.join(FIX, "step_history.tsv"), os.path.join(root, "pipeline_info", "step_history", "opsbox.tsv"))
    with open(os.path.join(FIX, "events.jsonl")) as f:
        text = f.read().replace("@WORK@", work).replace("@OUT@", out)
    events = [json.loads(line) for line in text.splitlines() if line.strip()] + list(extra_events)
    run_dir = None
    if run_all:
        run_dir = os.path.join(root, SAMPLE)
        os.makedirs(os.path.join(run_dir, "logs"), exist_ok=True)
        os.makedirs(os.path.join(run_dir, "nextflow", "work"), exist_ok=True)
        shutil.copy(os.path.join(FIX, "run_status.tsv"), os.path.join(run_dir, "logs", "run_status.tsv"))
        os.utime(os.path.join(run_dir, "logs", "run_status.tsv"), (T0, T0))
    return root, run_dir, events, work, out


def cfg_for(root, **kw):
    c = types.SimpleNamespace(root=root, steps=run_state.STEPS_TSV, eta_box="opsbox", eta_cores=8,
                              disk_warn_percent=None, mem_warn_percent=None, stall_minutes=60)
    for k, v in kw.items():
        setattr(c, k, v)
    return c


def ops_at(mod, now, run_all=True, extra=(), physical=64 * G, disk=0.50, **kw):
    tmp = tempfile.mkdtemp(prefix="pgp-ops-")
    try:
        root, run_dir, events, work, out = make_tree(tmp, run_all, extra)
        steps = run_state.read_steps()
        status = run_state.read_run_status(os.path.join(run_dir, "logs", "run_status.tsv")) if run_dir else None
        cur = mod.current_run(events)
        state = run_state.build(SAMPLE, steps, status, None, [], run_state.event_tasks(cur, SAMPLE), cur, None, now)
        with disk_at(disk):
            o = mod.run_ops(run_dir or os.path.join(root, SAMPLE), state, events, cfg_for(root, **kw),
                            now=now, steps=steps, physical=physical)
        return o, {"tmp": tmp, "root": root, "work": work, "out": out, "state": state}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def kinds(o):
    return {a["kind"]: a["level"] for a in o["alerts"]}


def series(text):
    """{metric{labels}} names (no values) of a Prometheus text."""
    return {line.rsplit(" ", 1)[0] for line in text.splitlines() if line and not line.startswith("#")}


def finished(code):
    return {"event": "pgp_run_finished", "exit": code, "utcTime": "2026-01-02T10:06:00Z"}


def behaviour_checks(mod, quiet=False):
    """Checks 1 to 9; returns how many failed."""
    global FAILS
    before = FAILS
    out = io.StringIO()
    with contextlib.redirect_stdout(out) if quiet else contextlib.nullcontext():
        early, _ = ops_at(mod, T0 + 330)
        h = early["hours_left"] or {}
        check("1. half a minute into step 11: not late, not stalled, the hours left has both ends",
              "late" not in kinds(early) and "stalled" not in kinds(early) and not early["run"]["stalled"]
              and h.get("low_s") is not None and h.get("high_s") is not None and h["low_s"] <= h["high_s"],
              early)
        late, _ = ops_at(mod, T0 + 302 + 7200)
        check("2. two hours later and quiet: step 11 is late and the run is stalled",
              late["hours_left"]["late_steps"] == ["11"] and late["hours_left"]["high_s"] is None
              and kinds(late).get("late") == "warning" and kinds(late).get("stalled") == "warning", late["alerts"])
        inside, _ = ops_at(mod, T0 + 302 + 90, stall_minutes=1)
        check("3. quiet past --stall-minutes but step 11 still inside its usual time: not stalled",
              inside["run"]["seconds_since_activity"] >= 60 and "stalled" not in kinds(inside)
              and "late" not in kinds(inside), (inside["run"], inside["alerts"]))
        calm, _ = ops_at(mod, T0 + 302 + 7200, stall_minutes=300)
        check("4. --stall-minutes 300 keeps a 2-hour quiet run from being stalled",
              "stalled" not in kinds(calm) and "late" in kinds(calm), calm["alerts"])
        failed_task = {"runName": "ops_demo", "event": "process_completed", "utcTime": "2026-01-02T10:06:00Z",
                       "trace": {"task_id": 4, "hash": "ab/000004", "name": f"CLINICAL:ROH ({SAMPLE})",
                                 "process": "CLINICAL:ROH", "tag": SAMPLE, "status": "FAILED", "exit": 1,
                                 "attempt": 1, "complete": (T0 + 360) * 1000}}
        f1, _ = ops_at(mod, T0 + 400, extra=[failed_task])
        f2, _ = ops_at(mod, T0 + 400, extra=[finished(1)])
        f0, _ = ops_at(mod, T0 + 400, extra=[finished(0)])
        check("5. a failed task fails the run; exit 1 fails it; exit 0 does not",
              f1["run"]["failed"] and kinds(f1).get("failed") == "critical" and f2["run"]["failed"]
              and f2["run"]["exit"] == 1 and not f0["run"]["failed"] and "failed" not in kinds(f0),
              (f1["run"], f2["run"], f0["run"]))
        m = early["memory"]
        check("6. one retry after exit 137, and 7.6 of 8 GB asked is a hint",
              m["retries"] == 1 and m["near_request"] == 1 and kinds(early).get("memory") == "warning",
              (m, early["alerts"]))
        small, _ = ops_at(mod, T0 + 330, physical=16 * G)
        check("7. max_memory 30 GB on a 16 GB host warns, on a 64 GB host it does not",
              small["memory"]["max_memory_above_ram"] and not early["memory"]["max_memory_above_ram"]
              and "max_memory" in [a["text"] for a in small["alerts"] if a["kind"] == "memory"][0], small["memory"])
        full, _ = ops_at(mod, T0 + 330, disk=0.95)
        full96, _ = ops_at(mod, T0 + 330, disk=0.95, disk_warn_percent="96")
        check("8. a disk 95% used warns at 90, not at --disk-warn-percent 96",
              kinds(full).get("disk") == "warning" and "disk" not in kinds(full96) and "disk" not in kinds(early)
              and full["disk"] and abs(full["disk"][0]["used_ratio"] - 0.95) < 0.001, (full["disk"], full96["alerts"]))
        bare, _ = ops_at(mod, T0 + 330, run_all=False)
        a = series(mod.metrics_text([(SAMPLE, early)], 0.9))
        b = series(mod.metrics_text([(SAMPLE, bare)], 0.9))
        check("9. a bare 'nextflow run' gives the same series as the run-all.sh run", a == b and len(a) > 15,
              sorted(a ^ b))
    return FAILS - before


def path_and_label_checks(mod):
    o, env = ops_at(mod, T0 + 330)
    blob = json.dumps(o)
    leaked = [k for k in ("tmp", "root", "work", "out") if env[k] in blob] + (["/"] if re.search(r'"/[^"]*"', blob) else [])
    text = mod.metrics_text([(SAMPLE, o)], 0.9, 64 * G)
    labels = {k for m in re.finditer(r"\{([^}]*)\}", text) for k in re.findall(r'(\w+)="', m.group(1))}
    bad = [line for line in text.splitlines() if line and not line.startswith("#")
           and not re.fullmatch(r'pgp_[a-z_]+(\{[^}]*\})? -?[0-9.e+-]+', line)]
    check("10. no path in the result; the metrics carry only the labels run and fs, in a valid format",
          not leaked and labels == {"run", "fs"} and not bad, (leaked, labels, bad[:3]))


def auth_checks(mod):
    def req(password=None, sid=None, header=None):
        r = types.SimpleNamespace(live=types.SimpleNamespace(password_hash=password), sid=sid,
                                  headers={"Authorization": header} if header else {})
        return r
    ok = (mod.metrics_allowed(req(), None)
          and not mod.metrics_allowed(req("h"), None)
          and not mod.metrics_allowed(req("h", header="Bearer s3cret"), None)
          and mod.metrics_allowed(req("h", header="Bearer s3cret"), "s3cret")
          and not mod.metrics_allowed(req("h", header="Bearer wrong"), "s3cret")
          and not mod.metrics_allowed(req("h", header="Basic s3cret"), "s3cret")
          and mod.metrics_allowed(req("h", sid="abc"), None))

    def unreadable():
        raise FileNotFoundError("no token file")
    try:
        mod.metrics_allowed(req("h", header="Bearer s3cret"), unreadable)
        raised = False
    except OSError:
        raised = True
    lazy = (mod.metrics_allowed(req("h", sid="abc"), unreadable) and mod.metrics_allowed(req(), unreadable)
            and mod.metrics_allowed(req("h", header="Bearer s3cret"), lambda: "s3cret"))
    check("11. /metrics: open without a password; with one, a login or the bearer token only; the token "
          "file is read only when a token is needed, and an unreadable one is an error", ok and raised and lazy,
          (ok, raised, lazy))


def doc_and_rule_checks(mod):
    with open(os.path.join(REPO, "docs", "alerts.md")) as f:
        doc = f.read()
    with open(os.path.join(REPO, "assets", "monitoring", "pgp-alerts.yml")) as f:
        rules = f.read()
    missing = [n for n in mod.METRICS if f"`{n}`" not in doc]
    used = set(re.findall(r"\bpgp_[a-z_]+\b", "\n".join(line.split("#")[0] for line in rules.splitlines()
                                                        if re.match(r"\s*expr:", line) or line.startswith(" " * 10))))
    unknown = sorted(n for n in used if n not in mod.METRICS)
    check("12. every metric is in docs/alerts.md; every metric the rules use exists",
          not missing and not unknown and used, (missing, unknown))


def write_check(mod):
    tmp = tempfile.mkdtemp(prefix="pgp-ops-")
    try:
        root, run_dir, events, _, _ = make_tree(tmp)
        hist = os.path.join(root, "pipeline_info", "step_history", "opsbox.tsv")
        with open(os.path.join(root, "pipeline_info", "trace_2026-01-02_10-00-00.txt"), "w") as f:
            f.write("task_id\thash\tname\tstatus\texit\n1\tzz/1\tROH (OPS-001)\tCOMPLETED\t0\n")
        before = {p: os.path.getmtime(os.path.join(dp, p)) for dp, _, fs in os.walk(root) for p in fs}
        with open(hist) as f:
            text = f.read()
        steps = run_state.read_steps()
        st = run_state.build(SAMPLE, steps, run_state.read_run_status(os.path.join(run_dir, "logs", "run_status.tsv")),
                             None, [], run_state.event_tasks(events, SAMPLE), events, None, T0 + 330)
        mod.run_ops(run_dir, st, events, cfg_for(root), now=T0 + 330, steps=steps, physical=64 * G)
        after = {p: os.path.getmtime(os.path.join(dp, p)) for dp, _, fs in os.walk(root) for p in fs}
        with open(hist) as f:
            same = f.read() == text
        check("13. the step history is read, never written", same and before == after, (sorted(set(after) ^ set(before))))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p


def server_check():
    if not os.path.isfile(os.path.join(BIN, "pgp_live.py")):
        print("[NOTE] 14. bin/pgp_live.py is not in this tree yet: the server check runs once it is")
        return
    tmp = tempfile.mkdtemp(prefix="pgp-ops-")
    proc = None
    try:
        root, run_dir, events, _, _ = make_tree(tmp)
        state_dir = os.path.join(tmp, "state")
        os.makedirs(os.path.join(state_dir, "runs", SAMPLE))
        with open(os.path.join(state_dir, "runs", SAMPLE, "events.jsonl"), "w") as f:
            f.write("".join(json.dumps(e) + "\n" for e in events))
        port, iport = free_port(), free_port()
        env = dict(os.environ, PGP_BOX="opsbox")
        env.pop("PGP_LIVE_PASSWORD_HASH", None)
        proc = subprocess.Popen([sys.executable, os.path.join(BIN, "pgp_live.py"), "--root", root, "--state-dir",
                                 state_dir, "--port", str(port), "--ingest-port", str(iport), "--quiet",
                                 "--stall-minutes", "30"], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        base = f"http://127.0.0.1:{port}"
        st = text = None
        for _ in range(100):
            try:
                with urllib.request.urlopen(base + f"/runs/{SAMPLE}/state.json", timeout=5) as r:
                    st = json.load(r)
                with urllib.request.urlopen(base + "/metrics", timeout=5) as r:
                    text, ctype = r.read().decode(), r.headers.get("Content-Type")
                break
            except OSError:
                time.sleep(0.1)
        ops = (st or {}).get("plugins", {}).get("ops")
        check("14. pgp_live.py loads the plugin: state.json has plugins.ops, /metrics answers",
              ops is not None and "alerts" in ops and text is not None and "version=0.0.4" in (ctype or "")
              and f'pgp_run_running{{run="{SAMPLE}"}} 1' in text and "pgp_runs 1" in text,
              (st and st.get("plugins"), (text or "")[:400], proc.poll() is not None and proc.stderr.read()))
    finally:
        if proc:
            proc.terminate()
            proc.wait(timeout=10)
        shutil.rmtree(tmp, ignore_errors=True)


# Planted mistakes: each must make behaviour_checks fail.
MUTANTS = [
    ("stall ignores the quiet time", "stalled = running and quiet is not None and quiet >= stall_s and not covered",
     "stalled = running and quiet is not None and not covered"),
    ("stall ignores the running steps", "quiet >= stall_s and not covered", "quiet >= stall_s"),
    ("late steps dropped", 'late = list(rem["late_steps"]) if rem else []', "late = []"),
    ("exit code ignored", "(exit_code not in (None, 0)) or error", "error"),
    ("disk warn ignored", 'disk_warn = [d for d in snap["disk"] if d["warn"]]', "disk_warn = []"),
    ("retries not counted", 'retries = sum(t["oom_retry"] for t in snap["tasks"])', "retries = 0"),
    ("bare run measures no folder", "else _abs_dir(m_work)", "else None"),
]


def mutant_checks():
    with open(os.path.join(BIN, "pgp_live_ops.py")) as f:
        src = f.read()
    tmp = tempfile.mkdtemp(prefix="pgp-ops-mut-")
    try:
        for name, old, new in MUTANTS:
            if old not in src:
                check(f"mutant '{name}': its line is in pgp_live_ops.py", False, old)
                continue
            p = os.path.join(tmp, "pgp_live_ops_mut.py")
            with open(p, "w") as f:
                f.write(src.replace(old, new, 1))
            global FAILS
            keep = FAILS
            n = behaviour_checks(load(p, "pgp_live_ops_mut"), quiet=True)
            FAILS = keep
            check(f"mutant '{name}' fails a check ({n} failed)", n > 0)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    mod = load(os.path.join(BIN, "pgp_live_ops.py"))
    behaviour_checks(mod)
    path_and_label_checks(mod)
    auth_checks(mod)
    doc_and_rule_checks(mod)
    write_check(mod)
    server_check()
    mutant_checks()
    print(f"\n{'FAILED: ' + str(FAILS) + ' check(s)' if FAILS else 'All checks passed.'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
