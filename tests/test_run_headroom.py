#!/usr/bin/env python3
"""bin/run_headroom.py on synthetic disks, synthetic tasks (tests/fixtures/run_headroom) and a fake host.

Checks:
  1. a disk 91% used warns at the default 90, one 89% used does not;
  2. --disk-warn-percent 80 (and PGP_DISK_WARN_PERCENT=80) warns at 85%;
  3. a work folder not created yet is measured on its nearest existing parent;
  4. two folders on one filesystem are reported once, on two filesystems twice;
  5. a filesystem with f_files 0 gives no inode ratio;
  6. 7.6 of 8 GB asked warns, 6 of 8 does not;
  7. '7.6 GB' (a trace) and raw bytes (a weblog message) give the same ratio;
  8. a trace with no memory column gives no ratio;
  9. exit 137 then attempt 2 counts one retry, and a task still running shows its request only;
 10. max_memory 64 GB on a fake 32 GB host warns, on a fake 128 GB host does not;
 11. the snapshot holds no '/' in any string value, and takes its folders from the weblog metadata;
 12. MEMORY_EXITS are the exits conf/base.config's errorStrategy retries;
 13. a real statvfs on a temp folder gives a used ratio in 0..1 and the command exits 0.
Each check is then run against a planted mistake (a mutant of run_headroom.py or a
broken config) and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_run_headroom.py
"""
import contextlib
import importlib.util
import io
import os
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import run_eta  # noqa: E402
import run_headroom  # noqa: E402
import run_state  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "run_headroom")
G = 1024 ** 3
FAILS = 0
TMP = None


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def mutant(**repl):
    """A fresh copy of bin/run_headroom.py with the named values replaced."""
    spec = importlib.util.spec_from_file_location("run_headroom_mutant", os.path.join(REPO, "bin", "run_headroom.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for k, v in repl.items():
        setattr(mod, k, v(mod) if getattr(v, "wants_mod", False) else v)
    return mod


def with_mod(f):
    f.wants_mod = True
    return f


class FakeVfs:
    def __init__(self, used_pct, files=1000, ffree=900):
        self.f_frsize = self.f_bsize = 4096
        self.f_blocks = 1000
        self.f_bfree = self.f_bavail = 1000 - used_pct * 10
        self.f_files, self.f_ffree = files, ffree


@contextlib.contextmanager
def fake_disks(m, vfs, devs=None):
    """Every path statvfs-es to VFS (or VFS[path]); st_dev is DEVS[path] or 1."""
    old = m._statvfs, m._dev
    m._statvfs = lambda p: vfs[p] if isinstance(vfs, dict) else vfs
    m._dev = lambda p: (devs or {}).get(p, 1)
    try:
        yield
    finally:
        m._statvfs, m._dev = old


def run_main(m, argv, env=None):
    out, saved = io.StringIO(), dict(os.environ)
    os.environ.update(env or {})
    try:
        with contextlib.redirect_stdout(out):
            rc = m.main(argv)
    finally:
        os.environ.clear()
        os.environ.update(saved)
    return rc, out.getvalue()


def events():
    return run_state.read_events(os.path.join(FIX, "events.jsonl"))


def by_process(tasks):
    return {t["process"]: t for t in tasks}


def strings(o):
    if isinstance(o, str):
        yield o
    elif isinstance(o, dict):
        for v in o.values():
            yield from strings(v)
    elif isinstance(o, list):
        for v in o:
            yield from strings(v)


def base_config_exits(text):
    m = re.search(r"errorStrategy\s*=\s*\{\s*task\.exitStatus in \[([\d,\s]+)\]", text)
    return {x.strip() for x in m.group(1).split(",")} if m else None


# --- the checks: each returns {name: (ok, detail)} ------------------------------------------

def disk_checks(m):
    out = {}
    d = tempfile.mkdtemp(dir=TMP)
    with fake_disks(m, FakeVfs(91)):
        hi = m.disk([("work", d)])
    with fake_disks(m, FakeVfs(89)):
        lo = m.disk([("work", d)])
    out["91% warns at the default, 89% does not"] = (
        len(hi) == 1 and hi[0]["warn"] and abs(hi[0]["used_ratio"] - 0.91) < 1e-9 and len(lo) == 1 and not lo[0]["warn"],
        (hi, lo))

    with fake_disks(m, FakeVfs(85)):
        _, plain = run_main(m, ["--work-dir", d])
        _, arg = run_main(m, ["--work-dir", d, "--disk-warn-percent", "80"])
        _, env = run_main(m, ["--work-dir", d], {"PGP_DISK_WARN_PERCENT": "80"})
    out["--disk-warn-percent 80 and PGP_DISK_WARN_PERCENT=80 warn at 85%"] = (
        "WARN" not in plain.splitlines()[0] and "WARN: at or above 80%" in arg and "WARN: at or above 80%" in env,
        (plain, arg, env))

    seen = []
    with fake_disks(m, FakeVfs(50)):
        old = m._statvfs
        m._statvfs = lambda p: seen.append(p) or old(p)
        miss = m.disk([("work", os.path.join(d, "not", "yet"))])
    out["a missing work folder is measured on its parent"] = (
        len(miss) == 1 and miss[0]["missing"] and seen == [d], (miss, seen))

    a, b = os.path.join(d, "a"), os.path.join(d, "b")
    os.makedirs(a)
    os.makedirs(b)
    with fake_disks(m, FakeVfs(50), {a: 7, b: 7}):
        one = m.disk([("work", a), ("results", b)])
    with fake_disks(m, FakeVfs(50), {a: 7, b: 8}):
        two = m.disk([("work", a), ("results", b)])
    out["one st_dev reports once, two report twice"] = (
        [e["fs"] for e in one] == ["work and results"] and [e["fs"] for e in two] == ["work", "results"], (one, two))

    with fake_disks(m, FakeVfs(50, files=0, ffree=0)):
        noi = m.disk([("work", d)])
    with fake_disks(m, FakeVfs(50)):
        withi = m.disk([("work", d)])
    out["f_files 0 gives no inode ratio"] = (
        noi[0]["inodes_used_ratio"] is None and abs(withi[0]["inodes_used_ratio"] - 0.1) < 1e-9, (noi, withi))

    real = m.disk([("work", d)])
    rc, text = run_main(m, ["--work-dir", d, "--outdir", d])
    out["a real statvfs on a temp folder gives a ratio in 0..1; the command exits 0"] = (
        len(real) == 1 and real[0]["used_ratio"] is not None and 0 <= real[0]["used_ratio"] <= 1
        and real[0]["size_bytes"] > 0 and rc == 0 and text.startswith("disk work and results: "), (real, rc, text))
    return out


def memory_checks(m):
    out = {}
    trace = run_eta.read_trace(os.path.join(FIX, "trace.txt"))
    tt = by_process(m.task_memory(trace))
    out["7.6 of 8 GB warns, 6 of 8 does not"] = (
        tt["FASTP"]["warn"] and abs(tt["FASTP"]["ratio"] - 0.95) < 1e-6
        and not tt["VEP"]["warn"] and abs(tt["VEP"]["ratio"] - 0.75) < 1e-6, tt)

    ev = by_process(m.task_memory(m.event_rows(events())))
    out["'7.6 GB' and raw bytes give the same ratio"] = (
        tt["FASTP"]["ratio"] is not None and abs(tt["FASTP"]["ratio"] - ev["FASTP"]["ratio"]) < 1e-6
        and ev["FASTP"]["step"] == "01b", (tt["FASTP"], ev["FASTP"]))

    nomem = m.task_memory(run_eta.read_trace(os.path.join(FIX, "trace_no_memory.txt")))
    out["no memory column gives no ratio"] = (
        len(nomem) == 1 and nomem[0]["ratio"] is None and not nomem[0]["warn"], nomem)

    out["exit 137 then attempt 2 counts one retry; a running task shows its request only"] = (
        ev["VEP"]["oom_retry"] == 1 and ev["VEP"]["attempt"] == 2 and abs(ev["VEP"]["ratio"] - 9 / 16) < 1e-6
        and ev["FASTP"]["oom_retry"] == 0
        and ev["DEEPVARIANT"]["state"] == "running" and ev["DEEPVARIANT"]["ratio"] is None
        and ev["DEEPVARIANT"]["request_bytes"] == 16 * G and ev["DEEPVARIANT"]["peak_rss_bytes"] is None, ev)

    small = m.host_memory("64.GB", physical=32 * G)
    big = m.host_memory("64.GB", physical=128 * G)
    out["64 GB max on a 32 GB host warns, on 128 GB not"] = (
        small["warn"] and small["max_memory_bytes"] == 64 * G and not big["warn"], (small, big))
    return out


def snapshot_checks(m):
    out = {}
    d = tempfile.mkdtemp(dir=TMP)
    work, res = os.path.join(d, "work"), os.path.join(d, "results")
    os.makedirs(work)
    os.makedirs(res)
    meta = {"event": "started", "metadata": {"workflow": {"workDir": work},
                                              "parameters": {"outdir": res, "max_memory": "64.GB"}}}
    with fake_disks(m, FakeVfs(50), {work: 1, res: 2}):
        snap = m.snapshot([meta] + events(), physical=32 * G,
                          trace_rows=run_eta.read_trace(os.path.join(FIX, "trace.txt")))
    bad = [s for s in strings(snap) if "/" in s]
    out["the snapshot holds no '/' and takes its folders from the metadata"] = (
        not bad and [e["fs"] for e in snap["disk"]] == ["work", "results"] and snap["host"]["warn"]
        and snap["warn_used_ratio"] == 0.9 and snap["mem_warn_ratio"] == 0.9
        and {t["process"] for t in snap["tasks"]} == {"FASTP", "VEP", "DEEPVARIANT"}, (bad, snap))
    return out


def all_checks(m):
    out = {}
    for f in (disk_checks, memory_checks, snapshot_checks):
        out.update(f(m))
    return out


def mutants():
    """{check name: a module with a planted mistake that check must catch}."""
    ms = {}
    ms["91% warns at the default, 89% does not"] = mutant(DEFAULT_WARN_PERCENT=95)

    @with_mod
    def ignore_percent(mod):
        orig = mod.warn_percent
        return lambda value=None, env=None: orig(None, None)
    ms["--disk-warn-percent 80 and PGP_DISK_WARN_PERCENT=80 warn at 85%"] = mutant(warn_percent=ignore_percent)
    ms["a missing work folder is measured on its parent"] = mutant(nearest_existing=os.path.abspath)

    @with_mod
    def no_merge(mod):
        orig = mod.disk
        return lambda paths, warn=None: [e for p in paths for e in orig([p], warn)]
    ms["one st_dev reports once, two report twice"] = mutant(disk=no_merge)

    @with_mod
    def inodes_always(mod):
        orig = mod._measure

        def measure(path):
            used, free, size, inodes = orig(path)
            return used, free, size, inodes if inodes is not None else 0.0
        return measure
    ms["f_files 0 gives no inode ratio"] = mutant(_measure=inodes_always)
    ms["a real statvfs on a temp folder gives a ratio in 0..1; the command exits 0"] = mutant(
        _statvfs=lambda p: FakeVfs(150))

    @with_mod
    def strict_above(mod):
        orig = mod.task_memory

        def tm(rows, warn=None, steps=None):
            out = orig(rows, warn, steps)
            for t in out:
                t["warn"] = t["ratio"] is not None and t["ratio"] > 0.99
            return out
        return tm
    ms["7.6 of 8 GB warns, 6 of 8 does not"] = mutant(task_memory=strict_above)

    @with_mod
    def decimal_peak(mod):
        orig = mod._row

        def row(r, i):
            t = orig(r, i)
            if "memory_bytes" not in r and t["peak"] is not None:
                t["peak"] = t["peak"] * 1000 ** 3 // 1024 ** 3  # weblog bytes read as decimal units
            return t
        return row
    ms["'7.6 GB' and raw bytes give the same ratio"] = mutant(_row=decimal_peak)

    @with_mod
    def zero_request_is_8g(mod):
        orig = mod._row

        def row(r, i):
            t = orig(r, i)
            t["memory"] = t["memory"] or 8 * G
            return t
        return row
    ms["no memory column gives no ratio"] = mutant(_row=zero_request_is_8g)
    ms["exit 137 then attempt 2 counts one retry; a running task shows its request only"] = mutant(
        MEMORY_EXITS={"143"})
    ms["64 GB max on a 32 GB host warns, on 128 GB not"] = mutant(
        host_memory=lambda max_memory=None, physical=None: {"memory_bytes": physical, "max_memory_bytes": None,
                                                            "warn": False})

    @with_mod
    def path_label(mod):
        orig = mod.disk
        return lambda paths, warn=None: orig([(p or label, p) for label, p in paths], warn)
    ms["the snapshot holds no '/' and takes its folders from the metadata"] = mutant(disk=path_label)
    return ms


def main():
    global TMP
    TMP = tempfile.mkdtemp(prefix="test_run_headroom.")
    try:
        base = open(os.path.join(REPO, "conf", "base.config")).read()
        got = base_config_exits(base)
        check("MEMORY_EXITS are the exits conf/base.config retries", got == run_headroom.MEMORY_EXITS, got)
        results = all_checks(run_headroom)
        for name, (ok, detail) in results.items():
            check(name, ok, detail)

        print("\nNegative controls: each check against a planted mistake must fail")
        broken = base.replace("[143,137,", "[143,")
        check("fails on a planted mistake: base.config without exit 137",
              base_config_exits(broken) != run_headroom.MEMORY_EXITS, "passed")
        ms = mutants()
        check("every check has a planted mistake", set(ms) == set(results), set(results) ^ set(ms))
        for name, mod in ms.items():
            res = all_checks(mod)
            check(f"fails on a planted mistake: {name}", not res[name][0], "the check passed on the mutant")
        same = all_checks(mutant())
        check("an unchanged copy passes (the controls test the mistake, not the copy)",
              all(ok for ok, _ in same.values()), [k for k, (ok, _) in same.items() if not ok])
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
