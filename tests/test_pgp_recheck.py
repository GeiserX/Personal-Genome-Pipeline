#!/usr/bin/env python3
"""bin/pgp_recheck.py on synthetic sample folders (tests/fixtures/recheck), with no network:
the upstream answers come from an offline file or from a server on 127.0.0.1.

Checks:
  1. ClinVar: the run used 20260901 and NCBI has 20261005: refresh-available, and
     the steps it changes are 06, every step downstream of it (13 reads ClinVar too,
     20 only for a table of its own, so 12 is not marked) and the report 24;
  2. a ClinVar file installed here that is newer than the one the run used: rerun-ready;
  3. VEP: upstream 117 with the newest pipeline release pinning 117 is
     in-pipeline-release (with that tag); with the release at 116 it is upstream-only;
  4. a newer PGS harmonised file upstream is upstream-only and changes steps 25, 26, 24;
  5. a failed lookup (an offline 500, and a 127.0.0.1 server answering 500 for the
     release and the PGS lookups) is error, never no-newer-release, and the command
     exits 1; the online path prints the one host it asked;
  6. state.json is replaced whole (no temporary file left) and log.tsv only grows;
  7. --rerun refuses (exit 3) while run_status.tsv has no finished_utc;
  8. the comparison after a rerun names a predicted step that was not redone, a
     redone step that was not predicted, and a predicted script step as not redone;
  9. --launch-dir prints the last command of .nextflow/history with -resume;
 10. the live-site plugin (bin/pgp_live_recheck.py) lets the browser fetch
     recheck/state.json, and its state hook gives each step its badges ("newer
     ClinVar available", "redone <date>") and {} when there is no check.
Each check is then run against a planted mistake and must fail there.

Run: python3 tests/test_pgp_recheck.py
"""
import contextlib
import gzip
import http.server
import importlib.util
import io
import json
import os
import shutil
import sys
import tempfile
import threading

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(REPO, "tests", "fixtures", "recheck")
TOOL = os.path.join(REPO, "bin", "pgp_recheck.py")
FAILS = 0
CLINVAR_STEPS = ["06", "13", "20", "23", "24", "30", "31", "38", "39"]
# The pins "here" are the run's own: the version rows of the fixture manifest
# (the image name stays in the fixture, never as a literal in a test).
with open(os.path.join(FIX, "run_manifest.tsv")) as _f:
    THIS_PINS = {r[1]: r[2] for r in (line.rstrip("\n").split("\t") for line in _f) if r[0] == "version"}


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def load():
    spec = importlib.util.spec_from_file_location("pgp_recheck_under_test", TOOL)
    m = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(m)
    finally:
        sys.dont_write_bytecode = keep
    return m


def clinvar_file(path, date):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(f"##fileformat=VCFv4.1\n##fileDate={date}\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def genome(tmp, installed_clinvar="20260901", finished=True):
    g = os.path.join(tmp, "genome")
    s = os.path.join(g, "S1")
    os.makedirs(os.path.join(s, "logs"), exist_ok=True)
    shutil.copy(os.path.join(FIX, "run_manifest.tsv"), os.path.join(s, "run_manifest.tsv"))
    clinvar_file(os.path.join(g, "clinvar", "clinvar_pathogenic_chr.vcf.gz"), installed_clinvar)
    with open(os.path.join(s, "logs", "run_status.tsv"), "w") as f:
        f.write("meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\n"
                "meta\tcommand\techo replayed\nmeta\tswitches\tGENOME_DIR TOOLS\n")
        if finished:
            f.write("meta\tfinished_utc\t2026-09-02 12:00:00 UTC\n")
    return g


def answers(m, path, **over):
    data = json.load(open(os.path.join(FIX, "offline.json")))
    for k, v in over.items():
        data[k.replace("__", ":")] = v
    p = os.path.join(os.path.dirname(path), "answers.json")
    json.dump(data, open(p, "w"))
    return m.Offline(p)


def by_item(state):
    return {i["item"]: i for i in state["items"]}


def all_checks(m):
    out = {}
    tmp = tempfile.mkdtemp(prefix="recheck.", dir=os.path.join(REPO, "tmp") if os.path.isdir(os.path.join(REPO, "tmp")) else None)
    steps = m.run_state.read_steps(os.path.join(REPO, "assets", "steps.tsv"))
    try:
        g = genome(tmp)
        st = m.check("S1", g, answers(m, g), steps, THIS_PINS, "3.66.0")
        it = by_item(st)
        cv = it["clinvar_file_date"]
        out["ClinVar newer upstream is refresh-available, with 06 and its downstream"] = (
            cv["class"] == "refresh-available" and cv["steps"] == CLINVAR_STEPS
            and st["affected_steps"] == CLINVAR_STEPS, (cv, st["affected_steps"]))

        g2 = genome(os.path.join(tmp, "b"), installed_clinvar="20261005")
        st2 = m.check("S1", g2, answers(m, g2), steps, THIS_PINS, "3.66.0")
        out["a newer ClinVar installed here is rerun-ready"] = (
            by_item(st2)["clinvar_file_date"]["class"] == "rerun-ready"
            and st2["rerun_ready_steps"] == CLINVAR_STEPS, by_item(st2)["clinvar_file_date"])

        rel = {"tag": "v9.9.0", "pins": {"VEP_CACHE_RELEASE": "117"}}
        v1 = by_item(m.check("S1", g, answers(m, g, release=rel), steps, THIS_PINS, "3.66.0"))["VEP_CACHE_RELEASE"]
        rel116 = {"tag": "v9.9.0", "pins": {"VEP_CACHE_RELEASE": "116"}}
        v2 = by_item(m.check("S1", g, answers(m, g, release=rel116), steps, THIS_PINS, "3.66.0"))["VEP_CACHE_RELEASE"]
        out["VEP 117 in the newest release is in-pipeline-release, else upstream-only"] = (
            v1["class"] == "in-pipeline-release" and v1["in_pgp_release"] == {"tag": "v9.9.0", "pin": "117"}
            and "v9.9.0" in v1["action"] and v2["class"] == "upstream-only", (v1, v2))

        p = it["pgs:PGS000001"]
        out["a newer PGS file upstream is upstream-only, steps 24 25 26"] = (
            p["class"] == "upstream-only" and p["steps"] == ["24", "25", "26"], p)

        bad = by_item(m.check("S1", g, answers(m, g, PHARMCAT_IMAGE={"error": "HTTP 500"}), steps, THIS_PINS,
                              "3.66.0"))["PHARMCAT_IMAGE"]
        srv, port = serve_500()
        try:
            fr = m.load_freshness()
            fr.RETRIES = 0
            os.environ["FRESHNESS_GITHUB_API"] = f"http://127.0.0.1:{port}"
            err = io.StringIO()
            on = m.Online(fr, {}, "o/r", pgs_api=f"http://127.0.0.1:{port}/rest/score", raw=f"http://127.0.0.1:{port}",
                          out=err)
            st3 = m.check("S1", g, Mixed(on, answers(m, g)), steps, THIS_PINS, "3.66.0")
        finally:
            srv.shutdown()
            os.environ.pop("FRESHNESS_GITHUB_API", None)
        it3 = by_item(st3)
        with contextlib.redirect_stdout(io.StringIO()):
            rc = m.main(["--sample", "S1", "--genome-dir", g, "--offline",
                         write(tmp, "a500.json", dict(json.load(open(os.path.join(FIX, "offline.json"))),
                                                      PHARMCAT_IMAGE={"error": "HTTP 500"}))])
        out["a failed lookup is error and the command exits 1"] = (
            bad["class"] == "error" and "500" in bad["error"] and rc == 1
            and it3["pgs:PGS000001"]["class"] == "error" and it3["VEP_CACHE_RELEASE"]["class"] == "error"
            and "500" in it3["pgs:PGS000001"]["error"] and on.hosts == [f"127.0.0.1:{port}"]
            and f"asking 127.0.0.1:{port}" in err.getvalue(), (bad, rc, it3.get("pgs:PGS000001"), on.hosts))

        d = os.path.join(g, "S1", "recheck")
        with contextlib.redirect_stdout(io.StringIO()):
            m.main(["--sample", "S1", "--genome-dir", g, "--offline", os.path.join(FIX, "offline.json")])
            n1 = sum(1 for _ in open(os.path.join(d, "log.tsv")))
            first = open(os.path.join(d, "log.tsv")).read()
            m.main(["--sample", "S1", "--genome-dir", g, "--offline", os.path.join(FIX, "offline.json")])
        log = open(os.path.join(d, "log.tsv")).read()
        state = json.load(open(os.path.join(d, "state.json")))
        out["state.json is replaced whole and log.tsv only grows"] = (
            log.startswith(first) and log.count("\n") > n1 and sorted(os.listdir(d)) == ["log.tsv", "state.json"]
            and state["sample"] == "S1" and log.startswith("utc\tevent\t"), (n1, log.count("\n"), os.listdir(d)))

        g4 = genome(os.path.join(tmp, "c"), finished=False)
        e = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(e):
            rc = m.main(["--sample", "S1", "--genome-dir", g4, "--offline", os.path.join(FIX, "offline.json"), "--rerun"])
        out["--rerun refuses while a run is going on"] = (rc == 3 and "going on" in e.getvalue(), (rc, e.getvalue()))

        rs = {"steps": [
            {"step": "06", "processes": {"CLINVAR_SCREEN": {"done": 1, "cached": 0}}},
            {"step": "13", "processes": {"VEP": {"done": 0, "cached": 1}}},
            {"step": "07", "processes": {"PHARMCAT": {"done": 1, "cached": 0}}},
            {"step": "14", "processes": {}},
        ]}
        cmp = m.compare_rerun(["06", "13", "14"], rs, steps)
        out["the rerun comparison names each difference"] = (
            cmp["redone"] == ["06", "07"] and cmp["predicted_not_redone"] == ["13"]
            and cmp["redone_not_predicted"] == ["07"] and cmp["not_redone_script_steps"] == ["14"], cmp)

        ld = os.path.join(tmp, "launch")
        os.makedirs(os.path.join(ld, ".nextflow"))
        with open(os.path.join(ld, ".nextflow", "history"), "w") as f:
            f.write("2026-09-01 10:00:00\t1h\tfirst_run\tOK\tabc\tsess1\tnextflow run main.nf --input a.csv\n"
                    "2026-09-02 10:00:00\t2h\tlast_run\tOK\tabc\tsess2\tnextflow run main.nf --input b.csv\n")
        out["--launch-dir gives the last command with -resume"] = (
            m.launch_dir_command(ld) == "nextflow run main.nf --input b.csv -resume", m.launch_dir_command(ld))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


class Mixed:
    """The online release and PGS lookups, the offline answers for the rest."""

    def __init__(self, on, off):
        self.on, self.off, self.hosts = on, off, on.hosts

    def release(self):
        return self.on.release()

    def pgs(self, pid):
        return self.on.pgs(pid)

    def clinvar(self):
        return self.off.clinvar()

    def pinned(self, *a):
        return self.off.pinned(*a)


def write(tmp, name, obj):
    p = os.path.join(tmp, name)
    json.dump(obj, open(p, "w"))
    return p


def serve_500():
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(500)
            self.end_headers()
        do_HEAD = do_GET

        def log_message(self, *a):
            pass
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, srv.server_address[1]


def mutants():
    ms = {}
    m = load()
    orig = m.affected_steps
    m.affected_steps = lambda steps, keys: orig([dict(s, reads=s["reads"].replace("clinvar_file_date,", ""))
                                                 for s in steps], keys)
    ms["ClinVar newer upstream is refresh-available, with 06 and its downstream"] = m
    m = load()
    m.installed_clinvar = lambda g: None
    ms["a newer ClinVar installed here is rerun-ready"] = m
    m = load()
    oc = m.classify

    def no_release(item, tag, pins, err, clinvar=False):
        return oc(item, tag, {}, err, clinvar)
    m.classify = no_release
    ms["VEP 117 in the newest release is in-pipeline-release, else upstream-only"] = m
    m = load()
    m.affected_steps = lambda steps, keys: [s["step"] for s in steps if any(m.key_matches(k, r) for k in keys
                                                                             for r in m.reads_of(s))]
    ms["a newer PGS file upstream is upstream-only, steps 24 25 26"] = m
    m = load()
    oc2 = m.classify
    m.classify = lambda item, *a, **k: "no-newer-release" if item.get("error") else oc2(item, *a, **k)
    ms["a failed lookup is error and the command exits 1"] = m
    m = load()

    def truncating(path, rows):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(m.LOG_HEADER)
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
    m.append_log = truncating
    ms["state.json is replaced whole and log.tsv only grows"] = m
    m = load()
    orc = m.run_state.collect

    def never_running(*a, **k):
        st = orc(*a, **k)
        st["run"]["state"] = "finished"
        return st
    m.run_state = type("RS", (), {"collect": staticmethod(never_running), "read_steps": staticmethod(m.run_state.read_steps),
                                  "read_run_status": staticmethod(m.run_state.read_run_status)})
    m.subprocess = type("SP", (), {"run": staticmethod(lambda *a, **k: type("R", (), {"returncode": 0})())})
    ms["--rerun refuses while a run is going on"] = m
    m = load()
    m.compare_rerun = lambda predicted, rs, steps: {"predicted": predicted, "redone": m.redone_steps(rs),
                                                    "predicted_not_redone": [], "redone_not_predicted": [],
                                                    "not_redone_script_steps": []}
    ms["the rerun comparison names each difference"] = m
    m = load()
    m.launch_dir_command = lambda d: "nextflow run main.nf --input b.csv"
    ms["--launch-dir gives the last command with -resume"] = m
    return ms


def load_plugin():
    spec = importlib.util.spec_from_file_location("pgp_live_recheck_test", os.path.join(REPO, "bin", "pgp_live_recheck.py"))
    p = importlib.util.module_from_spec(spec)
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        spec.loader.exec_module(p)
    finally:
        sys.dont_write_bytecode = keep
    return p


def plugin_ok(p):
    class App:
        def __init__(self):
            self.globs, self.hooks = [], {}

        def run_file(self, g):
            self.globs.append(g)

        def state_hook(self, name, fn):
            self.hooks[name] = fn
    app = App()
    p.register(app)
    tmp = tempfile.mkdtemp(prefix="recheck-live.")
    try:
        run = os.path.join(tmp, "S1")
        none = app.hooks["recheck"](run, {}, [])
        os.makedirs(os.path.join(run, "recheck"))
        json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar (step 06's file)", "class": "refresh-available", "steps": ["06", "24"]},
            {"item": "VEP_CACHE_RELEASE", "label": "Ensembl / VEP cache", "class": "no-newer-release", "steps": ["13"]}],
            "rerun": {"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": {"redone": ["06"]}}},
            open(os.path.join(run, "recheck", "state.json"), "w"))
        got = app.hooks["recheck"](run, {}, [])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = (app.globs == ["recheck/state.json"] and none == {}
          and got.get("steps") == {"06": ["newer ClinVar available", "redone 2026-10-05"], "24": ["newer ClinVar available"]}
          and got.get("banner") == "Database re-check of 2026-10-06T00:00:00Z: newer ClinVar available")
    return ok, (app.globs, none, got)


def main():
    ok, detail = plugin_ok(load_plugin())
    check("the live-site plugin serves the check and gives each step its badges", ok, detail)
    bad = load_plugin()
    bad.ACT = dict(bad.ACT, **{"refresh-available": "newer {name}"})
    check("fails on a planted mistake: the live-site plugin serves the check and gives each step its badges",
          not plugin_ok(bad)[0], "passed on the mutant")
    results = all_checks(load())
    for name, (ok, detail) in results.items():
        check(name, ok, detail)
    print("\nNegative controls: each check against a planted mistake must fail")
    for name, mod in mutants().items():
        check(f"fails on a planted mistake: {name}", not all_checks(mod)[name][0], "the check passed on the mutant")
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
