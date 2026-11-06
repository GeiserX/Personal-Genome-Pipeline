#!/usr/bin/env python3
"""bin/owngenome_recheck.py on synthetic sample folders (tests/fixtures/recheck), with no network:
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
 10. the live-site plugin (bin/owngenome_live_recheck.py) lets the browser fetch
     recheck/state.json, and its state hook gives each step its badges ("newer
     ClinVar available", "redone <date>") and {} when there is no check.
 11. a bumped pin of the opt-in steps 50 to 53 (LOCITYPER_PANEL_VERSION and
     LOCITYPER_IMAGE for 53, MOCHA_PANEL_VERSION and MOCHA_IMAGE for 51,
     VNTYPER_IMAGE for 52, MEGANE_IMAGE for 50) marks its step.
 12. check() itself names step 53 rerun-ready when the run used an older
     LOCITYPER_PANEL_VERSION than this checkout pins, and not when they match;
     its action names setup.sh --locityper-data, which the main setup leaves out.
 13. an image re-pinned at the same tag with a new digest is rerun-ready, for a pin
     with no upstream lookup (step 50's) and for a looked-up one (PharmCAT); the
     same reference on both sides is no-newer-release;
 14. a pin with no upstream lookup that the run's manifest lacks is unknown, never
     no-newer-release, names no step and does not make the command exit 1;
 15. a pin read only by a step the run did not run is left out: an opt-in step not in
     TOOLS (run-all.sh writes no row for it) and one recorded as skipped (no BAM);
     a status file with no step rows at all leaves every pin in;
 16. a malformed Last-Modified of a PGS file, and an offline release with no tag,
     are errors of that lookup, never a crash;
 17. --refresh-clinvar refuses (exit 3) while a run holds the ClinVar lock shared,
     and holds it exclusive while setup.sh refreshes;
 18. the ClinVar lock block of run-all.sh and of step 06 waits while a refresh holds
     the lock and goes on once it ends (needs flock; it fails in CI without it).
Each check is then run against a planted mistake and must fail there.

Run: python3 tests/test_owngenome_recheck.py
"""
import contextlib
import fcntl
import gzip
import http.server
import importlib.util
import io
import json
import os
import re
import shlex
import shutil
import sys
import tempfile
import threading

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(REPO, "tests", "fixtures", "recheck")
TOOL = os.path.join(REPO, "bin", "owngenome_recheck.py")
FAILS = 0
CLINVAR_STEPS = ["06", "13", "20", "23", "24", "27", "30", "31", "38", "39"]
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
    spec = importlib.util.spec_from_file_location("owngenome_recheck_under_test", TOOL)
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


def genome(tmp, installed_clinvar="20260901", finished=True, installed_pgs=None, run_clinvar=None,
           command="echo replayed"):
    g = os.path.join(tmp, "genome")
    s = os.path.join(g, "S1")
    os.makedirs(os.path.join(s, "logs"), exist_ok=True)
    shutil.copy(os.path.join(FIX, "run_manifest.tsv"), os.path.join(s, "run_manifest.tsv"))
    if run_clinvar is not None:  # the run's ClinVar row as given ("" drops it)
        rows = [x for x in open(os.path.join(s, "run_manifest.tsv")) if not x.startswith("data\tclinvar_file_date\t")]
        if run_clinvar:
            rows.append(f"data\tclinvar_file_date\t{run_clinvar} (clinvar/clinvar_pathogenic_chr.vcf.gz)\n")
        open(os.path.join(s, "run_manifest.tsv"), "w").writelines(rows)
    clinvar_file(os.path.join(g, "clinvar", "clinvar_pathogenic_chr.vcf.gz"), installed_clinvar)
    if installed_pgs:
        os.makedirs(os.path.join(g, "prs_scores"), exist_ok=True)
        with gzip.open(os.path.join(g, "prs_scores", "PGS000001.txt.gz"), "wt") as f:
            f.write(f"#HmPOS_date={installed_pgs}\n")
    with open(os.path.join(s, "logs", "run_status.tsv"), "w") as f:
        f.write("meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\n"
                f"meta\tcommand\t{command}\nmeta\tswitches\tGENOME_DIR TOOLS\n")
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
            v1["class"] == "in-pipeline-release" and v1["in_owngenome_release"] == {"tag": "v9.9.0", "pin": "117"}
            and "v9.9.0" in v1["action"] and v2["class"] == "upstream-only", (v1, v2))

        p = it["pgs:PGS000001"]
        out["a newer PGS file upstream is upstream-only, steps 24 25 26"] = (
            p["class"] == "upstream-only" and p["steps"] == ["24", "25", "26"], p)

        gp = genome(os.path.join(tmp, "pgs"), installed_pgs="2026-09-01")
        stp = m.check("S1", gp, answers(m, gp), steps, THIS_PINS, "3.66.0")
        out["a newer PGS file installed here is rerun-ready, steps 24 25 26"] = (
            by_item(stp)["pgs:PGS000001"]["class"] == "rerun-ready"
            and stp["rerun_ready_steps"] == ["24", "25", "26"], by_item(stp)["pgs:PGS000001"])

        try:
            gi = genome(os.path.join(tmp, "impossible"), run_clinvar="20260230")
            ci = by_item(m.check("S1", gi, answers(m, gi), steps, THIS_PINS, "3.66.0"))["clinvar_file_date"]
            out["an impossible ClinVar date of the run is error, not a crash"] = (
                ci["class"] == "error" and "which ClinVar file" in ci["error"], ci)
        except ValueError as e:
            out["an impossible ClinVar date of the run is error, not a crash"] = (False, e)
        gn = genome(os.path.join(tmp, "nodate"), run_clinvar="")
        cn = by_item(m.check("S1", gn, answers(m, gn), steps, THIS_PINS, "3.66.0"))["clinvar_file_date"]
        out["a run whose manifest has no ClinVar date is error"] = (cn["class"] == "error", cn)

        ph = by_item(m.check("S1", g, answers(m, g, PHARMCAT_IMAGE={"upstream": "latest"}), steps, THIS_PINS,
                             "3.66.0"))["PHARMCAT_IMAGE"]
        out["an answer that names no release number is error"] = (
            ph["class"] == "error" and "names no release" in ph["error"], ph)

        hl = by_item(m.check("S1", g, answers(m, g), steps, THIS_PINS, None))["hla_database"]
        out["an unread HLA_DB_RELEASE in this checkout is error"] = (
            hl["class"] == "error" and "HLA_DB_RELEASE" in hl["error"], hl)

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

        # --refresh-clinvar: another sample of the same GENOME_DIR is running, and it reads the same files.
        gs = genome(os.path.join(tmp, "shared"))
        shutil.copytree(os.path.join(gs, "S1", "logs"), os.path.join(gs, "S2", "logs"))
        with open(os.path.join(gs, "S2", "logs", "run_status.tsv"), "w") as f:
            f.write("meta\tstarted_epoch\t1000\nmeta\tstarted_utc\t2026-09-02 10:00:00 UTC\n")
        keep_sp, ran = m.subprocess, []
        m.subprocess = type("SP", (), {"run": staticmethod(lambda *a, **k: ran.append(a) or type("R", (), {"returncode": 0})())})
        e = io.StringIO()
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(e):
                rc = m.main(["--sample", "S1", "--genome-dir", gs, "--offline", os.path.join(FIX, "offline.json"),
                             "--rerun", "--refresh-clinvar"])
        finally:
            m.subprocess = keep_sp
        out["--refresh-clinvar refuses while another sample of GENOME_DIR runs"] = (
            rc == 3 and "S2" in e.getvalue() and not ran, (rc, e.getvalue(), ran))

        gl = genome(os.path.join(tmp, "lock"))
        os.makedirs(os.path.join(gl, "S1", "recheck"), exist_ok=True)
        e = io.StringIO()
        with open(os.path.join(gl, "S1", "recheck", "rerun.lock"), "a") as held:
            fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(e):
                rc = m.main(["--sample", "S1", "--genome-dir", gl, "--offline", os.path.join(FIX, "offline.json"), "--rerun"])
        out["--rerun refuses while another --rerun of the sample holds the lock"] = (
            rc == 3 and "another --rerun" in e.getvalue(), (rc, e.getvalue()))

        # The replay writes a manifest that used the newer ClinVar: the state after it must say so.
        newman = os.path.join(tmp, "new_manifest.tsv")
        with open(os.path.join(FIX, "run_manifest.tsv")) as f, open(newman, "w") as o:
            for x in f:
                o.write("data\tclinvar_file_date\t20261005 (clinvar/clinvar_pathogenic_chr.vcf.gz)\n"
                        if x.startswith("data\tclinvar_file_date\t") else x)
        gr = genome(os.path.join(tmp, "after"), installed_clinvar="20261005",
                    command="cp " + shlex.quote(newman) + " " + shlex.quote(os.path.join(tmp, "after", "genome", "S1", "run_manifest.tsv")))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = m.main(["--sample", "S1", "--genome-dir", gr, "--offline", os.path.join(FIX, "offline.json"), "--rerun"])
        sr = json.load(open(os.path.join(gr, "S1", "recheck", "state.json")))
        cr = by_item(sr)["clinvar_file_date"]
        out["after --rerun the check is taken again on the manifest the rerun wrote"] = (
            rc == 0 and cr["class"] == "no-newer-release" and cr["ran_with"] == "20261005"
            and (sr.get("rerun") or {}).get("rc") == 0, (rc, cr, sr.get("rerun")))

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

        new_checks(m, tmp, steps, out)

        ld = os.path.join(tmp, "launch")
        os.makedirs(os.path.join(ld, ".nextflow"))
        with open(os.path.join(ld, ".nextflow", "history"), "w") as f:
            f.write("2026-09-01 10:00:00\t1h\tfirst_run\tOK\tabc\tsess1\tnextflow run main.nf --input a.csv\n"
                    "2026-09-02 10:00:00\t2h\tlast_run\tOK\tabc\tsess2\tnextflow run main.nf --input b.csv\n")
        out["--launch-dir gives the last command with -resume"] = (
            m.launch_dir_command(ld) == "nextflow run main.nf --input b.csv -resume", m.launch_dir_command(ld))
        lq = os.path.join(tmp, "launch-q")
        os.makedirs(os.path.join(lq, ".nextflow"))
        with open(os.path.join(lq, ".nextflow", "history"), "w") as f:
            f.write("2026-09-02 10:00:00\t2h\tlast_run\tOK\tabc\tsess2\tnextflow run main.nf --input 'data -resume old.csv'\n")
        out["--launch-dir: an argument holding -resume is not the option"] = (
            m.launch_dir_command(lq) == "nextflow run main.nf --input 'data -resume old.csv' -resume", m.launch_dir_command(lq))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def manifest_add(g, rows):
    with open(os.path.join(g, "S1", "run_manifest.tsv"), "a") as f:
        f.writelines(f"version\t{k}\t{v}\n" for k, v in rows)


def step_rows(g, rows):
    """run_status.tsv step rows as run-all.sh writes them; a step it left out gets no row."""
    with open(os.path.join(g, "S1", "logs", "run_status.tsv"), "a") as f:
        f.writelines(f"step\t{k}\t{v}\n" for k, v in rows)


def new_checks(m, tmp, steps, out):
    # Image references built from the fixture's own, so no image name is a literal here.
    ref = THIS_PINS["PHARMCAT_IMAGE"]
    d1, d2 = ref + "@sha256:" + "a" * 64, ref + "@sha256:" + "b" * 64
    g = genome(os.path.join(tmp, "digest"))
    manifest_add(g, [("MEGANE_IMAGE", d1)])
    new = by_item(m.check("S1", g, answers(m, g), steps, dict(THIS_PINS, MEGANE_IMAGE=d2), "3.66.0"))
    same = by_item(m.check("S1", g, answers(m, g), steps, dict(THIS_PINS, MEGANE_IMAGE=d1), "3.66.0"))
    mg, ms = new.get("MEGANE_IMAGE", {}), same.get("MEGANE_IMAGE", {})
    out["a new digest at the same tag of a pin with no lookup is rerun-ready"] = (
        mg.get("class") == "rerun-ready" and "50" in mg.get("steps", []) and mg.get("installed_ref") == d2
        and ms.get("class") == "no-newer-release", (mg, ms))

    gp = genome(os.path.join(tmp, "digest-pharmcat"))
    rows = [x for x in open(os.path.join(gp, "S1", "run_manifest.tsv")) if not x.startswith("version\tPHARMCAT_IMAGE\t")]
    open(os.path.join(gp, "S1", "run_manifest.tsv"), "w").writelines(rows + [f"version\tPHARMCAT_IMAGE\t{d1}\n"])
    ph = by_item(m.check("S1", gp, answers(m, gp), steps, dict(THIS_PINS, PHARMCAT_IMAGE=d2), "3.66.0"))["PHARMCAT_IMAGE"]
    ph1 = by_item(m.check("S1", gp, answers(m, gp), steps, dict(THIS_PINS, PHARMCAT_IMAGE=d1), "3.66.0"))["PHARMCAT_IMAGE"]
    text = io.StringIO()
    m.print_check({"sample": "S1", "checked_utc": "x", "items": [ph], "affected_steps": []}, text)
    out["a new digest at the same tag of a looked-up image is rerun-ready"] = (
        ph["class"] == "rerun-ready" and ph1["class"] == "no-newer-release" and "same tag, new image" in text.getvalue(),
        (ph, ph1, text.getvalue()))

    gu = genome(os.path.join(tmp, "unknown"))
    st = m.check("S1", gu, answers(m, gu), steps, dict(THIS_PINS, MEGANE_IMAGE=d1), "3.66.0")
    mu = by_item(st).get("MEGANE_IMAGE", {})
    with contextlib.redirect_stdout(io.StringIO()):
        ans = write(tmp, "unknown.json", json.load(open(os.path.join(FIX, "offline.json"))))
        keep = m.check
        m.check = lambda s_, g_, a_, st_, pins, hla: keep(s_, g_, a_, st_, dict(pins, MEGANE_IMAGE=d1), hla)
        try:
            rc = m.main(["--sample", "S1", "--genome-dir", gu, "--offline", ans])
        finally:
            m.check = keep
    out["a pin with no lookup the run's manifest lacks is unknown, exit 0"] = (
        mu.get("class") == "unknown" and "50" not in st["affected_steps"] and "50" not in st["rerun_ready_steps"]
        and "older version" in mu.get("action", "") and rc == 0, (mu, st["affected_steps"], rc))

    # A default run (no TOOLS): run-all.sh writes rows for the steps it ran and none for step 50.
    gd = genome(os.path.join(tmp, "not-in-tools"))
    manifest_add(gd, [("MEGANE_IMAGE", d1)])
    step_rows(gd, [("06", "planned"), ("06", "ok")])
    nd = m.check("S1", gd, answers(m, gd), steps, dict(THIS_PINS, MEGANE_IMAGE=d2), "3.66.0")
    out["a pin read only by an opt-in step not in TOOLS is left out"] = (
        "MEGANE_IMAGE" not in by_item(nd) and "50" not in nd["rerun_ready_steps"], by_item(nd).get("MEGANE_IMAGE"))

    gs = genome(os.path.join(tmp, "skipped"))
    manifest_add(gs, [("MEGANE_IMAGE", d1)])
    step_rows(gs, [("06", "ok"), ("50", "skipped (no BAM)")])
    sk = m.check("S1", gs, answers(m, gs), steps, dict(THIS_PINS, MEGANE_IMAGE=d2), "3.66.0")
    gr = genome(os.path.join(tmp, "ran"))
    manifest_add(gr, [("MEGANE_IMAGE", d1)])
    step_rows(gr, [("06", "ok"), ("50", "planned"), ("50", "ok")])
    rn = by_item(m.check("S1", gr, answers(m, gr), steps, dict(THIS_PINS, MEGANE_IMAGE=d2), "3.66.0"))
    out["a pin read only by a step the run skipped is left out"] = (
        "MEGANE_IMAGE" not in by_item(sk) and "50" not in sk["rerun_ready_steps"]
        and rn.get("MEGANE_IMAGE", {}).get("class") == "rerun-ready", (by_item(sk).get("MEGANE_IMAGE"), rn.get("MEGANE_IMAGE")))

    srv, port = serve_pgs("not a date")
    try:
        fr = m.load_freshness()
        fr.RETRIES = 0
        on = m.Online(fr, {}, "o/r", pgs_api=f"http://127.0.0.1:{port}/rest/score", raw=f"http://127.0.0.1:{port}",
                      out=io.StringIO())
        try:
            gm = genome(os.path.join(tmp, "lastmod"))
            pg = by_item(m.check("S1", gm, Mixed(on, answers(m, gm)), steps, THIS_PINS, "3.66.0"))["pgs:PGS000001"]
            got = (pg["class"] == "error" and "Last-Modified" in pg.get("error", ""), pg)
        except Exception as e:  # noqa: BLE001  the crash this check is about
            got = (False, repr(e))
    finally:
        srv.shutdown()
    out["a malformed Last-Modified of a PGS file is error, not a crash"] = got

    try:
        go = genome(os.path.join(tmp, "notag"))
        so = m.check("S1", go, answers(m, go, release={"pins": {}}), steps, THIS_PINS, "3.66.0")
        got = ("error" in so["pipeline_release"] and "tag" in so["pipeline_release"]["error"], so["pipeline_release"])
    except Exception as e:  # noqa: BLE001
        got = (False, repr(e))
    out["an offline release with no tag is error, not a crash"] = got

    gl = genome(os.path.join(tmp, "cvlock"))
    lock = os.path.join(gl, "clinvar", ".refresh.lock")
    keep_sp, ran = m.subprocess, []
    m.subprocess = type("SP", (), {"run": staticmethod(lambda *a, **k: ran.append(a) or type("R", (), {"returncode": 0})())})
    e = io.StringIO()
    try:
        with open(lock, "a") as held:
            fcntl.flock(held, fcntl.LOCK_SH | fcntl.LOCK_NB)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(e):
                rc = m.main(["--sample", "S1", "--genome-dir", gl, "--offline", os.path.join(FIX, "offline.json"),
                             "--rerun", "--refresh-clinvar"])
    finally:
        m.subprocess = keep_sp
    out["--refresh-clinvar refuses while a run holds the ClinVar lock"] = (
        rc == 3 and ".refresh.lock" in e.getvalue() and not ran, (rc, e.getvalue(), ran))

    seen = []

    def setup_run(argv, *a, **k):
        if "--refresh" in argv:
            with open(lock, "a") as f:
                try:
                    fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
                    seen.append("free")
                except OSError:
                    seen.append("held")
        return type("R", (), {"returncode": 0})()
    m.subprocess = type("SP", (), {"run": staticmethod(setup_run)})
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = m.main(["--sample", "S1", "--genome-dir", gl, "--offline", os.path.join(FIX, "offline.json"),
                         "--rerun", "--refresh-clinvar"])
        with open(lock, "a") as f:  # let go before the replay's run-all.sh takes it shared
            fcntl.flock(f, fcntl.LOCK_SH | fcntl.LOCK_NB)
            after = True
    except OSError:
        after = False
    finally:
        m.subprocess = keep_sp
    out["--refresh-clinvar holds the ClinVar lock while setup.sh refreshes"] = (
        seen == ["held"] and after and rc == 0, (seen, after, rc))


def lock_block(path):
    """The ClinVar lock block of a script: from 'CVLOCK=' through its closing 'fi'."""
    x = re.search(r"^CVLOCK=.*?\n  fi\nfi\n", open(path).read(), re.S | re.M)
    return x.group(0) if x else ""


def lock_waits(block):
    """block, run in bash while this process holds the lock exclusive: it must wait, then go on
    once the lock is let go."""
    import subprocess
    import time
    tmp = tempfile.mkdtemp(prefix="recheck-cvlock.")
    try:
        g = os.path.join(tmp, "genome")
        os.makedirs(os.path.join(g, "clinvar"))
        lock = os.path.join(g, "clinvar", ".refresh.lock")
        with open(lock, "a") as held:
            fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
            p = subprocess.Popen(["bash", "-c", "set -euo pipefail\n" + block + "echo through\n"],
                                 env=dict(os.environ, G=g, GENOME_DIR=g), stdout=subprocess.PIPE, text=True)
            time.sleep(1.5)
            waited = p.poll() is None
        try:
            o, _ = p.communicate(timeout=20)
        except subprocess.TimeoutExpired:
            p.kill()
            o, _ = p.communicate()
        return (bool(block) and waited and p.returncode == 0 and "Waiting" in o and o.rstrip().endswith("through"),
                (waited, p.returncode, o))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def serve_pgs(last_modified):
    """A PGS Catalog score answer, and a HEAD of its file with this Last-Modified."""
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            port = self.server.server_address[1]
            body = json.dumps({"id": "PGS000001", "ftp_harmonized_scoring_files": {"GRCh38": {
                "positions": f"http://127.0.0.1:{port}/f/PGS000001_hmPOS_GRCh38.txt.gz"}}}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_HEAD(self):
            self.send_response(200)
            self.send_header("Last-Modified", last_modified)
            self.send_header("Content-Length", "0")
            self.end_headers()

        def log_message(self, *a):
            pass
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, srv.server_address[1]


def from_source(old, new):
    """bin/owngenome_recheck.py with old replaced by new (which must be in it), as a module."""
    import types
    src = open(TOOL).read()
    if src.count(old) != 1:
        raise AssertionError(f"planted mistake not in the copy: {old[:60]!r}")
    mod = types.ModuleType("owngenome_recheck_planted")
    mod.__file__ = TOOL
    keep, sys.dont_write_bytecode = sys.dont_write_bytecode, True
    try:
        exec(compile(src.replace(old, new), TOOL, "exec"), mod.__dict__)
    finally:
        sys.dont_write_bytecode = keep
    return mod


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
    m.running_samples = lambda genome_dir: []
    ms["--refresh-clinvar refuses while another sample of GENOME_DIR runs"] = m
    m = load()
    m.compare_rerun = lambda predicted, rs, steps: {"predicted": predicted, "redone": m.redone_steps(rs),
                                                    "predicted_not_redone": [], "redone_not_predicted": [],
                                                    "not_redone_script_steps": []}
    ms["the rerun comparison names each difference"] = m
    m = load()
    m.launch_dir_command = lambda d: "nextflow run main.nf --input b.csv"
    ms["--launch-dir gives the last command with -resume"] = m
    m = load()
    m.header_value = lambda *a: None
    ms["a newer PGS file installed here is rerun-ready, steps 24 25 26"] = m
    m = load()

    def date8_unchecked(text):
        x = re.search(r"(\d{4})-?(\d{2})-?(\d{2})", str(text or ""))
        return "".join(x.groups()) if x else None
    m.date8 = date8_unchecked
    ms["an impossible ClinVar date of the run is error, not a crash"] = m
    m = load()
    oc3 = m.classify

    def installed_as_ran(item, *a, **k):
        item.setdefault("ran_with", None)
        item["ran_with"] = item["ran_with"] or item.get("installed")
        return oc3(item, *a, **k)
    m.classify = installed_as_ran
    ms["a run whose manifest has no ClinVar date is error"] = m
    m = load()
    ov = m.ver
    m.ver = lambda text: ov(text) or (0,)
    ms["an answer that names no release number is error"] = m
    m = load()
    ock = m.check
    m.check = lambda sample, g, a, steps, pins, hla: ock(sample, g, a, steps, pins, hla or "3.66.0")
    ms["an unread HLA_DB_RELEASE in this checkout is error"] = m
    m = load()
    m.fcntl = type("F", (), {"LOCK_EX": 0, "LOCK_NB": 0, "flock": staticmethod(lambda *a: None)})
    m.subprocess = type("SP", (), {"run": staticmethod(lambda *a, **k: type("R", (), {"returncode": 0})())})
    ms["--rerun refuses while another --rerun of the sample holds the lock"] = m
    m = load()
    orm, seen = m.read_manifest, {}
    m.read_manifest = lambda path: seen.setdefault(path, orm(path))  # the manifest as first read: no check after
    ms["after --rerun the check is taken again on the manifest the rerun wrote"] = m
    m = load()

    def regex_resume(d):
        cmd = open(os.path.join(d, ".nextflow", "history")).read().splitlines()[-1].split("\t")[-1].strip()
        return cmd if re.search(r"(^|\s)-resume(\s|$)", cmd) else cmd + " -resume"
    m.launch_dir_command = regex_resume
    ms["--launch-dir: an argument holding -resume is not the option"] = m
    ms["a new digest at the same tag of a pin with no lookup is rerun-ready"] = from_source(
        '"rerun-ready" if ran_ref != inst_ref else', '"rerun-ready" if ran != inst else')
    ms["a new digest at the same tag of a looked-up image is rerun-ready"] = from_source(
        "and ran_ref != inst_ref:\n        return \"rerun-ready\"", "and False:\n        return \"rerun-ready\"")
    ms["a pin with no lookup the run's manifest lacks is unknown, exit 0"] = from_source(
        'cls = "unknown" if not', 'cls = "no-newer-release" if not')
    ms["a pin read only by an opt-in step not in TOOLS is left out"] = from_source(
        'if s["state"] in ("skipped", "not planned")}', 'if s["state"] in ("skipped",)}')
    m = load()
    m.steps_not_run = lambda *a: set()
    ms["a pin read only by a step the run skipped is left out"] = m
    ms["a malformed Last-Modified of a PGS file is error, not a crash"] = from_source(
        "        try:  # a malformed date is a failed lookup of this score, never a crash of the whole check\n"
        "            return email.utils.parsedate_to_datetime(lm).strftime(\"%Y%m%d\")\n"
        "        except (TypeError, ValueError, AttributeError, OverflowError):\n",
        "        return email.utils.parsedate_to_datetime(lm).strftime(\"%Y%m%d\")\n        if False:\n")
    ms["an offline release with no tag is error, not a crash"] = from_source(
        'tag, pins = e.get("tag"), e.get("pins") or {}\n        if not isinstance(tag, str) or not tag or not isinstance(pins, dict):',
        'tag, pins = e["tag"], e.get("pins") or {}\n        if not isinstance(pins, dict):')
    m = load()
    m.fcntl = type("F", (), {"LOCK_EX": 0, "LOCK_NB": 0, "flock": staticmethod(lambda *a: None)})
    ms["--refresh-clinvar refuses while a run holds the ClinVar lock"] = m
    ms["--refresh-clinvar holds the ClinVar lock while setup.sh refreshes"] = from_source(
        "            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)\n        except OSError:\n"
        "            print(f\"ERROR: a run or step 06",
        "            pass\n        except OSError:\n            print(f\"ERROR: a run or step 06")
    return ms


def load_plugin():
    spec = importlib.util.spec_from_file_location("owngenome_live_recheck_test", os.path.join(REPO, "bin", "owngenome_live_recheck.py"))
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


def plugin_odd_rerun_ok(p):
    """A rerun block that is not an object: the banner still comes, no badge from it."""
    tmp = tempfile.mkdtemp(prefix="recheck-live-odd.")
    try:
        run = os.path.join(tmp, "S1")
        os.makedirs(os.path.join(run, "recheck"))
        out = []
        for rr in ([1], {"finished_utc": "2026-10-05T10:00:00Z", "predicted_vs_redone": ["06"]}):
            json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [], "rerun": rr},
                      open(os.path.join(run, "recheck", "state.json"), "w"))
            try:
                out.append(p.recheck_state(run, {}, []))
            except (AttributeError, TypeError) as e:
                out.append(e)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = all(isinstance(g, dict) and g.get("steps") == {} and "nothing newer" in g.get("banner", "") for g in out)
    return ok, out


def plugin_stale_ok(p):
    """A check older than the run's run_manifest.tsv describes an earlier run: no banner, no badge.
    One written after the manifest is still shown."""
    tmp = tempfile.mkdtemp(prefix="recheck-live-stale.")
    try:
        run = os.path.join(tmp, "S1")
        os.makedirs(os.path.join(run, "recheck"))
        json.dump({"checked_utc": "2026-10-06T00:00:00Z", "items": [
            {"item": "clinvar_file_date", "label": "ClinVar", "class": "refresh-available", "steps": ["06"]}]},
            open(os.path.join(run, "recheck", "state.json"), "w"))
        out = []
        for written in ("2026-10-07 09:00:00 UTC", "2026-10-05 09:00:00 UTC"):
            with open(os.path.join(run, "run_manifest.tsv"), "w") as f:
                f.write(f"run\twritten_utc\t{written}\n")
            out.append(p.recheck_state(run, {}, []))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    ok = out[0] == {} and out[1].get("steps") == {"06": ["newer ClinVar available"]}
    return ok, out


def plugin_from_source(src, name):
    p = os.path.join(REPO, "tmp") if os.path.isdir(os.path.join(REPO, "tmp")) else tempfile.gettempdir()
    path = os.path.join(tempfile.mkdtemp(prefix="recheck-plugin.", dir=p), name + ".py")
    open(path, "w").write(src)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Steps 50-53 read pins of versions.env that change what they find: the data
# setup.sh builds or downloads with them, or the reference data inside the
# image. A bump of each must mark its step (steps.tsv's reads column).
OPT_IN_PINS = [("MEGANE_IMAGE", "50"), ("MOCHA_PANEL_VERSION", "51"), ("MOCHA_IMAGE", "51"),
               ("VNTYPER_IMAGE", "52"), ("LOCITYPER_PANEL_VERSION", "53"), ("LOCITYPER_IMAGE", "53")]


def pins_ok(m, steps):
    """Each pin of OPT_IN_PINS, bumped alone, makes affected_steps list exactly its step and the
    HTML report (no step reads the outputs of steps 50 to 53)."""
    missed = [(k, n, got) for k, n in OPT_IN_PINS
              for got in [m.affected_steps(steps, [k])] if sorted(got) != sorted([n, m.REPORT_STEP])]
    return not missed, missed


def check_marks_53(m, steps):
    """check() on a run that used an older LOCITYPER_PANEL_VERSION than this checkout pins puts 53 in
    rerun_ready_steps; the same pin on both sides puts nothing there."""
    tmp = tempfile.mkdtemp(prefix="recheck.", dir=os.path.join(REPO, "tmp") if os.path.isdir(os.path.join(REPO, "tmp")) else None)
    try:
        g = genome(tmp, installed_clinvar="20260901")
        with open(os.path.join(g, "S1", "run_manifest.tsv"), "a") as f:
            f.write("version\tLOCITYPER_PANEL_VERSION\tpanel-1.0\n")
        st = m.check("S1", g, answers(m, g), steps, dict(THIS_PINS, LOCITYPER_PANEL_VERSION="panel-2.0"), "3.66.0")
        same = m.check("S1", g, answers(m, g), steps, dict(THIS_PINS, LOCITYPER_PANEL_VERSION="panel-1.0"), "3.66.0")
        act = by_item(st).get("LOCITYPER_PANEL_VERSION", {}).get("action", "")
        return ("53" in st["rerun_ready_steps"] and "53" not in same["rerun_ready_steps"]
                and "./scripts/setup.sh --locityper-data $GENOME_DIR" in act,
                (st["rerun_ready_steps"], same["rerun_ready_steps"], act))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    m = load()
    steps = m.run_state.read_steps(m.STEPS_TSV)
    ok, detail = check_marks_53(m, steps)
    check("a run on an older LOCITYPER_PANEL_VERSION makes check() name step 53 rerun-ready, with setup.sh --locityper-data", ok, detail)
    blank = [dict(st, reads="") if st["step"] == "53" else st for st in steps]
    check("fails on a planted mistake: a run on an older LOCITYPER_PANEL_VERSION makes check() name step 53 "
          "(row 53 with no reads)", not check_marks_53(m, blank)[0], "passed on the mutant")
    ok, detail = pins_ok(m, steps)
    check("a bumped pin of steps 50-53 (LOCITYPER_PANEL_VERSION and the rest) marks its step", ok, detail)
    blank = [dict(st, reads="") if st["step"] == "53" else st for st in steps]
    check("fails on a planted mistake: a bumped pin of steps 50-53 marks its step (row 53 with no reads)",
          not pins_ok(m, blank)[0], "passed on the mutant")
    no_report = [st for st in steps if st["step"] != m.REPORT_STEP]
    check("fails on a planted mistake: a bumped pin of steps 50-53 marks its step (no report step in the catalogue)",
          not pins_ok(m, no_report)[0], "passed on the mutant")
    ok, detail = plugin_ok(load_plugin())
    check("the live-site plugin serves the check and gives each step its badges", ok, detail)
    bad = load_plugin()
    bad.ACT = dict(bad.ACT, **{"refresh-available": "newer {name}"})
    check("fails on a planted mistake: the live-site plugin serves the check and gives each step its badges",
          not plugin_ok(bad)[0], "passed on the mutant")
    ok, detail = plugin_odd_rerun_ok(load_plugin())
    check("the live-site plugin survives a rerun block that is not an object", ok, detail)
    src = open(os.path.join(REPO, "bin", "owngenome_live_recheck.py")).read()
    guard = "rr = rr if isinstance(rr, dict) else {}"
    check("planted mistake is in the copy: the live-site plugin survives a rerun block that is not an object", guard in src)
    unguarded = plugin_from_source(src.replace(guard, "rr = rr or {}"), "owngenome_live_recheck_unguarded")
    check("fails on a planted mistake: the live-site plugin survives a rerun block that is not an object",
          not plugin_odd_rerun_ok(unguarded)[0], "passed on the mutant")
    ok, detail = plugin_stale_ok(load_plugin())
    check("the live-site plugin leaves out a check older than the run's manifest", ok, detail)
    stale = "if utc14(manifest_written(run_dir)) > utc14(st.get(\"checked_utc\")):"
    check("planted mistake is in the copy: the live-site plugin leaves out a check older than the run's manifest",
          stale in src)
    unfiltered = plugin_from_source(src.replace(stale, "if False:"), "owngenome_live_recheck_unfiltered")
    check("fails on a planted mistake: the live-site plugin leaves out a check older than the run's manifest",
          not plugin_stale_ok(unfiltered)[0], "passed on the mutant")
    import shutil as _sh
    for script in ("scripts/run-all.sh", "scripts/06-clinvar-screen.sh"):
        block = lock_block(os.path.join(REPO, script))
        if not _sh.which("flock"):
            check(f"the ClinVar lock block of {script} waits for a refresh (flock is not installed: "
                  f"{'FAIL in CI' if os.environ.get('CI') else 'not run here'})", not os.environ.get("CI"))
            continue
        ok, detail = lock_waits(block)
        check(f"the ClinVar lock block of {script} waits for a refresh, then goes on", ok, detail)
        planted = re.sub(r"(?m)^    flock -s -n 7 .*$", "    :", block)
        check(f"planted mistake is in the copy: the ClinVar lock block of {script}", planted != block)
        check(f"fails on a planted mistake: the ClinVar lock block of {script} waits for a refresh (no flock call)",
              not lock_waits(planted)[0], "passed on the mutant")
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
