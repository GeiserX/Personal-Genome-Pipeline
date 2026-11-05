#!/usr/bin/env python3
"""build_live_demo.py: the live site's public demo, a static folder replaying an invented run.

  build_live_demo.py --out DIR [--genome-dir GDIR] [--records N] [--work WDIR]

writes DIR, which the docs workflow publishes at /demo/ of the docs site:

  index.html                    the list of runs, as the live site's / shows it
  assets/                       assets/live/ as it is: live.js, live.css, plugins/*.js
  runs/DEMO-001/index.html      the live page (assets/live/index.html) in replay mode
  runs/DEMO-001/frames.json     the page's state.json at each moment of the run
  runs/DEMO-001/...             the files the plugins point at (the charts, the
                                search index, the re-check), the report as
                                report.html and the files it links to
  runs/DEMO-FAMILY/             the family page of DEMO-001 and its invented partner
  family/DEMO-FAMILY/           the pedigree and the files the family page opens

No state is written by hand. The frames are bin/owngenome_live.py's own
Live.state(), with every plugin bin/owngenome_live_*.py loaded and each state
hook called, while the messages of tests/fixtures/live/demo-001.weblog.jsonl
(tests/demo/make_demo_events.py) are stored one by one, as the ingest port
stores them. The run folder is DEMO-001, the invented sample of the docs
pictures: GDIR/DEMO-001 when given (after tests/demo/retake-screenshots.sh ran
its steps), else written now by tests/demo/make_demo_sample.py with N records
(default 20000) and its report written by bin/render_report.py with step 24's
arguments. Steps 27, 28 and 36 run in containers and are not run here, so
without GDIR their sections stay missing in the report and on the page.

What the replay changes, and why:
  - step 24's files (the report, the charts, the search index, summary.json)
    are held back until the replay reaches the end of HTML_REPORT, as a run
    writes them then; the re-check (tests/fixtures/live/demo) is written one
    hour after the run ends, and its steps come from owngenome_recheck's own
    affected_steps;
  - the clock the plugins read is the time of the message just stored, so the
    run-health card shows the run as it was then, not as a run gone quiet;
    the memory and the disk it measures are an invented machine's (64 GiB of
    RAM, a 2 TiB disk 55% used), not those of the machine that builds the demo;
  - logs/run_status.tsv and run_manifest.tsv of GDIR's DEMO-001 describe
    another run, so they are left out;
  - the page cannot get a header from a static host, so it carries the
    server's Content-Security-Policy as a meta tag (frame-ancestors, which a
    meta tag cannot set, left out), and the site's absolute links (/runs/...,
    /family/...) become relative ones;
  - the genome viewer opens the reader's own reads on the reader's own
    computer, so the demo's viewer/ page says it is not in the demo.

Standard library only. tests/test_live_demo.py builds it and checks it.
"""
import argparse
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
import owngenome_live  # noqa: E402
import owngenome_recheck  # noqa: E402
import pedigree_svg  # noqa: E402
import run_headroom  # noqa: E402
import run_state  # noqa: E402

SAMPLE = "DEMO-001"
FAMILY = "DEMO-FAMILY"
PARTNER = "DEMO-002"
FIX = os.path.join(REPO, "tests", "fixtures", "live")
WEBLOG = os.path.join(FIX, "demo-001.weblog.jsonl")
EXTRA = os.path.join(FIX, "demo")
ASSETS = owngenome_live.ASSETS
REPORT_PROCESS = "HTML_REPORT"
HELD = (f"{SAMPLE}_report.html", f"{SAMPLE}_charts.json", f"{SAMPLE}_search_index.json", "summary.json")
OTHER_RUN = ("logs/run_status.tsv", "run_manifest.tsv")
RECHECK_AFTER_S = 3600
# The invented machine the run-health card measures: 64 GiB of RAM and a 2 TiB disk, 55% used.
DEMO_RAM = 64 << 30
DEMO_DISK = os.statvfs_result((4096, 4096, 512 << 20, 230 << 20, 230 << 20, 128 << 20, 100 << 20, 100 << 20, 0, 255))
INTERVAL_MS = 600
# The server's policy; a meta tag cannot set frame-ancestors (the browser ignores it there and warns).
CSP = "; ".join(p for p in owngenome_live.CSP.split("; ") if not p.startswith("frame-ancestors"))


def utc(epoch):
    return datetime.fromtimestamp(epoch, tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def epoch(text):
    return datetime.strptime(text, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()


# --- the run folder -------------------------------------------------------------------------

def report(d):
    """bin/render_report.py with the arguments scripts/24-html-report.sh gives it."""
    subprocess.run([sys.executable, os.path.join(BIN, "render_report.py"), "--assets", os.path.join(REPO, "assets"),
                    "--sample", SAMPLE, "--sample-dir", d, "--json", os.path.join(d, "summary.json"),
                    "-o", os.path.join(d, f"{SAMPLE}_report.html"),
                    "--search-index", os.path.join(d, f"{SAMPLE}_search_index.json"),
                    "--charts-json", os.path.join(d, f"{SAMPLE}_charts.json")],
                   check=True, stdout=subprocess.DEVNULL)


def prepare(root, genome_dir=None, records=20000):
    """root/DEMO-001 with its report, the family folder and the held-back files; returns the held folder."""
    d = os.path.join(root, SAMPLE)
    if genome_dir:
        src = os.path.join(genome_dir, SAMPLE)
        if not os.path.isdir(src):
            raise SystemExit(f"ERROR: {src} not found: run tests/demo/retake-screenshots.sh, or leave out --genome-dir")
        shutil.copytree(src, d, symlinks=False)
    else:
        subprocess.run([sys.executable, os.path.join(HERE, "make_demo_sample.py"), "--genome-dir", root,
                        "--records", str(records)], check=True, stdout=subprocess.DEVNULL)
    if not os.path.isfile(os.path.join(d, f"{SAMPLE}_report.html")):
        report(d)
    for rel in OTHER_RUN:
        p = os.path.join(d, rel)
        if os.path.isfile(p):
            os.remove(p)
    held = os.path.join(root, ".held")
    os.makedirs(held)
    for f in HELD:
        if os.path.isfile(os.path.join(d, f)):
            shutil.move(os.path.join(d, f), os.path.join(held, f))
    fam = os.path.join(root, FAMILY)
    shutil.copytree(os.path.join(EXTRA, FAMILY), fam)
    couple = os.path.join(d, "couple", f"{SAMPLE}_couple_carrier.tsv")
    if os.path.isfile(couple):
        a, b = sorted((SAMPLE, PARTNER))
        shutil.copy(couple, os.path.join(fam, f"couple_{a}_{b}.tsv"))
    return held


def write_recheck(root, checked):
    with open(os.path.join(EXTRA, SAMPLE, "recheck", "state.json")) as f:
        st = json.load(f)
    steps = run_state.read_steps()
    st["checked_utc"] = checked
    for i in st["items"]:
        i["steps"] = owngenome_recheck.affected_steps(steps, [i["item"]]) if i["class"] != "no-newer-release" else []
    d = os.path.join(root, SAMPLE, "recheck")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "state.json"), "w") as f:
        json.dump(st, f, indent=1)


# --- the frames -----------------------------------------------------------------------------

def relink(v):
    """The site's absolute links as links relative to a run page (runs/<run>/)."""
    if isinstance(v, dict):
        return {k: relink(x) for k, x in v.items()}
    if isinstance(v, list):
        return [relink(x) for x in v]
    if isinstance(v, str) and re.match(r"^/(?:runs|family)/", v):
        return "../.." + v
    return v


def frame_key(st):
    """A new frame each time a step ends, the run's state changes or a plugin's state does (the
    run-health numbers aside: they move with each message). A frame per message would make
    frames.json tens of megabytes; each frame still shows the steps running at its moment."""
    return json.dumps([st["run"]["state"], [s["step"] for s in st["steps"] if s["state"] in ("done", "failed")],
                       {k: v for k, v in st["plugins"].items() if k != "ops"}], sort_keys=True, default=str)


def replay(root, held, drop_hook=None):
    """(frames of DEMO-001, the family page's state, the app) with the clock at each message."""
    app = owngenome_live.App()
    cfg = owngenome_live.configure(["--root", root, "--state-dir", os.path.join(root, ".owngenome-live"),
                                    "--eta-cores", "16", "--eta-box", "demo", "--quiet"], app)
    if drop_hook:  # tests/test_live_demo.py's planted mistake: one plugin's state left out
        app.state_hooks = [h for h in app.state_hooks if h[0] != drop_hook]
    live = owngenome_live.Live(app, cfg, None)
    with open(WEBLOG) as f:
        events = [owngenome_live.cut(json.loads(line)) for line in f if line.strip()]
    clock = [epoch(events[0]["utcTime"])]
    real = time.time, run_headroom._statvfs, run_headroom.physical_memory
    time.time = lambda: clock[0]
    run_headroom._statvfs = lambda path: DEMO_DISK
    run_headroom.physical_memory = lambda: DEMO_RAM
    frames, last_key = [], None

    def snap(t):
        nonlocal last_key
        clock[0] = epoch(t)
        st = live.state(live.run(SAMPLE))
        st["generated_utc"] = t
        key = frame_key(st)
        if key != last_key:
            last_key = key
            st["plugins"] = relink(st["plugins"])
            frames.append({"t": t, "state": st})
    try:
        for i, ev in enumerate(events):
            live.store.append(SAMPLE, ev)
            tr = ev.get("trace") or {}
            if ev["event"] == "process_completed" and run_state.task_name(tr.get("name", ""))[0] == REPORT_PROCESS:
                for f in os.listdir(held):
                    shutil.move(os.path.join(held, f), os.path.join(root, SAMPLE, f))
            # One look per moment: the messages of one second (a failed attempt and its resubmission)
            # are stored together before the state is read.
            if i == len(events) - 1 or events[i + 1].get("utcTime") != ev["utcTime"]:
                snap(ev["utcTime"])
        checked = utc(epoch(events[-1]["utcTime"]) + RECHECK_AFTER_S)
        write_recheck(root, checked)
        snap(checked)
        fam = live.state(live.run(FAMILY))
        fam["generated_utc"] = checked
        fam["plugins"] = relink(fam["plugins"])
    finally:
        time.time, run_headroom._statvfs, run_headroom.physical_memory = real
    return frames, {"t": checked, "state": fam}, app


# --- the static folder ----------------------------------------------------------------------

def with_csp(page):
    meta = f'<meta http-equiv="Content-Security-Policy" content="{html.escape(CSP, quote=False)}">'
    return page.replace('<meta charset="utf-8">', '<meta charset="utf-8">\n' + meta, 1)


def run_page(name, plugin_js):
    with open(os.path.join(ASSETS, "index.html")) as f:
        page = f.read()
    tags = "\n".join(f'<script src="../../assets/plugins/{html.escape(f)}" defer></script>' for f in plugin_js)
    page = (page.replace("{{assets}}", "../../assets").replace("{{run}}", html.escape(name))
            .replace("{{mode}}", "replay").replace("{{plugins}}", tags)
            .replace('<a href="/">All runs</a>', '<a href="../../">All runs</a>'))
    return with_csp(page)


def plain_page(title, body, css):
    return with_csp('<!doctype html>\n<html lang="en" data-theme="dark"><head><meta charset="utf-8">'
                    '<meta name="viewport" content="width=device-width, initial-scale=1">'
                    f'<title>{html.escape(title)}</title><link rel="stylesheet" href="{css}"></head>'
                    f'<body><header><h1>{html.escape(title)}</h1></header><main>{body}</main></body></html>\n')


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"))


def write_site(out, root, frames, fam):
    if os.path.exists(out):
        shutil.rmtree(out)
    shutil.copytree(ASSETS, os.path.join(out, "assets"), ignore=shutil.ignore_patterns("index.html", "login.html"))
    plugin_dir = os.path.join(ASSETS, "plugins")
    plugin_js = sorted(f for f in os.listdir(plugin_dir) if owngenome_live.ASSET_RE.match(f) and f.endswith(".js")) \
        if os.path.isdir(plugin_dir) else []

    runs = [(SAMPLE, "sample"), (FAMILY, "family")]
    items = "\n".join(f'<li><a href="runs/{n}/">{n}</a> <span class="kind">{k}</span></li>' for n, k in runs)
    with open(os.path.join(out, "index.html"), "w") as f:
        f.write(plain_page("Runs", '<p class="muted">A demo of the live site on DEMO-001, an invented sample: '
                           'the run, the results and the family are made up. Nothing here is a real person.</p>'
                           f'<ul class="runs">{items}</ul>', "assets/live.css"))

    sd, d = os.path.join(out, "runs", SAMPLE), os.path.join(root, SAMPLE)
    os.makedirs(sd)
    with open(os.path.join(sd, "index.html"), "w") as f:
        f.write(run_page(SAMPLE, plugin_js))
    write_json(os.path.join(sd, "frames.json"), {"interval_ms": INTERVAL_MS, "autoplay": True, "frames": frames})
    files = {f"{SAMPLE}_report.html": "report.html"}
    # The files the report links to, as the live server lists them (Live.run_files reads no state).
    files.update({rel: rel for rel in owngenome_live.Live.run_files(None, {"dir": d, "name": SAMPLE})})
    for v in frames[-1]["state"]["plugins"].values():
        if isinstance(v, dict) and isinstance(v.get("file"), str):
            files[v["file"]] = v["file"]
    for src, dst in files.items():
        p = os.path.join(d, src)
        if os.path.isfile(p):
            os.makedirs(os.path.dirname(os.path.join(sd, dst)), exist_ok=True)
            shutil.copy(p, os.path.join(sd, dst))
    os.makedirs(os.path.join(sd, "viewer"))
    with open(os.path.join(sd, "viewer", "index.html"), "w") as f:
        f.write(plain_page("Genome viewer: not in the demo",
                           "<p>On your own live site this link opens the genome viewer on your own reads, on your "
                           "own computer. The demo has no reads to show.</p><p><a href=\"../\">Back to the run</a></p>",
                           "../../../assets/live.css"))

    fd = os.path.join(out, "runs", FAMILY)
    os.makedirs(fd)
    with open(os.path.join(fd, "index.html"), "w") as f:
        f.write(run_page(FAMILY, plugin_js))
    write_json(os.path.join(fd, "frames.json"), {"interval_ms": INTERVAL_MS, "frames": [fam]})
    famdir, ffiles = os.path.join(out, "family", FAMILY), os.path.join(root, FAMILY)
    os.makedirs(famdir)
    order, members, _ = pedigree_svg.load(ffiles)
    with open(os.path.join(famdir, "pedigree.svg"), "w") as f:
        f.write(pedigree_svg.draw(order, members, FAMILY))
    for name in (fam["state"]["plugins"].get("family") or {}).get("files") or []:
        shutil.copy(os.path.join(ffiles, name), os.path.join(famdir, name))


def build(out, genome_dir=None, records=20000, work=None, drop_hook=None):
    """Write the demo into out; returns {frames, hooks, family} for the tests."""
    tmp = None
    if work:
        os.makedirs(work, exist_ok=True)
        if os.listdir(work):
            raise SystemExit(f"ERROR: --work {work} is not empty")
        root = work
    else:
        tmp = tempfile.mkdtemp(prefix="owngenome-live-demo.")
        root = tmp
    try:
        held = prepare(root, genome_dir, records)
        frames, fam, app = replay(root, held, drop_hook)
        write_site(out, root, frames, fam)
        return {"frames": frames, "family": fam, "hooks": [n for n, _ in app.state_hooks]}
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True, help="the folder to write (replaced whole)")
    ap.add_argument("--genome-dir", help="a folder holding DEMO-001 (default: write one now)")
    ap.add_argument("--records", type=int, default=20000, help="variant records of a DEMO-001 written now")
    ap.add_argument("--work", help="an empty folder for the run folder, kept afterwards (default: a temporary one)")
    a = ap.parse_args(argv)
    r = build(a.out, a.genome_dir, a.records, a.work)
    size = sum(os.path.getsize(os.path.join(p, f)) for p, _, fs in os.walk(a.out) for f in fs)
    print(f"Demo written: {a.out} ({len(r['frames'])} frames, {size // 1024} KB; plugins {', '.join(r['hooks'])})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
