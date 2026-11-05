#!/usr/bin/env python3
"""The live site's public demo (tests/demo/build_live_demo.py), built and checked.

Checks:
  1. tests/fixtures/live/demo-001.weblog.jsonl is what tests/demo/make_demo_events.py
     writes from assets/steps.tsv now;
  2. every state hook the plugins bin/owngenome_live_*.py register has a value
     in the last frame of DEMO-001's page (the family page for the family
     plugin's own view), none of them an error;
  3. the frames follow the run: running at the first, finished at the last with
     every replayed step done; the charts and the search index appear only
     once the report step has ended, the re-check only after the run;
  4. the run-health card measures the invented machine, never the machine that
     builds the demo, and no frame holds a path of that machine;
  5. every link a frame gives is relative and reaches a file of the demo;
  6. every page the build writes carries the server's Content-Security-Policy
     as a meta tag, and no page or script of the demo loads anything from
     outside;
  7. in a headless Chrome (skipped without Chrome, except in CI), served as a
     static host serves it (no policy header): the runs list, the replay to its
     last frame, the charts, a search, the re-check banner, the run-health
     card, the family link and its pedigree, the viewer page, and no policy
     violation.
Checks 1, 2, 3 and 5 also run against a planted mistake and must fail there (2:
a build with one plugin's state hook dropped): a check that cannot fail is no check.

Run: python3 tests/test_live_demo.py
"""
import copy
import json
import os
import shutil
import sys
import tempfile
import threading
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "demo"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import build_live_demo  # noqa: E402
import make_demo_events  # noqa: E402
import owngenome_live  # noqa: E402
import run_state  # noqa: E402
from test_owngenome_live import Chrome, asset_problems, find_chrome  # noqa: E402

FAILS = []
S, FAM = build_live_demo.SAMPLE, build_live_demo.FAMILY


def check(desc, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + desc + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS.append(desc)


def control(desc, still_passes):
    """A planted mistake: the check must fail on it."""
    check(f"control: {desc}", not still_passes, "the check passed on the planted mistake")


def registered_hooks():
    app = owngenome_live.App()
    owngenome_live.configure(["--root", tempfile.gettempdir()], app)
    return [n for n, _ in app.state_hooks]


def missing_hooks(built, hooks):
    """The hooks with no value (or an error) on the demo's pages."""
    last = built["frames"][-1]["state"]["plugins"]
    fam = built["family"]["state"]["plugins"]
    bad = []
    for h in hooks:
        v = last.get(h)
        if h == "family":
            v = v if (v or {}).get("kind") == "member" and (fam.get(h) or {}).get("kind") == "family" else None
        if not v or (isinstance(v, dict) and "error" in v):
            bad.append(h)
    return bad


def replayed_steps():
    procs = set()
    with open(build_live_demo.WEBLOG) as f:
        for line in f:
            tr = json.loads(line).get("trace") or {}
            if tr.get("name"):
                procs.add(run_state.task_name(tr["name"])[0])
    return sorted(s["step"] for s in run_state.read_steps() if set(s["processes"]) & procs)


def follow_problems(frames, steps):
    """What in the frames does not follow the run."""
    out = []
    first, last = frames[0]["state"], frames[-1]["state"]
    if first["run"]["state"] != "running":
        out.append(f"first frame {first['run']['state']}")
    if last["run"]["state"] != "finished":
        out.append(f"last frame {last['run']['state']}")
    st = {s["step"]: s["state"] for s in last["steps"]}
    out += [f"step {n} {st.get(n)} at the end" for n in steps if st.get(n) != "done"]
    for fr in frames:
        p, s24 = fr["state"]["plugins"], {s["step"]: s["state"] for s in fr["state"]["steps"]}.get("24")
        for k in ("charts", "search"):
            if bool((p.get(k) or {}).get("file")) != (s24 == "done"):
                out.append(f"{k} at {fr['t']} with step 24 {s24}")
        if p.get("recheck") and fr is not frames[-1]:
            out.append(f"re-check at {fr['t']}, before the run ended")
    if not (last["plugins"].get("recheck") or {}).get("steps"):
        out.append("no re-check badge at the end")
    return out


def links(v):
    if isinstance(v, dict):
        for k, x in v.items():
            if k in ("url", "pedigree", "segments", "report", "table") and isinstance(x, str):
                yield x
            else:
                yield from links(x)
    elif isinstance(v, list):
        for x in v:
            yield from links(x)


def link_problems(out, built):
    """Links of the frames that are absolute or reach no file of the demo."""
    bad = []
    pages = [(S, fr["state"]) for fr in built["frames"]] + [(FAM, built["family"]["state"])]
    for run, st in pages:
        base = os.path.join(out, "runs", run)
        for u in links(st["plugins"]):
            if u.startswith("/"):
                bad.append(f"{run}: absolute {u}")
                continue
            p = os.path.normpath(os.path.join(base, u))
            if not (os.path.isfile(p) or os.path.isfile(os.path.join(p, "index.html"))):
                bad.append(f"{run}: {u} reaches nothing")
        for k, v in st["plugins"].items():
            if isinstance(v, dict) and isinstance(v.get("file"), str) and not os.path.isfile(os.path.join(base, v["file"])):
                bad.append(f"{run}: {k} file {v['file']} not in the demo")
    return bad


def located_term(out):
    """A term of the demo's search index with a mention at a locus (its hit links the viewer)."""
    with open(os.path.join(out, "runs", S, f"{S}_search_index.json")) as f:
        terms = json.load(f)["terms"]
    return next((t for t, items in sorted(terms.items())
                 if any(m.get("locus") for it in items for m in it.get("mentions") or [])), None)


class Static(SimpleHTTPRequestHandler):
    """A static host: files only, no Content-Security-Policy header."""

    def log_message(self, *a):
        pass


def browser_checks(out):
    chrome = find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the demo in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    handler = type("H", (Static,), {"__init__": lambda self, *a, **k: SimpleHTTPRequestHandler.__init__(
        self, *a, directory=out, **k)})
    srv = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_address[1]}"
    b = Chrome(chrome, os.path.join(os.path.dirname(out), "chrome"))
    try:
        b.go(f"{base}/")
        runs = b.wait("document.querySelectorAll('ul.runs a').length === 2 && "
                      "Array.from(document.querySelectorAll('ul.runs a')).map(function (a) { return a.getAttribute('href'); }).join(' ')")
        check("browser: the runs list links the sample and the family", runs == f"runs/{S}/ runs/{FAM}/", runs)
        b.go(f"{base}/runs/{S}/")
        first = b.wait("document.documentElement.dataset.pgpFrame && document.documentElement.dataset.pgpStatus === 'ready' "
                       "&& document.querySelectorAll('#map g.box').length")
        b.eval("var s = document.getElementById('replay-slider'); s.value = s.max; "
               "s.dispatchEvent(new Event('input')); true")
        n = b.wait("document.documentElement.dataset.pgpFrame === document.getElementById('replay-slider').max "
                   "&& document.documentElement.dataset.pgpCharts")
        boxes = b.eval("document.querySelectorAll('#map g.box').length")
        check("browser: the replay draws the step map, which grows as steps are planned",
              isinstance(first, int) and first >= 1 and isinstance(boxes, int) and boxes > 30, (first, boxes))
        check("browser: at the last frame the charts are drawn from the demo's charts file", n and int(n) > 0, n)
        state = b.eval("document.getElementById('run-state').textContent")
        check("browser: the last frame shows the run finished", state == "finished", state)
        term = located_term(out)
        ready = b.wait("var i = document.querySelector('input[type=search]'); i && !i.disabled && 'enabled'")
        check("browser: the search box opens on the demo's index", ready == "enabled",
              b.eval("Array.from(document.querySelectorAll('section.plugin')).map(function (s) { "
                     "return s.dataset.plugin + ': ' + s.textContent.slice(0, 120); }).join(' | ')"))
        b.eval("var i = document.querySelector('input[type=search]'); i.value = %s; "
               "i.dispatchEvent(new Event('input')); true" % json.dumps(term))
        hits = b.wait("!document.querySelector('input[type=search]').disabled && document.documentElement.dataset.pgpSearch")
        check(f"browser: a search for {term} finds it in the demo's index", hits and int(hits) > 0, (term, hits))
        cards = b.eval("Array.from(document.querySelectorAll('section.plugin')).map(function (s) { "
                       "return s.dataset.plugin + '=' + (s.hidden ? '' : s.textContent.length); }).join(' ')")
        banner = b.eval("(document.querySelector('section.plugin[data-plugin=recheck]') || {}).textContent || ''")
        check("browser: the re-check card shows the newer ClinVar", "newer ClinVar" in banner, banner[:300])
        ops = b.eval("(document.querySelector('section.plugin[data-plugin=ops]') || {}).textContent || ''")
        check("browser: the run-health card says the run ended and names the retried task",
              "retried" in ops, ops[:400])
        famlink = b.eval("var a = document.querySelector('section.plugin[data-plugin=family] a'); a && a.getAttribute('href')")
        check("browser: the family card links the family page", famlink == f"../../runs/{FAM}/", (famlink, cards))
        viewer = b.eval("var a = document.querySelector('a.viewer'); a && a.href")
        if viewer:
            b.go(viewer)
            text = b.wait("document.body && document.body.textContent.indexOf('not in the demo') >= 0 && 'ok'")
            check("browser: a viewer link opens the page that says the viewer is not in the demo", text == "ok", viewer)
        else:
            check("browser: the charts or the search give a viewer link", False, cards)
        b.go(f"{base}/runs/{FAM}/")
        ped = b.wait("var i = document.querySelector('img.pedigree'); i && i.complete && i.naturalWidth > 0 && 'drawn'")
        check("browser: the family page shows the pedigree", ped == "drawn",
              b.eval("var i = document.querySelector('img.pedigree'); i && i.src"))
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the family page", not blocked, blocked)
        b.go(f"{base}/runs/{S}/")
        b.wait("document.documentElement.dataset.pgpStatus === 'ready'")
        b.eval("var s = document.getElementById('replay-slider'); s.value = s.max; s.dispatchEvent(new Event('input')); true")
        b.wait("document.documentElement.dataset.pgpCharts")
        time.sleep(1)
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the sample page", not blocked, blocked)
        csp = b.eval("var m = document.querySelector('meta[http-equiv=Content-Security-Policy]'); m && m.content")
        check("browser: the page runs under the meta policy", csp == build_live_demo.CSP, csp)
    finally:
        b.close()
        srv.shutdown()
        srv.server_close()


def main():
    work = tempfile.mkdtemp(prefix="test_live_demo.")
    try:
        # 1. the weblog fixture
        check("1. the weblog fixture is what make_demo_events.py writes now", make_demo_events.main(["--check"]) == 0)
        changed = os.path.join(work, "changed.jsonl")
        with open(build_live_demo.WEBLOG) as f:
            lines = f.readlines()
        with open(changed, "w") as f:
            f.writelines(lines[:-1])
        control("a fixture missing its last message", make_demo_events.main(["--check", "--out", changed]) == 0)

        out = os.path.join(work, "demo")
        built = build_live_demo.build(out, work=os.path.join(work, "root"))
        frames = built["frames"]
        print(f"built: {len(frames)} frames, hooks {built['hooks']}")

        # 2. every plugin's state hook
        hooks = registered_hooks()
        check("2. the build calls every registered state hook", sorted(built["hooks"]) == sorted(hooks), (built["hooks"], hooks))
        bad = missing_hooks(built, hooks)
        check(f"2. every state hook ({', '.join(hooks)}) has a value on the demo's pages", not bad and len(hooks) >= 5, bad)
        dropped = build_live_demo.build(os.path.join(work, "demo-dropped"), work=os.path.join(work, "root-dropped"),
                                        drop_hook=hooks[0])
        control(f"a build that drops the {hooks[0]} hook", not missing_hooks(dropped, hooks))

        # 3. the frames follow the run
        steps = replayed_steps()
        bad = follow_problems(frames, steps)
        check(f"3. the frames follow the run ({len(frames)} frames, {len(steps)} steps)", not bad and len(steps) > 30, bad)
        early = copy.deepcopy(frames)
        early[0]["state"]["plugins"]["charts"] = early[-1]["state"]["plugins"]["charts"]
        control("charts shown before the report step ended", not follow_problems(early, steps))

        # 4. the invented machine, no path of this one
        ops = frames[-1]["state"]["plugins"]["ops"]
        disk = ops.get("disk") or [{}]
        check("4. the run-health card measures the invented machine",
              ops["memory"]["host_bytes"] == build_live_demo.DEMO_RAM
              and disk[0].get("size_bytes") == build_live_demo.DEMO_DISK.f_blocks * build_live_demo.DEMO_DISK.f_frsize,
              (ops["memory"], disk))
        with open(os.path.join(out, "runs", S, "frames.json")) as f:
            text = f.read()
        leaks = [p for p in (work, os.path.realpath(work), REPO, os.path.expanduser("~")) if p in text]
        check("4. no frame holds a path of the machine that built it", not leaks, leaks)

        # 5. links
        bad = link_problems(out, built)
        check("5. every link of the frames is relative and reaches a file of the demo", not bad, bad)
        planted = copy.deepcopy(built)
        planted["family"]["state"]["plugins"]["family"]["pedigree"] = f"/family/{FAM}/pedigree.svg"
        control("an absolute pedigree link", not link_problems(out, planted))

        # 6. the policy on every page, nothing from outside
        pages, texts = [], {}
        for p, _, fs in os.walk(out):
            for f in fs:
                rel = os.path.relpath(os.path.join(p, f), out)
                if f.endswith((".html", ".js", ".css")) and rel != os.path.join("runs", S, "report.html"):
                    with open(os.path.join(p, f)) as fh:
                        texts[rel] = fh.read()
                    if f.endswith(".html"):
                        pages.append(rel)
        meta = f'<meta http-equiv="Content-Security-Policy" content="{build_live_demo.CSP}">'
        nometa = [p for p in pages if meta not in texts[p]]
        check(f"6. every page the build writes carries the policy ({len(pages)} pages)", not nometa and len(pages) >= 4, nometa)
        check("6. the policy is the server's, frame-ancestors aside",
              build_live_demo.CSP + "; frame-ancestors 'none'" == owngenome_live.CSP, build_live_demo.CSP)
        bad = asset_problems(texts)
        check("6. no page or script of the demo loads anything from outside", not bad, bad)

        # 7. in a browser
        browser_checks(out)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
