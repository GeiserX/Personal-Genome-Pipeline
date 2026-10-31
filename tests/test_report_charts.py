#!/usr/bin/env python3
"""The report's charts (bin/report_charts.py) and where bin/render_report.py
puts them, on the invented sections of tests/fixtures/report_charts/ and on
DEMO-001, the invented demo sample.

Checks, each with a red control that feeds it a broken input:
  1. ClinVar: one bar per star 4 to 0 with the by_stars count, and one grey
     bar with the low-confidence count; no green in the chart; with no hit a
     grey frame says that no hit is not an all-clear;
  2. PRS: no strip without a panel (adjusted false), only the report's note;
     with a panel a score with no percentile, or one of the other sex, says
     "no percentile" and gets no dot;
  3. fine-grained ancestry: the groups are render_report.ancestry_fine_rows,
     and each dot lies inside its range line;
  4. genome: a run of homozygosity whose region cannot be read stops the
     chart; a run under 5 MB is not drawn; ClinVar hits on chrM and on alt
     contigs are counted in not_drawn; ROH_MIN_MB is collect_summary's;
  5. mito: one dot per heteroplasmic site, the poly-C call marked apart;
  6. no caveat sentence is written in report_charts.py: they come from
     render_report's *_NOTE constants;
  7. the HTML report: no <script>, every href="#..." of the page resolves to
     an id, each chart sits in its section below that section's caveat, and a
     hit past the 50-row table links to the section;
  8. --charts-json: both report paths (modules/local/html_report/main.nf and
     scripts/24-html-report.sh) pass it with <sample>_charts.json beside the
     report, and two renders of the demo give the same file apart from
     generated_utc;
  9. run_delta.py sees no change between a summary without
     heteroplasmic_sites and one with it;
 10. the live page (bin/pgp_live.py with bin/pgp_live_charts.py and
     assets/live/plugins/charts.js): the plugin serves <run>_charts.json and
     names it in state.json; in a headless Chrome, signed in, the five charts
     appear, a click on a ClinVar tick opens the panel with its caveat and a
     viewer/#<locus> link, the page fetched nothing off its origin and the
     policy blocked nothing. Red control: the same page without charts.js
     fails the check. Without bin/pgp_live.py on the branch this part says
     SKIP; without Chrome it says SKIP, and FAIL in CI.
"""
import contextlib
import copy
import io
import json
import os
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402
import report_charts  # noqa: E402
import run_delta  # noqa: E402

FAILS = 0
GREENS = ("#d5f5e3", "#1d3b2a", "#8ee0ae", "#196f3d", "green")


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def fixture():
    with open(os.path.join(REPO, "tests", "fixtures", "report_charts", "sections.json")) as f:
        return {"sample": "S", "generated_utc": "2026-01-01 00:00 UTC", "sections": json.load(f)["sections"]}


def charts(s):
    return {c["id"]: c for c in render_report.chart_specs(s)}


# --- the check functions, each also run on a broken input -----------------------------------

def clinvar_problems(sp, data):
    out = []
    want = [int(data["by_stars"][k]) for k in ("4", "3", "2", "1", "0")] + [int(data["low_confidence"])]
    got = [m["value"] for m in sp["marks"]]
    if got != want:
        out.append(f"bars {got} != by_stars and low_confidence {want}")
    if sp["marks"][-1]["tone"] != "muted":
        out.append("the low-confidence bar is not grey")
    if any(m["tone"] not in ("attention", "muted") for m in sp["marks"]):
        out.append("a bar is neither yellow nor grey")
    low = sp["svg"].lower()
    out += [f"green {g} in the chart" for g in GREENS if g in low]
    return out


def prs_problems(sp, scores):
    """A dot per score with a percentile of the declared sex, 'no percentile' for the others."""
    out = []
    for m, r in zip(sp["marks"], scores):
        placeable = bool(r["percentile"]) and r["sex_check"] != render_report.PRS_NOT_APPLICABLE
        if placeable and m["value"] != float(r["percentile"]):
            out.append(f"{r['pgs_id']}: value {m['value']} is not the percentile {r['percentile']}")
        if not placeable and m["value"] is not None:
            out.append(f"{r['pgs_id']}: placed at {m['value']} with no percentile of the declared sex")
    dots = len(re.findall(r"<circle ", sp["svg"]))
    want = sum(1 for m in sp["marks"] if m["value"] is not None)
    if dots != want:
        out.append(f"{dots} dots for {want} percentiles")
    if sp["svg"].count(">no percentile</text>") != len(sp["marks"]) - want:
        out.append("not every score without a percentile says 'no percentile'")
    return out


def whisker_problems(svg):
    """Each mark with a range: its dot between the two ends of its line."""
    out = []
    for g in re.findall(r"<a [^>]*>.*?</a>", svg):
        ln = re.search(r'class="whisker" x1="([-0-9.]+)" y1="[-0-9.]+" x2="([-0-9.]+)"', g)
        dot = re.search(r'<circle cx="([-0-9.]+)"', g)
        if ln and dot and not float(ln.group(1)) <= float(dot.group(1)) <= float(ln.group(2)):
            out.append(f"dot {dot.group(1)} outside {ln.group(1)}..{ln.group(2)}")
    return out


def note_literals(src):
    """Caveat sentences written into a source file instead of passed in."""
    found = []
    for k, v in render_report.__dict__.items():
        if k.endswith("_NOTE") and isinstance(v, str) and len(v) > 30 and v[:30] in src:
            found.append(k)
    if re.search(r"^[A-Z_]+_NOTE\s*=", src, re.M):
        found.append("a *_NOTE constant")
    return found


def href_problems(page):
    ids = set(re.findall(r'\sid="([^"]+)"', page))
    return sorted({h for h in re.findall(r'href="#([^"]*)"', page) if h not in ids})


# --- the checks -----------------------------------------------------------------------------

def unit_checks():
    s = fixture()
    c = charts(s)
    cv = c["clinvar-stars"]
    check("1. ClinVar: bars equal by_stars, then the low-confidence count, yellow and grey only",
          clinvar_problems(cv, s["sections"]["clinvar"]["data"]) == [], clinvar_problems(cv, s["sections"]["clinvar"]["data"]))
    keep = report_charts.TONE["attention"]
    report_charts.TONE["attention"] = "#d5f5e3"
    try:
        bad = charts(s)["clinvar-stars"]
    finally:
        report_charts.TONE["attention"] = keep
    check("1. red control: a green bar is caught", clinvar_problems(bad, s["sections"]["clinvar"]["data"]) != [])
    tampered = copy.deepcopy(cv)
    tampered["marks"][1]["value"] += 1
    check("1. red control: a bar off by one is caught",
          clinvar_problems(tampered, s["sections"]["clinvar"]["data"]) != [])
    z = copy.deepcopy(s)
    z["sections"]["clinvar"]["data"].update(by_stars={k: 0 for k in "43210"}, hits=[], count=0, low_confidence=0)
    zs = charts(z)["clinvar-stars"]["svg"]
    check("1. no hit: a grey frame that says no hit is not an all-clear",
          "No hit is not an all-clear" in zs and "<a " not in zs and not any(g in zs.lower() for g in GREENS), zs)
    past = [m for m in c["genome"]["marks"] if m.get("kind") == "hit"]
    hits = s["sections"]["clinvar"]["data"]["hits"]
    row = {f"{h['chrom']}:{h['pos']}": i for i, h in enumerate(hits)}
    want = ["clinvar-hits" if row[m["locus"]] < 50 else "clinvar" for m in past]
    check("7. a hit past the 50-row table links to the section, the others to the table",
          [m["anchor"] for m in past] == want and "clinvar" in want, [m["anchor"] for m in past][-5:])

    scores = s["sections"]["prs"]["data"]["scores"]
    pr = c["prs"]
    check("2. PRS with a panel: a dot per percentile, 'no percentile' for '' and the other sex",
          pr["kind"] == "strip" and prs_problems(pr, scores) == [], prs_problems(pr, scores))
    raw = copy.deepcopy(pr)
    raw["marks"][1]["value"] = float(scores[1]["score"])
    raw["svg"] = report_charts.svg(raw)
    check("2. red control: the raw Score_SUM drawn as a dot is caught", prs_problems(raw, scores) != [])
    other = copy.deepcopy(pr)
    other["marks"][2]["value"] = float(scores[2]["percentile"])
    other["svg"] = report_charts.svg(other)
    check("2. red control: the other sex's percentile drawn is caught", prs_problems(other, scores) != [])
    na = copy.deepcopy(s)
    na["sections"]["prs"]["data"]["adjusted"] = False
    nc = charts(na)["prs"]
    check("2. PRS without a panel: no strip, the report's own note",
          nc["kind"] == "text" and nc["marks"] == [] and "<circle" not in nc["svg"]
          and nc["text"] == render_report.prs_note(na["sections"]["prs"]["data"]), nc["kind"])

    an = c["ancestry-fine"]
    d = s["sections"]["ancestry_fine"]["data"]
    want = [r["group"] for lv, _ in render_report.ANCESTRY_FINE_LEVELS for r in render_report.ancestry_fine_rows(d, lv)]
    check("3. ancestry: the groups are ancestry_fine_rows (a 0.4% group left out, a small regional one with a wide range kept)",
          [m["label"] for m in an["marks"]] == want and "Invented region B" in want and "Invented group C" not in want,
          [m["label"] for m in an["marks"]])
    check("3. ancestry: every dot inside its range line", whisker_problems(an["svg"]) == [] and "whisker" in an["svg"],
          whisker_problems(an["svg"]))
    off = copy.deepcopy(an)
    off["marks"][0]["value"] = 99.0
    off["svg"] = report_charts.svg(off)
    check("3. red control: a dot moved outside its range is caught", whisker_problems(off["svg"]) != [])

    g = c["genome"]
    roh = [m for m in g["marks"] if m.get("kind") == "roh"]
    nd = {x["what"]: x["count"] for x in g["not_drawn"]}
    check("4. genome: the two runs of 5 MB or more drawn, the 4.9 MB one not",
          [m["locus"] for m in roh] == ["chr3:10000000-22300000", "chr15:40000000-46100000"]
          and nd.get("runs of homozygosity under 5 MB") == 1, (roh, nd))
    check("4. genome: the chrM hit and the alt-contig hit are counted in not_drawn, not drawn",
          nd.get("ClinVar hits on chrM") == 1 and nd.get("ClinVar hits on other contigs (alt, unplaced)") == 1
          and "Not drawn:" in g["svg"]
          and not any(m["chrom"] in ("chrM", "chr22_KI270879v1_alt") for m in g["marks"]), nd)
    check("4. genome: the 5 MB floor is collect_summary.ROH_MIN_BP",
          report_charts.ROH_MIN_MB * 1e6 == collect_summary.ROH_MIN_BP)
    for region in ("chr3:10000000_22300000", "3:1-2", "chr3:20-10"):
        bad = copy.deepcopy(s)
        bad["sections"]["roh"]["data"]["autosomal_over_5mb"][0]["region"] = region
        try:
            charts(bad)
            stopped = False
        except ValueError:
            stopped = True
        check(f"4. an unreadable run of homozygosity ({region}) stops the chart", stopped)
    try:
        charts(s)
        ran = True
    except ValueError:
        ran = False
    check("4. control: the readable regions do not stop it", ran)

    mi = c["mito"]
    check("5. mito: one dot per heteroplasmic site, the poly-C call drawn apart",
          mi["svg"].count('class="site"') == 2 and mi["svg"].count('class="polyc"') == 1, mi["marks"])
    old = copy.deepcopy(s)
    del old["sections"]["mito"]["data"]["heteroplasmic_sites"]
    om = charts(old)["mito"]
    check("5. mito: a summary without the site list draws none and says how many it leaves out",
          'class="site"' not in om["svg"] and om["not_drawn"][0]["count"] == 2, om["not_drawn"])

    with open(os.path.join(REPO, "bin", "report_charts.py")) as f:
        src = f.read()
    check("6. no caveat sentence or *_NOTE constant in report_charts.py", note_literals(src) == [], note_literals(src))
    check("6. red control: a copied caveat is caught",
          note_literals(src + "\nX = '" + render_report.ROH_NOTE + "'\n") == ["ROH_NOTE"])
    check("6. each chart's caveat is render_report's text",
          all(sp["caveat"] == getattr(render_report, sp["caveat_key"]) for sp in c.values()))


def demo_checks(work):
    sys.path.insert(0, os.path.join(REPO, "tests", "demo"))
    import make_demo_sample
    g = f"{work}/demo"
    os.makedirs(g)
    make_demo_sample.traits = lambda d: None
    with contextlib.redirect_stdout(io.StringIO()):
        make_demo_sample.main(["--genome-dir", g, "--records", "20000"])
    sd = f"{g}/DEMO-001"
    outs = []
    for i in (1, 2):
        with contextlib.redirect_stdout(io.StringIO()):
            rc = render_report.main(["--sample", "DEMO-001", "--sample-dir", sd, "--json", f"{work}/s{i}.json",
                                     "-o", f"{work}/r{i}.html", "--charts-json", f"{work}/c{i}.json"])
        with open(f"{work}/c{i}.json") as f:
            outs.append(json.load(f))
        check(f"8. render {i} of the demo exits 0", rc == 0)
    a, b = (dict(x, generated_utc=None) for x in outs)
    check("8. charts.json of two renders is the same apart from generated_utc", a == b and outs[0]["charts"])
    b2 = copy.deepcopy(b)
    b2["charts"][0]["marks"][0]["value"] = -1
    check("8. red control: one changed value is caught", a != b2)
    check("8. charts.json carries its schema version and the five charts of the demo",
          outs[0]["schema_version"] == 1 and [c["id"] for c in outs[0]["charts"]]
          == ["clinvar-stars", "genome", "mito", "ancestry-fine", "prs"], [c["id"] for c in outs[0]["charts"]])

    with open(os.path.join(REPO, "modules", "local", "html_report", "main.nf")) as f:
        nf = f.read()
    with open(os.path.join(REPO, "scripts", "24-html-report.sh")) as f:
        sh = f.read()
    check("8. the Nextflow module writes and emits ${id}_charts.json beside the report",
          "--charts-json ${id}_charts.json" in nf and 'path("${meta.id}_charts.json")' in nf
          and "${meta.id}_charts.json" in nf.split("stub:")[1])
    check("8. step 24 writes ${SAMPLE}_charts.json beside the report",
          '--charts-json "/genome/${SAMPLE}/${SAMPLE}_charts.json"' in sh)

    # 7. the offline page, with the fixture's sections on the demo (53 hits, ROH, poly-C)
    with open(f"{work}/s1.json") as f:
        s = json.load(f)
    s["sections"].update(fixture()["sections"])
    page = render_report.html_report(s)
    check("7. the HTML report has no <script", "<script" not in page.lower())
    check("7. every href=\"#...\" resolves to an id of the page", href_problems(page) == [], href_problems(page))
    check("7. red control: a link to a missing id is caught",
          href_problems(page.replace('id="clinvar-hits"', 'id="gone"')) == ["clinvar-hits"])
    for sid, key in (("clinvar", "CLINVAR_FREQ_NOTE"), ("roh", "ROH_NOTE"), ("mito", "MITO_NOTE"),
                     ("ancestry-fine", "ANCESTRY_FINE_NOTE"), ("prs", "PRS_ANCESTRY_NOTE")):
        sec = re.search(rf'<section class="sec[^"]*" id="{sid}">(.*?)\n  </section>', page, re.S).group(1)
        cav = sec.find(render_report.E(getattr(render_report, key)))
        fig = sec.find('<figure class="chart">')
        check(f"7. section {sid}: one chart, below its caveat",
              sec.count('<figure class="chart">') == 1 and -1 < cav < fig, (cav, fig))
    check("7. the hits table carries the id the charts link to, and its rows keep their shape",
          page.count('<a id="clinvar-hits"></a>') == 1 and len(re.findall(r"<tr><td>chr", page)) >= 50)

    # 9. run_delta with the additive field
    prev = copy.deepcopy(s)
    del prev["sections"]["mito"]["data"]["heteroplasmic_sites"]
    check("9. run_delta: no change for the added heteroplasmic_sites", run_delta.compare(prev, {}, s, {}) == [],
          run_delta.compare(prev, {}, s, {}))
    moved = copy.deepcopy(s)
    moved["sections"]["clinvar"]["data"]["hits"][0]["significance"] = "Benign"
    check("9. red control: a reclassified hit is a change", run_delta.compare(prev, {}, moved, {}) != [])


def page_problems(b, base, run):
    """What is wrong with the live page's charts after a sign-in, from a click on a ClinVar tick."""
    import test_pgp_live as tl
    b.go(f"{base}/runs/{run}/")
    if not b.wait("location.pathname === '/login' && !!document.getElementById('password')"):
        return ["no login form"]
    b.eval(f"document.getElementById('password').value = {json.dumps(tl.PASSWORD)}; "
           "document.querySelector('form').submit(); true")
    if not b.wait(f"location.pathname === '/runs/{run}/' && document.documentElement.dataset.pgpStatus === 'ready'"):
        return ["the run page did not load"]
    out = []
    n = b.wait("document.documentElement.dataset.pgpCharts", 15)
    if n != "5":
        out.append(f"{n} charts shown, not 5")
    else:
        b.eval("document.querySelector('figure[data-chart=\"genome\"] a[data-mark^=\"hit-\"]')"
               ".dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true})); true")
        side = b.wait("(function () { var d = document.querySelector('.chart-detail'); if (!d || d.hidden) { "
                      "return null; } var v = d.querySelector('a.viewer'); return JSON.stringify({text: d.textContent, "
                      "viewer: v ? v.getAttribute('href') : '', at: location.pathname}); })()", 10)
        side = json.loads(side) if side else {}
        if render_report.CLINVAR_FREQ_NOTE not in side.get("text", ""):
            out.append(f"the panel has no caveat: {side}")
        if not re.fullmatch(r"viewer/#chr[0-9XY]+:[0-9]+", side.get("viewer", "")):
            out.append(f"no viewer/#<locus> link: {side}")
        if side.get("at") != f"/runs/{run}/":
            out.append(f"the click left the page: {side}")
    off = b.eval("performance.getEntriesByType('resource').map(function (e) { return e.name; })"
                 f".filter(function (u) {{ return u.indexOf({json.dumps(base + '/')}) !== 0; }})")
    if off:
        out.append(f"fetched off the origin: {off}")
    blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
    if blocked:
        out.append(f"the policy blocked: {blocked}")
    return out


def live_checks(work, charts_file):
    import pgp_live_charts

    class FakeApp:
        globs, hooks = [], []

        def run_file(self, g):
            self.globs.append(g)

        def state_hook(self, name, fn):
            self.hooks.append((name, fn))
    fa = FakeApp()
    pgp_live_charts.register(fa)
    d = os.path.join(work, "plugin", "RUNX")
    os.makedirs(d)
    before = pgp_live_charts.charts_state(d, {}, [])
    shutil.copy(charts_file, os.path.join(d, "RUNX_charts.json"))
    after = pgp_live_charts.charts_state(d, {}, [])
    check("10. the plugin serves *_charts.json and names the run's file in state.json, once it exists",
          fa.globs == ["*_charts.json"] and [n for n, _ in fa.hooks] == ["charts"] and before == {}
          and after.get("file") == "RUNX_charts.json" and isinstance(after.get("mtime"), int), (fa.globs, before, after))

    if not os.path.isfile(os.path.join(REPO, "bin", "pgp_live.py")):
        print("[SKIP] 10. the live page in Chrome: bin/pgp_live.py is not on this branch")
        return
    sys.path.insert(0, os.path.join(REPO, "tests"))
    import pgp_live
    import test_pgp_live as tl
    chrome = tl.find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "10. no Chrome or Chromium to load the live page in")
        if ci:
            check("10. Chrome in CI", False)
        return
    root = tl.make_root(work)
    shutil.copy(charts_file, os.path.join(root, tl.RUN, f"{tl.RUN}_charts.json"))
    pw = os.path.join(work, "password")
    with open(pw, "w") as f:
        f.write(pgp_live.hash_password(tl.PASSWORD) + "\n")
    bare = os.path.join(work, "assets-without-charts")
    shutil.copytree(os.path.join(REPO, "assets", "live"), bare)
    os.remove(os.path.join(bare, "plugins", "charts.js"))
    for tag, assets in (("live", os.path.join(REPO, "assets", "live")), ("control", bare)):
        argv = ["--root", root, "--state-dir", os.path.join(work, f"state-{tag}"), "--port", "0", "--ingest-port",
                "0", "--password-file", pw, "--assets", assets, "--quiet"]
        with contextlib.redirect_stderr(io.StringIO()):
            srv = tl.serve(argv, pgp_live.App())
        b = tl.Chrome(chrome, os.path.join(work, f"chrome-{tag}"))
        try:
            problems = page_problems(b, f"http://127.0.0.1:{srv[1].server_address[1]}", tl.RUN)
        finally:
            b.close()
            tl.stop(srv)
        if tag == "live":
            check("10. live page: five charts; a click on a ClinVar tick opens the panel with its caveat and a "
                  "viewer/# link; nothing fetched off the origin; nothing blocked", problems == [], problems)
        else:
            check("10. red control: the same page without charts.js is caught", problems != [], problems)


def main():
    unit_checks()
    with tempfile.TemporaryDirectory() as work:
        demo_checks(work)
        live_checks(work, f"{work}/c1.json")
    print(f"\n{FAILS} failed" if FAILS else "\nAll checks passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
