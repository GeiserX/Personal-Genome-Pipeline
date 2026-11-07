#!/usr/bin/env python3
"""The OwnGenome release in both reports (bin/render_report.py).

A run_manifest.tsv with a `run pipeline_version` row puts "OwnGenome <version>"
next to the pipeline commit in the text report and in the HTML header. Red
control: the same manifest without the row prints no release, so the check
sees the row and not some other text.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RENDER = os.path.join(ROOT, "bin", "render_report.py")


def render(rows):
    with tempfile.TemporaryDirectory() as d:
        sample = os.path.join(d, "S1")
        os.makedirs(sample)
        with open(os.path.join(sample, "run_manifest.tsv"), "w") as f:
            for r in rows:
                f.write("\t".join(r) + "\n")
        txt, html = os.path.join(d, "r.txt"), os.path.join(d, "r.html")
        subprocess.run([sys.executable, RENDER, "--sample", "S1", "--sample-dir", sample,
                        "--json", os.path.join(d, "summary.json"), "-o", txt, "-o", html],
                       check=True, stdout=subprocess.DEVNULL)
        with open(txt) as a, open(html) as b:
            return a.read(), b.read()


base = [("run", "written_utc", "2026-11-06 12:00:00 UTC"), ("run", "written_by", "run-all.sh"),
        ("run", "git_commit", "0123456789abcdef0123456789abcdef01234567")]
fails = 0


def check(ok, what):
    global fails
    print(("PASS " if ok else "FAIL ") + what)
    fails += not ok


txt, html = render(base + [("run", "pipeline_version", "1.0.0")])
check("(OwnGenome 1.0.0, pipeline commit 0123456789abcdef0123456789abcdef01234567," in txt,
      "the text report names the release next to the commit")
check("OwnGenome release: 1.0.0 &nbsp;|&nbsp;" in html and "Pipeline commit: 0123456789abcdef" in html,
      "the HTML header names the release next to the commit")

txt, html = render(base)
check("OwnGenome 1." not in txt and "(pipeline commit 0123456789abcdef" in txt,
      "red control: without the row the text report names no release")
check("OwnGenome release:" not in html, "red control: without the row the HTML header names no release")

sys.exit(1 if fails else 0)
