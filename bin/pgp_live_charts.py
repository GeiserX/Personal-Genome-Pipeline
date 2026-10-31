"""pgp_live_charts.py: the live page's charts, a plugin of bin/pgp_live.py.

Step 24 (and the Nextflow HTML_REPORT process) writes <sample>_charts.json
beside the report: the report's five charts as data and as SVG
(bin/report_charts.py). This plugin lets the signed-in browser fetch that file
from the run folder (/runs/<run>/<sample>_charts.json) and tells the page, in
state.json under plugins.charts, whether it is there and when it was written,
so assets/live/plugins/charts.js fetches it once and again only when the
report is written anew. The file holds findings: like the report, it is only
served behind the password.
"""
import os


def charts_state(run_dir, state, events):
    """{file, mtime} of the run's <run>_charts.json, or {} while there is none."""
    if not run_dir:
        return {}
    name = f"{os.path.basename(os.path.normpath(run_dir))}_charts.json"
    p = os.path.join(run_dir, name)
    try:
        return {"file": name, "mtime": int(os.path.getmtime(p))} if os.path.isfile(p) else {}
    except OSError:
        return {}


def register(app):
    app.run_file("*_charts.json")
    app.state_hook("charts", charts_state)
