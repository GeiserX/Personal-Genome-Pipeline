"""owngenome_live_search.py: the live page's search box, a plugin of bin/owngenome_live.py.

Step 24 (and the Nextflow HTML_REPORT process) writes <sample>_search_index.json
beside the report: every gene symbol and rsID the results mention, with the
sections that name them (bin/results_index.py). This plugin lets the
signed-in browser fetch that file from the run folder
(/runs/<run>/<sample>_search_index.json) and tells the page, in state.json
under plugins.search, whether it is there and when it was written.
assets/live/plugins/search.js fetches it once (and again only when the report
is written anew) and filters it in the browser: typing sends no request, so a
query is in no request and no log. The file names genes of the results, so,
like the report, it is only served behind the password.
"""
import os


def search_state(run_dir, state, events):
    """{file, mtime (ns)} of the run's <run>_search_index.json, or {} while there is none."""
    if not run_dir:
        return {}
    name = f"{os.path.basename(os.path.normpath(run_dir))}_search_index.json"
    p = os.path.join(run_dir, name)
    try:
        return {"file": name, "mtime": os.stat(p).st_mtime_ns} if os.path.isfile(p) else {}
    except OSError:
        return {}


def register(app):
    app.run_file("*_search_index.json")
    app.state_hook("search", search_state)
