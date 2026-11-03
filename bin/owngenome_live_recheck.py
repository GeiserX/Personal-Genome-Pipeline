"""owngenome_live_recheck.py: the database re-check on the live page, a plugin of bin/owngenome_live.py.

bin/owngenome_recheck.py writes <run>/recheck/state.json: which databases the run
used have a newer release, which steps that would change, and the last rerun.
This plugin lets the signed-in browser fetch that file
(/runs/<run>/recheck/state.json) and puts a short form of it in state.json
under plugins.recheck: a banner line and, per step, the badges the step map
shows ("newer ClinVar available", "redone 2026-10-05").
assets/live/plugins/recheck.js draws them. Release names and step numbers
only: the check holds no result.
"""
import json
import os
import re

ACT = {"rerun-ready": "a rerun picks up the newer {name}", "refresh-available": "newer {name} available"}


def short(label):
    """'ClinVar (step 06's file)' -> 'ClinVar'."""
    return re.sub(r"\s*\(.*\)$", "", label or "").strip() or "database"


def utc14(text):
    """'2026-09-02 12:00:00 UTC' or '2026-09-02T12:00:00Z' as '20260902120000'; '' without one."""
    return "".join(c for c in str(text or "") if c.isdigit())[:14]


def manifest_written(run_dir):
    """The written_utc row of <run_dir>/run_manifest.tsv, or '' without one."""
    try:
        with open(os.path.join(run_dir, "run_manifest.tsv"), errors="replace") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 3 and c[0] == "run" and c[1] == "written_utc":
                    return c[2]
    except OSError:
        pass
    return ""


def recheck_state(run_dir, state, events):
    """{file, mtime (ns), checked_utc, banner, steps: {step: [badge]}}, or {} while there is no check."""
    if not run_dir:
        return {}
    p = os.path.join(run_dir, "recheck", "state.json")
    try:
        mtime = os.stat(p).st_mtime_ns
        with open(p) as f:
            st = json.load(f)
        items = [i for i in st.get("items") or [] if isinstance(i, dict)]
    except (OSError, ValueError, AttributeError):
        return {}
    # A check older than the run's run_manifest.tsv describes an earlier run: what it calls newer may
    # be what this run used. Not shown, as in the reports (bin/run_delta.py read_recheck).
    if utc14(manifest_written(run_dir)) > utc14(st.get("checked_utc")):
        return {}
    steps = {}
    newer = []
    for i in items:
        c = i.get("class")
        if c in ACT:
            text = ACT[c].format(name=short(i.get("label") or i.get("item")))
            newer.append(text)
            for n in i.get("steps") or []:
                steps.setdefault(n, []).append(text)
        elif c in ("in-pipeline-release", "upstream-only", "error"):
            newer.append(f"{short(i.get('label') or i.get('item'))}: {c}")
    rr = st.get("rerun")
    rr = rr if isinstance(rr, dict) else {}
    done = str(rr.get("finished_utc") or "")[:10]
    pvr = rr.get("predicted_vs_redone")
    for n in ((pvr.get("redone") if isinstance(pvr, dict) else None) or []) if done else []:
        steps.setdefault(n, []).append(f"redone {done}")
    banner = (f"Database re-check of {st.get('checked_utc') or 'unknown date'}: "
              + ("; ".join(newer) if newer else "nothing newer than this run"))
    return {"file": "recheck/state.json", "mtime": mtime, "checked_utc": st.get("checked_utc"),
            "banner": banner, "steps": steps}


def register(app):
    app.run_file("recheck/state.json")
    app.state_hook("recheck", recheck_state)
