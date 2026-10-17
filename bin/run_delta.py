#!/usr/bin/env python3
"""run_delta.py: what changed since the previous run of the same sample.

scripts/generate-report.sh runs it after the text report, before it keeps a
copy of this run's summary.json and run_manifest.tsv under
<sample>/history/<UTC time>/. It compares this run with the newest copy there
and writes <out>.txt, <out>.html and <out>.json. Standard library only, so it
runs in the plain python image.

  run_delta.py --sample S --summary summary.json --manifest run_manifest.tsv
               [--previous history/<stamp>] --out <sample_dir>/<sample>_delta

What it compares:
  inputs   the run manifest (bin/write_manifest.sh): versions.env values, the
           image each *_IMAGE names and the digest it resolved to, the data
           rows (ClinVar file date, VEP cache, CPSR bundle, pypgx bundle, HLA
           database, ancestry panel, PGS headers), the pipeline commit and the
           declared sex. When and by which script it was written are left out.
  sections a report section that ran in one run and not in the other.
  ClinVar  hits added or dropped (by chrom:pos ref>alt), a hit whose
           significance changed (reclassified), whose review status or stars
           changed (re-starred), or whose genotype changed.
  CPSR     the count of each classification.
  PGx      each CPIC gene's diplotype and phenotype, pypgx's CYP2D6 call and
           the CYP2D6 consensus.
  PRS      scores added or dropped, and a percentile that moved by
           PRS_BAND points or more, or appeared or went away.
For each result change it lists the changed inputs that feed that section, or
says that the inputs it tracks are the same and the call itself differs.
With no --previous it writes that there is no previous run on record.
Exits 2 when an input cannot be read.
"""
import argparse
import html
import json
import os
import sys

PRS_BAND = 5.0  # percentile points

# Changed inputs that can explain a change in a section: an input whose name
# contains one of these (case-insensitive). The pipeline commit explains any.
WHY_TOKENS = {
    "clinvar": ["clinvar", "DEEPVARIANT", "BCFTOOLS", "MINIMAP2", "BWA", "GATK", "FASTP"],
    "cpsr": ["PCGR", "VEP", "DEEPVARIANT"],
    "cpic": ["PHARMCAT", "PYPGX", "CYRIUS", "DEEPVARIANT", "BCFTOOLS"],
    "pypgx": ["PYPGX", "pypgx_bundle"],
    "cyp2d6": ["PYPGX", "pypgx_bundle", "CYRIUS", "PHARMCAT"],
    "prs": ["pgs:", "PGSC", "PLINK2", "ancestry_panel", "DEEPVARIANT"],
}
IGNORED_RUN_KEYS = ("written_utc", "written_by")


def read_manifest(path):
    """{(section, key): value} of a run_manifest.tsv; {} when there is none."""
    rows = {}
    if not path or not os.path.isfile(path):
        return rows
    with open(path, errors="replace") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3 and not c[0].startswith("#"):
                rows[(c[0], c[1])] = "\t".join(c[2:])
    return rows


def manifest_items(rows):
    """The manifest as {item: value}, one entry per thing that can change.
    An image is keyed by its *_IMAGE name, its value is 'ref (digest)', so a
    new tag and a new digest of the same tag are each one change."""
    out = {}
    for (sec, key), val in rows.items():
        if sec == "run":
            if key not in IGNORED_RUN_KEYS:
                out[key] = val
        elif sec == "version":
            if key.endswith("_IMAGE"):
                out[key] = f"{val} ({rows.get(('image', val), 'not in the manifest')})"
            else:
                out[key] = val
        elif sec == "data":
            out[key] = val
    return out


def change(area, kind, item, before, after, why=None):
    c = {"area": area, "kind": kind, "item": item, "before": before, "after": after}
    if why is not None:
        c["why"] = why
    return c


def diff_inputs(prev, cur):
    a, b = manifest_items(prev), manifest_items(cur)
    a = b  # RED CONTROL: ignore the previous manifest
    out = []
    for k in sorted(set(a) | set(b)):
        if a.get(k) != b.get(k):
            out.append(change("inputs", "input", k, a.get(k, "absent"), b.get(k, "absent")))
    return out


def why_for(area, inputs):
    tokens = [t.lower() for t in WHY_TOKENS.get(area, [])] + ["git_commit"]
    hit = [c for c in inputs if any(t in c["item"].lower() for t in tokens)]
    if not hit:
        return "the inputs this delta tracks for it are the same; the call itself differs"
    return "changed inputs: " + "; ".join(f"{c['item']} {c['before']} -> {c['after']}" for c in hit)


def section(s, key):
    return ((s.get("sections") or {}).get(key)) or {}


def has_data(sec):
    return sec.get("state") in ("ok", "stale")


def diff_states(p, c):
    out = []
    keys = sorted(set(p.get("sections") or {}) | set(c.get("sections") or {}))
    for k in keys:
        a, b = section(p, k), section(c, k)
        sa, sb = a.get("state", "absent"), b.get("state", "absent")
        if sa != sb:
            out.append(change(b.get("title") or a.get("title") or k, "section", k, sa, sb))
    return out


def hit_key(h):
    return f"{h.get('chrom')}:{h.get('pos')} {h.get('ref')}>{h.get('alt')}"


def diff_clinvar(p, c):
    a, b = section(p, "clinvar"), section(c, "clinvar")
    if not (has_data(a) and has_data(b)):
        return []
    ha = {hit_key(h): h for h in a.get("data", {}).get("hits", [])}
    hb = {hit_key(h): h for h in b.get("data", {}).get("hits", [])}
    out = []
    for k in sorted(set(ha) | set(hb)):
        x, y = ha.get(k), hb.get(k)
        gene = (y or x).get("gene", ".")
        item = f"{k} ({gene})"
        if x is None:
            out.append(change("ClinVar", "added", item, "absent", f"{y.get('significance')}, {y.get('stars')} stars"))
        elif y is None:
            out.append(change("ClinVar", "dropped", item, f"{x.get('significance')}, {x.get('stars')} stars", "absent"))
        else:
            if x.get("significance") != y.get("significance"):
                out.append(change("ClinVar", "reclassified", item, x.get("significance"), y.get("significance")))
            if x.get("stars") != y.get("stars") or x.get("review_status") != y.get("review_status"):
                out.append(change("ClinVar", "re-starred", item,
                                  f"{x.get('stars')} stars ({x.get('review_status')})",
                                  f"{y.get('stars')} stars ({y.get('review_status')})"))
            if x.get("genotype") != y.get("genotype"):
                out.append(change("ClinVar", "genotype", item, x.get("genotype"), y.get("genotype")))
    return out


def diff_cpsr(p, c):
    a, b = section(p, "cpsr"), section(c, "cpsr")
    if not (has_data(a) and has_data(b)):
        return []
    ca, cb = a.get("data", {}).get("classification") or {}, b.get("data", {}).get("classification") or {}
    return [change("CPSR", "classification count", k, ca.get(k, 0), cb.get(k, 0))
            for k in sorted(set(ca) | set(cb)) if ca.get(k, 0) != cb.get(k, 0)]


def diff_cpic(p, c):
    a, b = section(p, "cpic"), section(c, "cpic")
    if not (has_data(a) and has_data(b)):
        return []
    ga = {g.get("gene"): g for g in a.get("data", {}).get("genes", [])}
    gb = {g.get("gene"): g for g in b.get("data", {}).get("genes", [])}
    out = []
    for k in sorted(set(ga) | set(gb), key=str):
        x, y = ga.get(k), gb.get(k)
        fx = f"{x.get('diplotype')} ({x.get('phenotype')})" if x else "absent"
        fy = f"{y.get('diplotype')} ({y.get('phenotype')})" if y else "absent"
        if fx != fy:
            out.append(change("PGx (CPIC)", "call", k, fx, fy))
    return out


def diff_cyp2d6(p, c):
    out = []
    a, b = section(p, "pypgx"), section(c, "pypgx")
    if has_data(a) and has_data(b):
        x, y = a.get("data", {}).get("cyp2d6"), b.get("data", {}).get("cyp2d6")
        if x != y:
            out.append(change("PGx (pypgx)", "call", "CYP2D6", x, y))
    x = ((p.get("cyp2d6") or {}).get("consensus") or {}).get("result")
    y = ((c.get("cyp2d6") or {}).get("consensus") or {}).get("result")
    if x != y:
        out.append(change("PGx (CYP2D6 consensus)", "call", "CYP2D6", x or "none", y or "none"))
    return out


def pct(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def diff_prs(p, c):
    a, b = section(p, "prs"), section(c, "prs")
    if not (has_data(a) and has_data(b)):
        return []
    sa = {r.get("pgs_id"): r for r in a.get("data", {}).get("scores", [])}
    sb = {r.get("pgs_id"): r for r in b.get("data", {}).get("scores", [])}
    out = []
    for k in sorted(set(sa) | set(sb), key=str):
        x, y = sa.get(k), sb.get(k)
        item = f"{k} ({(y or x).get('condition', '')})"
        if x is None or y is None:
            out.append(change("PRS", "added" if x is None else "dropped", item,
                              "absent" if x is None else "scored", "scored" if x is None else "absent"))
            continue
        px, py = pct(x.get("percentile")), pct(y.get("percentile"))
        if (px is None) != (py is None) or (px is not None and abs(px - py) >= PRS_BAND):
            out.append(change("PRS", "percentile", item,
                              x.get("percentile") or "none", y.get("percentile") or "none"))
    return out


def compare(prev_summary, prev_manifest, summary, manifest):
    inputs = diff_inputs(prev_manifest, manifest)
    results = []
    for area, fn in (("clinvar", diff_clinvar), ("cpsr", diff_cpsr), ("cpic", diff_cpic),
                     ("cyp2d6", diff_cyp2d6), ("prs", diff_prs)):
        for ch in fn(prev_summary, summary):
            ch["why"] = why_for(area, inputs)
            results.append(ch)
    return inputs + diff_states(prev_summary, summary) + results


def load_json(path):
    with open(path) as f:
        return json.load(f)


def render_text(d):
    w = [f"What changed since the previous run: {d['sample']}", ""]
    if d["previous"] is None:
        w.append("No previous run of this sample on record (its history/ folder is empty).")
        w.append("This run is kept as the baseline the next run is compared with.")
        return "\n".join(w) + "\n"
    w.append(f"Previous run: history/{d['previous']} (summary of {d['previous_generated'] or 'unknown date'})")
    w.append(f"This run:     summary of {d['current_generated'] or 'unknown date'}")
    w.append("")
    if not d["changes"]:
        w.append("No change: the same inputs and the same results as the previous run.")
        return "\n".join(w) + "\n"
    w.append(f"{len(d['changes'])} change(s):")
    area = None
    for ch in d["changes"]:
        if ch["area"] != area:
            area = ch["area"]
            w.extend(["", f"## {area}"])
        w.append(f"  {ch['kind']}: {ch['item']}: {ch['before']} -> {ch['after']}")
        if ch.get("why"):
            w.append(f"    why: {ch['why']}")
    return "\n".join(w) + "\n"


def render_html(d):
    e = html.escape
    body = [f"<h1>What changed since the previous run: {e(d['sample'])}</h1>"]
    if d["previous"] is None:
        body.append("<p>No previous run of this sample on record (its history/ folder is empty). "
                    "This run is kept as the baseline the next run is compared with.</p>")
    else:
        body.append(f"<p>Previous run: history/{e(d['previous'])} (summary of {e(d['previous_generated'] or 'unknown date')}). "
                    f"This run: summary of {e(d['current_generated'] or 'unknown date')}.</p>")
        if not d["changes"]:
            body.append("<p>No change: the same inputs and the same results as the previous run.</p>")
        else:
            body.append(f"<p>{len(d['changes'])} change(s).</p>")
            body.append("<table><tr><th>Area</th><th>Change</th><th>Item</th><th>Before</th><th>Now</th><th>Why</th></tr>")
            for ch in d["changes"]:
                body.append("<tr>" + "".join(f"<td>{e(str(ch.get(k, '')))}</td>"
                                             for k in ("area", "kind", "item", "before", "after", "why")) + "</tr>")
            body.append("</table>")
    style = ("body{font-family:sans-serif;margin:2em;background:#111;color:#ddd}"
             "table{border-collapse:collapse}td,th{border:1px solid #444;padding:4px 8px;text-align:left;vertical-align:top}")
    return (f"<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><title>{e(d['sample'])}: delta</title>"
            f"<style>{style}</style></head><body>\n" + "\n".join(body) + "\n</body></html>\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--summary", required=True, help="this run's summary.json")
    ap.add_argument("--manifest", help="this run's run_manifest.tsv")
    ap.add_argument("--previous", help="the previous run's history folder (summary.json, run_manifest.tsv)")
    ap.add_argument("--out", required=True, help="output prefix: writes OUT.txt, OUT.html, OUT.json")
    a = ap.parse_args(argv)
    try:
        cur = load_json(a.summary)
        prev = load_json(os.path.join(a.previous, "summary.json")) if a.previous else None
    except (OSError, ValueError) as e:
        print(f"ERROR: cannot read a summary: {e}", file=sys.stderr)
        return 2
    d = {"sample": a.sample, "previous": os.path.basename(os.path.normpath(a.previous)) if a.previous else None,
         "previous_generated": (prev or {}).get("generated_utc"), "current_generated": cur.get("generated_utc"),
         "prs_band": PRS_BAND, "changes": []}
    if prev is not None:
        d["changes"] = compare(prev, read_manifest(os.path.join(a.previous, "run_manifest.tsv")),
                               cur, read_manifest(a.manifest))
    for ext, text in (("txt", render_text(d)), ("html", render_html(d)), ("json", json.dumps(d, indent=1) + "\n")):
        with open(f"{a.out}.{ext}", "w") as f:
            f.write(text)
    sys.stdout.write(render_text(d))
    return 0


if __name__ == "__main__":
    sys.exit(main())
