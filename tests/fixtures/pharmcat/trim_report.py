#!/usr/bin/env python3
"""trim_report.py: cut a PharmCAT report.json down to a test fixture.

  python3 tests/fixtures/pharmcat/trim_report.py REPORT.json OUT.json [MAX_ANNOTATIONS]

Keeps every gene with all its listed diplotypes (sourceDiplotypes, or
recommendationDiplotypes when that is all there is) and relatedDrugs, and for
each drug of each guidance source up to MAX_ANNOTATIONS (default 4) of the
annotations PharmCAT matched, with only the fields bin/pgx_parse.py reads.
Writes compact JSON. See README.md for how each fixture was made.
"""
import json
import sys

src, dst = sys.argv[1], sys.argv[2]
MAX_ANN = int(sys.argv[3]) if len(sys.argv) > 3 else 4
d = json.load(open(src))


def dip(x):
    if not isinstance(x, dict): return x
    out = {k: x[k] for k in ("gene", "label", "phenotypes", "activityScore", "outsidePhenotype") if k in x}
    for a in ("allele1", "allele2"):
        if a in x:
            out[a] = {k: x[a][k] for k in ("gene", "name", "function") if x[a] and k in x[a]} if x[a] else x[a]
    return out


genes = {}
for g, v in d["genes"].items():
    genes[g] = {k: v[k] for k in ("geneSymbol", "chr", "phased", "callSource") if k in v}
    genes[g]["relatedDrugs"] = [{"name": r["name"], "id": r.get("id")} for r in v.get("relatedDrugs") or []]
    # sourceDiplotypes is what pgx_parse.py reads; recommendationDiplotypes only
    # when sourceDiplotypes is absent, so it is kept only then.
    k = "sourceDiplotypes" if v.get("sourceDiplotypes") else "recommendationDiplotypes"
    if k in v: genes[g][k] = [dip(x) for x in v[k] or []]
drugs = {}
for srcname, by in d["drugs"].items():
    for name, rep in by.items():
        gls = []
        for gl in rep.get("guidelines") or []:
            anns = []
            for a in (gl.get("annotations") or [])[:MAX_ANN]:
                anns.append({"phenotypes": a.get("phenotypes"), "classification": a.get("classification"),
                             "drugRecommendation": a.get("drugRecommendation"),
                             "genotypes": [{"diplotypes": [{"gene": x.get("gene"), "label": x.get("label")}
                                                           for x in (gt.get("diplotypes") or []) if isinstance(x, dict)]}
                                           for gt in (a.get("genotypes") or [])]})
            if anns: gls.append({"name": gl.get("name"), "source": gl.get("source"), "annotations": anns})
        if gls:
            drugs.setdefault(srcname, {})[name] = {"name": rep.get("name"), "id": rep.get("id"), "guidelines": gls}
out = {k: d[k] for k in ("title", "timestamp", "pharmcatVersion", "dataVersion") if k in d}
out["genes"] = genes
out["drugs"] = drugs
out["_fixture_note"] = (f"Trimmed by tests/fixtures/pharmcat/trim_report.py: every gene, and up to {MAX_ANN} matched "
                        "annotations per drug guideline, with only the fields bin/pgx_parse.py reads.")
json.dump(out, open(dst, "w"), separators=(",", ":"))
