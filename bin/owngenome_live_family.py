"""owngenome_live_family.py: the family page of the live site (a bin/owngenome_live.py plugin).

A family folder under GENOME_DIR (one holding members.tsv or pedigree.tsv, as
steps 38 and 39 and the Nextflow family workflow leave it) is listed as a run
of kind "family". Its page shows the pedigree, drawn from that file at request
time by bin/pedigree_svg.py, each child's chromosome map
(<child>_segments.svg), the relationship check, the shared segments and the
couple carrier screen, as the steps wrote them. A sample's page links to the
family it is in.

  GET /family/<family>/<file>   pedigree.svg (drawn), or one of the family's
                                result files listed in state.json; any other
                                name is 404

state.json gets plugins.family: {"kind": "family", ...} on a family's page,
{"kind": "member", "families": [...]} on the page of a sample in a family,
null otherwise.

The plugin computes nothing about genetics: it reads the files the steps
wrote, counts their rows and draws the pedigree the file declares.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pedigree_svg  # noqa: E402

NAME = r"[A-Za-z0-9][A-Za-z0-9._-]*"
NAME_RE = re.compile("^%s$" % NAME)
MARKERS = ("members.tsv", "pedigree.tsv")
# The family's result files the page may open, besides the drawn pedigree and the
# per-child and per-couple files (named from the members).
FIXED = ("family_report.md", "ibd_report.md", "relationships.tsv", "ibd_pairs.tsv", "family_plan.tsv",
         "members.tsv", "pedigree.tsv")
TYPES = {".svg": "image/svg+xml", ".md": "text/plain; charset=utf-8", ".tsv": "text/plain; charset=utf-8"}


def inside(path, top):
    real, top = os.path.realpath(path), os.path.realpath(top)
    return real == top or real.startswith(top.rstrip(os.sep) + os.sep)


def is_family(d):
    return bool(d) and os.path.isdir(d) and any(os.path.isfile(os.path.join(d, m)) for m in MARKERS)


def family_dirs(root):
    """[(name, dir)] of the family folders under root, sorted."""
    try:
        names = sorted(os.listdir(root))
    except OSError:
        return []
    out = []
    for n in names:
        d = os.path.join(root, n)
        if not n.startswith(".") and NAME_RE.match(n) and inside(d, root) and is_family(d):
            out.append((n, d))
    return out


def read_tsv(path):
    """Rows of a tab-separated file with a header line, as dicts; [] when it is not there."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path, errors="replace") as f:
        head = f.readline().rstrip("\r\n").lstrip("#").split("\t")
        for line in f:
            line = line.rstrip("\r\n")
            if line.strip() and not line.startswith("#"):
                vals = line.split("\t")
                rows.append(dict(zip(head, vals + [""] * (len(head) - len(vals)))))
    return rows


def members_of(d):
    order, members, _ = pedigree_svg.load(d)
    if not order:
        raise ValueError("the pedigree names no one")
    return order, members


def couple_names(order, members):
    """couple_<a>_<b> of each pair of partners, a and b sorted by name as step 39 names them."""
    out = []
    for a, b in pedigree_svg.couples(order, members):
        if members[a]["partner"] == b or members[b]["partner"] == a:
            x, y = sorted((a, b))
            out.append((x, y, "couple_%s_%s" % (x, y)))
    return out


def served_files(d, order, members):
    """The file names of the family folder the page may open: those that exist."""
    names = list(FIXED) + ["%s_segments.svg" % s for s in order]
    for _, _, base in couple_names(order, members):
        names += [base + ".md", base + ".tsv"]
    return sorted(n for n in set(names) if os.path.isfile(os.path.join(d, n)) and inside(os.path.join(d, n), d))


def family_state(name, d):
    try:
        order, members = members_of(d)
    except (OSError, ValueError) as e:
        return {"kind": "family", "error": "the pedigree could not be read: %s" % e}
    files = served_files(d, order, members)
    url = "/family/%s/" % name

    def link(f):
        return url + f if f in files else None

    people = [{"sample": s, "sex": members[s]["sex"], "father": members[s]["father"],
               "mother": members[s]["mother"], "partner": members[s]["partner"],
               "sequenced": members[s]["sequenced"]} for s in order]
    children = [{"child": p["sample"], "father": p["father"], "mother": p["mother"],
                 "segments": link("%s_segments.svg" % p["sample"])}
                for p in people if p["father"] or p["mother"]]
    rel = [{"a": r.get("sample_a"), "b": r.get("sample_b"), "declared": r.get("declared"),
            "verdict": r.get("verdict"), "reason": r.get("reason")}
           for r in read_tsv(os.path.join(d, "relationships.tsv"))]
    ibd = [{"a": r.get("sample_a"), "b": r.get("sample_b"), "declared": r.get("declared"),
            "shared_cm": r.get("shared_cm"), "degree": r.get("degree")}
           for r in read_tsv(os.path.join(d, "ibd_pairs.tsv"))]
    couples = []
    for a, b, base in couple_names(order, members):
        rows = read_tsv(os.path.join(d, base + ".tsv"))
        results = [r.get("result") for r in rows]
        couples.append({"a": a, "b": b, "report": link(base + ".md"), "table": link(base + ".tsv"),
                        "screened": bool(rows) or base + ".md" in files,
                        "shared": results.count("shared"), "x_linked": results.count("x_linked")})
    plan = [{"analysis": r.get("analysis"), "members": r.get("members"), "decision": r.get("decision"),
             "reason": r.get("reason")} for r in read_tsv(os.path.join(d, "family_plan.tsv"))]
    return {"kind": "family", "name": name, "pedigree": url + "pedigree.svg", "members": people,
            "children": children, "relationships": rel, "ibd_pairs": ibd, "couples": couples, "plan": plan,
            "reports": [{"name": f, "url": url + f} for f in files if f.endswith(".md")],
            "files": files}


def member_state(root, sample):
    fams = []
    for name, d in family_dirs(root):
        try:
            order, _ = members_of(d)
        except (OSError, ValueError):
            continue
        if sample in order:
            fams.append({"name": name, "url": "/runs/%s/" % name})
    return {"kind": "member", "families": fams} if fams else None


def register(app):
    def root():
        return os.path.realpath(app.config.root)

    def runs(app_):
        return [{"name": n, "dir": d, "kind": "family"} for n, d in family_dirs(root())]

    def hook(run_dir, state, events):
        if not run_dir:
            return None
        name = os.path.basename(os.path.normpath(run_dir))
        if is_family(run_dir) and inside(run_dir, root()):
            return family_state(name, run_dir)
        return member_state(root(), name)

    def serve(req, m):
        name, f = m.group(1), m.group(2)
        d = req.run_dir(name)
        if not d or not is_family(d) or not inside(d, root()):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        try:
            order, members = members_of(d)
        except (OSError, ValueError):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        if f == "pedigree.svg":
            app.log("family", name, f, 200)
            return req.send_body(200, TYPES[".svg"], pedigree_svg.draw(order, members, name))
        p = os.path.join(d, f)
        if f not in served_files(d, order, members) or not inside(p, d):
            app.log("family", name, f, 404)
            return req.send_body(404, "text/plain; charset=utf-8", "not found\n")
        with open(p, "rb") as fh:
            body = fh.read()
        app.log("family", name, f, 200)
        return req.send_body(200, TYPES.get(os.path.splitext(f)[1], "text/plain; charset=utf-8"), body)

    app.run_kind(runs)
    app.state_hook("family", hook)
    app.route("GET", r"/family/(%s)/(%s)" % (NAME, NAME), serve)
