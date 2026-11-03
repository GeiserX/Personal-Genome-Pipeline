#!/usr/bin/env python3
"""bin/owngenome_live_family.py and assets/live/plugins/family.js on tests/fixtures/families_live.

Checks:
  1. the family folders (family_a with members.tsv, family_b with pedigree.tsv
     only) are listed as runs of kind "family"; the sample folder child_a is
     not;
  2. family_a's state: its members, each child with its parents and its
     chromosome map (child_a has one, child_b none), the relationship check,
     the shared segments, the couple screen counted from step 39's table, the
     plan, and exactly the result files that exist;
  3. child_a's state links to family_a; a sample in no family gets null;
  4. GET /family/<family>/pedigree.svg is the drawing of the declared parents;
     a listed file is served as written; notes.txt (in the folder, written by
     no step), a name with '..', a symlink out of the folder, a sample's
     folder and an unknown family are 404;
  5. family.js loads nothing from outside, sets no innerHTML and has no eval;
     node --check reads it when node is installed;
  6. with bin/owngenome_live.py present (the live-server change), the same through
     the real server: the run list, state.json's plugins.family, the route
     behind the login, and a family page that loads the script; in a headless
     Chrome (skipped without Chrome, except in CI) the card, the pedigree and a
     child's map load and the page's policy blocks nothing.
Checks 2 and 4 are run again on a copy of the fixture with two parents swapped
and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_owngenome_live_family.py
"""
import http.client
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import types
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import owngenome_live_family  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "families_live")
JS = os.path.join(REPO, "assets", "live", "plugins", "family.js")
NS = "{http://www.w3.org/2000/svg}"
FAILS = []

CHILDREN = [{"child": "child_a", "father": "father_a", "mother": "mother_a",
             "segments": "/family/family_a/child_a_segments.svg"},
            {"child": "child_b", "father": "father_a", "mother": "mother_a", "segments": None},
            {"child": "father_a", "father": "gf_a", "mother": "gm_a", "segments": None},
            {"child": "mother_a", "father": None, "mother": "gm_m", "segments": None}]
FILES = ["child_a_segments.svg", "couple_father_a_mother_a.md", "couple_father_a_mother_a.tsv", "family_plan.tsv",
         "family_report.md", "ibd_pairs.tsv", "ibd_report.md", "members.tsv", "relationships.tsv"]
DESCENT = {"child_a": ("father_a", "mother_a"), "child_b": ("father_a", "mother_a"), "father_a": ("gf_a", "gm_a"),
           "mother_a": (".", "gm_m")}


def check(name, ok, detail=""):
    print(("[ OK ] " if ok else "[FAIL] ") + name + ("" if ok else "  %s" % (detail,)))
    if not ok:
        FAILS.append(name)


def control(name, passed):
    print(("[ OK ] control: " if not passed else "[FAIL] control: ") + name +
          (" was caught" if not passed else " was NOT caught"))
    if passed:
        FAILS.append("control: " + name)


class FakeApp:
    """The plugin interface of bin/owngenome_live.py's App, as docs/live-site.md describes it."""

    def __init__(self, root):
        self.config = types.SimpleNamespace(root=root)
        self.routes, self.hooks, self.kinds, self.logs = [], [], [], []

    def route(self, method, pattern, fn, auth=True):
        self.routes.append((method, re.compile(pattern), fn, auth))

    def state_hook(self, name, fn):
        self.hooks.append((name, fn))

    def run_kind(self, fn):
        self.kinds.append(fn)

    def run_file(self, glob):
        raise AssertionError("the family plugin serves its files through its own route")

    def log(self, action, run="-", section="-", status="-"):
        self.logs.append((action, run, section, status))

    def runs(self):
        out = {}
        for fn in self.kinds:
            for r in fn(self):
                out[r["name"]] = r
        return out

    def get(self, path):
        req = FakeReq(self)
        for m, rx, fn, auth in self.routes:
            mt = rx.fullmatch(path)
            if mt and m == "GET":
                fn(req, mt)
                return req.code, req.ctype, req.body, auth
        return 404, None, b"", None


class FakeReq:
    def __init__(self, app):
        self.app, self.code, self.ctype, self.body = app, None, None, b""

    def run_dir(self, name):
        r = self.app.runs().get(name)
        if r:
            return r["dir"]
        d = os.path.join(self.app.config.root, name)
        return d if os.path.isdir(d) else None  # a sample run, as bin/owngenome_live.py lists one

    def send_body(self, code, ctype, body, headers=()):
        self.code, self.ctype, self.body = code, ctype, body.encode() if isinstance(body, str) else body

    def send_json(self, obj, code=200):
        self.send_body(code, "application/json", json.dumps(obj))


def descents(svg):
    out = {}
    for g in ET.fromstring(svg).iter(NS + "g"):
        if g.get("class") == "descent":
            for c in g.get("data-children").split():
                out[c] = (g.get("data-father"), g.get("data-mother"))
    return out


def state_checks(app, root):
    hook = dict(app.hooks)["family"]
    s = hook(os.path.join(root, "family_a"), {}, [])
    res = [("family_a's children with their parents and maps", s.get("children") == CHILDREN, s.get("children")),
           ("family_a's members: seven, gm_m not sequenced",
            [(m["sample"], m["sequenced"]) for m in s.get("members", [])] ==
            [("child_a", True), ("child_b", True), ("father_a", True), ("mother_a", True), ("gf_a", True),
             ("gm_a", True), ("gm_m", False)], s.get("members"))]
    return s, res


def route_checks(app):
    code, ctype, body, auth = app.get("/family/family_a/pedigree.svg")
    ok = code == 200 and ctype == "image/svg+xml" and auth is True
    return [("the pedigree is drawn from the declared parents, behind the login",
             ok and descents(body) == DESCENT, (code, ctype, auth, ok and descents(body)))]


def main():
    work = tempfile.mkdtemp(prefix="pgp-live-family-", dir=os.environ.get("TMPDIR"))
    try:
        root = os.path.join(work, "genome")
        shutil.copytree(FIX, root, symlinks=True)
        os.makedirs(os.path.join(root, "loner", "logs"))
        app = FakeApp(root)
        owngenome_live_family.register(app)

        # 1. the runs
        runs = app.runs()
        check("family_a and family_b are runs of kind family; child_a is not",
              sorted(runs) == ["family_a", "family_b"] and {r["kind"] for r in runs.values()} == {"family"}, runs)

        # 2. the family's state
        s, res = state_checks(app, root)
        for n, ok, d in res:
            check(n, ok, d)
        check("the relationship check as step 38 wrote it",
              [(r["a"], r["b"], r["verdict"]) for r in s["relationships"]] ==
              [("child_a", "father_a", "ok"), ("child_a", "child_b", "ok")], s["relationships"])
        check("the shared segments as step 38 wrote them",
              s["ibd_pairs"] == [{"a": "child_a", "b": "child_b", "declared": "full siblings", "shared_cm": "2550.0",
                                  "degree": "first (full siblings)"}], s["ibd_pairs"])
        check("the couple screen counted from step 39's table (1 shared, 1 X-linked), the grandparents not screened",
              s["couples"] == [{"a": "father_a", "b": "mother_a",
                                "report": "/family/family_a/couple_father_a_mother_a.md",
                                "table": "/family/family_a/couple_father_a_mother_a.tsv", "screened": True,
                                "shared": 1, "x_linked": 1},
                               {"a": "gf_a", "b": "gm_a", "report": None, "table": None, "screened": False,
                                "shared": 0, "x_linked": 0}], s["couples"])
        check("the plan rows", [r["analysis"] for r in s["plan"]] == ["joint_call", "couple_carrier"], s["plan"])
        check("exactly the result files that exist are listed", s["files"] == FILES, s["files"])
        check("family_b (pedigree.tsv only): a couple not screened yet",
              dict(app.hooks)["family"](os.path.join(root, "family_b"), {}, [])["couples"] ==
              [{"a": "partner_x", "b": "partner_y", "report": None, "table": None, "screened": False, "shared": 0,
                "x_linked": 0}])

        # 3. a sample's page
        hook = dict(app.hooks)["family"]
        check("child_a's page links to family_a", hook(os.path.join(root, "child_a"), {}, []) ==
              {"kind": "member", "families": [{"name": "family_a", "url": "/runs/family_a/"}]})
        check("a sample in no family gets null", hook(os.path.join(root, "loner"), {}, []) is None)

        # 4. the route
        for n, ok, d in route_checks(app):
            check(n, ok, d)
        code, ctype, body, _ = app.get("/family/family_a/couple_father_a_mother_a.md")
        with open(os.path.join(FIX, "family_a", "couple_father_a_mother_a.md"), "rb") as f:
            check("a listed file is served as written", code == 200 and body == f.read() and ctype.startswith("text/plain"),
                  (code, ctype))
        code, ctype, body, _ = app.get("/family/family_a/child_a_segments.svg")
        check("a child's map is served as an SVG", code == 200 and ctype == "image/svg+xml", (code, ctype))
        for path, why in (("/family/family_a/notes.txt", "a file no step writes"),
                          ("/family/family_a/..", "a name of dots"),
                          ("/family/child_a/pedigree.svg", "a sample's folder"),
                          ("/family/nobody/pedigree.svg", "an unknown family")):
            code = app.get(path)[0]
            check("404 for %s" % why, code == 404, code)
        rep = os.path.join(root, "family_a", "ibd_report.md")
        os.replace(rep, rep + ".moved")
        os.symlink(os.path.join(REPO, "README.md"), rep)
        code = app.get("/family/family_a/ibd_report.md")[0]
        listed = "ibd_report.md" in hook(os.path.join(root, "family_a"), {}, [])["files"]
        check("404 for a result file that is a symlink out of the folder, and it is not listed",
              code == 404 and not listed, (code, listed))
        os.remove(rep)
        os.replace(rep + ".moved", rep)
        check("every read is in the access log", ("family", "family_a", "pedigree.svg", 200) in app.logs and
              ("family", "family_a", "notes.txt", 404) in app.logs, app.logs)
        empty = os.path.join(root, "family_empty")
        os.makedirs(empty)
        with open(os.path.join(FIX, "family_a", "members.tsv")) as f, \
                open(os.path.join(empty, "members.tsv"), "w") as g:
            g.write(f.readline())
        es = hook(empty, {}, [])
        code = app.get("/family/family_empty/pedigree.svg")[0]
        check("a members file with a header and no rows: an error on the page, 404 for the pedigree",
              "error" in (es or {}) and code == 404, (es, code))
        shutil.rmtree(empty)

        # Planted mistake: child_a's father swapped for its grandfather in the copy.
        bad = os.path.join(work, "swapped")
        shutil.copytree(FIX, bad)
        p = os.path.join(bad, "family_a", "members.tsv")
        with open(p) as f:
            text = f.read()
        with open(p, "w") as f:
            f.write(text.replace("child_a\tfemale\tfather_a\tmother_a", "child_a\tfemale\tgf_a\tmother_a"))
        app2 = FakeApp(bad)
        owngenome_live_family.register(app2)
        control("a swapped parent in the state", state_checks(app2, bad)[1][0][1])
        control("a swapped parent in the drawing", route_checks(app2)[0][1])

        # 5. the page script
        with open(JS) as f:
            js = f.read()
        urls = re.findall(r"(?:https?:)?//[A-Za-z0-9.-]+\.[a-z]{2,}", js)
        check("family.js loads nothing from outside, sets no innerHTML, has no eval",
              not urls and "innerHTML" not in js and not re.search(r"\beval\s*\(|new\s+Function\s*\(", js), urls)
        node = shutil.which("node")
        if node:
            r = subprocess.run([node, "--check", JS], capture_output=True, text=True, timeout=60)
            check("node --check reads family.js", r.returncode == 0, r.stderr)
        elif os.environ.get("CI"):
            check("node is on the CI runner", False)
        else:
            print("[SKIP] no node to syntax-check family.js")

        # 6. the real server
        real_server(root)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print("%d check(s) failed: %s" % (len(FAILS), ", ".join(FAILS)))
        return 1
    print("all checks passed")
    return 0


def real_server(root):
    import owngenome_live  # on main: no skip, so a renamed server fails here instead of passing quietly
    state = os.path.join(os.path.dirname(root), "state")
    live, web, ingest = owngenome_live.start(["--root", root, "--state-dir", state, "--port", "0", "--ingest-port", "0",
                                              "--quiet"])
    threading.Thread(target=web.serve_forever, daemon=True).start()
    port = web.server_address[1]

    def get(path):
        c = http.client.HTTPConnection("127.0.0.1", port, timeout=30)
        c.request("GET", path)
        r = c.getresponse()
        return r.status, r.getheader("Content-Type"), r.read()

    try:
        code, _, body = get("/")
        check("the run list shows family_a as a family", code == 200 and b"family_a" in body and b"family" in body, code)
        code, _, body = get("/runs/family_a/state.json")
        fam = json.loads(body).get("plugins", {}).get("family") if code == 200 else None
        check("state.json carries plugins.family for family_a",
              fam and fam["kind"] == "family" and fam["children"] == CHILDREN, (code, fam))
        code, _, body = get("/runs/child_a/state.json")
        fam = json.loads(body).get("plugins", {}).get("family") if code == 200 else None
        check("state.json of child_a links to family_a", fam and fam["families"][0]["name"] == "family_a", (code, fam))
        code, ctype, body = get("/family/family_a/pedigree.svg")
        check("the real server draws the pedigree", code == 200 and ctype == "image/svg+xml"
              and descents(body) == DESCENT, (code, ctype))
        code, _, _ = get("/family/family_a/notes.txt")
        check("the real server does not serve notes.txt", code == 404, code)
        code, _, body = get("/runs/family_a/")
        check("the family page loads family.js", code == 200 and b"/assets/plugins/family.js" in body, code)
        browser_checks(port, os.path.dirname(root))
    finally:
        web.shutdown()
        ingest.server_close()
        web.server_close()


def browser_checks(port, work):
    sys.path.insert(0, os.path.join(REPO, "tests"))
    import test_owngenome_live as tpl
    chrome = tpl.find_chrome()
    if not chrome:
        ci = os.environ.get("CI") == "true"
        print(("[FAIL] " if ci else "[SKIP] ") + "no Chrome or Chromium to load the family page in")
        if ci:
            FAILS.append("no Chrome in CI")
        return
    b = tpl.Chrome(chrome, os.path.join(work, "chrome"))
    try:
        base = "http://127.0.0.1:%d" % port
        b.go(base + "/runs/family_a/")
        card = b.wait("var c = document.querySelector('section[data-plugin=family]'); "
                      "c && !c.hidden && c.querySelector('img.pedigree') && c.textContent")
        check("browser: the family card shows on family_a's page", card and "Couple carrier screen" in card, card)
        loaded = b.wait("var i = document.querySelector('img.pedigree'); i && i.complete && i.naturalWidth > 0 && "
                        "i.naturalWidth", 10)
        check("browser: the pedigree image loads", loaded, loaded)
        couple = b.eval("var p = document.querySelector('p[data-couple=\\'father_a mother_a\\']'); p && p.textContent")
        check("browser: the couple line counts 1 shared and 1 X-linked and links the report",
              couple and "1 gene(s) where both carry a finding, 1 X-linked" in couple
              and b.eval("!!document.querySelector('p[data-couple=\\'father_a mother_a\\'] a')"), couple)
        b.eval("document.querySelector('details[data-child=child_a]').open = true; true")
        seg = b.wait("var i = document.querySelector('details[data-child=child_a] img'); "
                     "i && i.complete && i.naturalWidth > 0 && i.getAttribute('src')", 10)
        check("browser: child_a's chromosome map loads", seg == "/family/family_a/child_a_segments.svg", seg)
        blocked = b.eval("document.documentElement.dataset.pgpBlocked || ''")
        check("browser: the policy blocked nothing on the family page", not blocked, blocked)
        b.go(base + "/runs/child_a/")
        href = b.wait("var a = document.querySelector('section[data-plugin=family] a'); a && a.getAttribute('href')")
        check("browser: child_a's page links to its family", href == "/runs/family_a/", href)
    finally:
        b.close()


if __name__ == "__main__":
    sys.exit(main())
