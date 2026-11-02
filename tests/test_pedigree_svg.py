#!/usr/bin/env python3
"""bin/pedigree_svg.py on tests/fixtures/families_live.

Checks, reading the drawing back from its data- attributes and coordinates:
  1. family_a (members.tsv, three generations): every person is drawn once,
     with the sex and sequenced state of the file; each generation is one
     row, parents above their children;
  2. each line of descent starts at its two parents (the middle of their
     partner line) and ends on top of each child it names, and the parents
     it names are the ones this test expects (written here, not read from
     the fixture);
  3. family_b (pedigree.tsv only): the couple is drawn joined, both filled;
  4. a person whose sex is not known is a diamond; the command line writes
     the same SVG and exits 2 on a folder without a pedigree.
Checks 1 and 2 are run again on a copy of the fixture with two parents
swapped and must fail there: a check that cannot fail is no check.

Run: python3 tests/test_pedigree_svg.py
"""
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import pedigree_svg  # noqa: E402

FIX = os.path.join(REPO, "tests", "fixtures", "families_live")
NS = "{http://www.w3.org/2000/svg}"
FAILS = []

# What family_a declares, written out by hand: child -> (father, mother).
EXPECTED = {"child_a": ("father_a", "mother_a"), "child_b": ("father_a", "mother_a"),
            "father_a": ("gf_a", "gm_a"), "mother_a": (".", "gm_m")}
SEX = {"child_a": "female", "child_b": "male", "father_a": "male", "mother_a": "female", "gf_a": "male",
       "gm_a": "female", "gm_m": "female"}
SEQUENCED = {"gm_m": "no"}


def check(name, ok, detail=""):
    print(("[ OK ] " if ok else "[FAIL] ") + name + ("" if ok else "  %s" % (detail,)))
    if not ok:
        FAILS.append(name)


def control(name, passed):
    """A check run on a planted mistake: it must fail there."""
    print(("[ OK ] control: " if not passed else "[FAIL] control: ") + name +
          (" was caught" if not passed else " was NOT caught"))
    if passed:
        FAILS.append("control: " + name)


def parse(svg):
    root = ET.fromstring(svg)
    people = {}
    for g in root.iter(NS + "g"):
        if g.get("class") == "person":
            shape = [c for c in g if c.tag in (NS + "rect", NS + "circle", NS + "polygon")][0]
            if shape.tag == NS + "rect":
                cx, cy = float(shape.get("x")) + 18, float(shape.get("y")) + 18
            elif shape.tag == NS + "circle":
                cx, cy = float(shape.get("cx")), float(shape.get("cy"))
            else:
                pts = [tuple(map(float, p.split(","))) for p in shape.get("points").split()]
                cx, cy = pts[0][0], (pts[0][1] + pts[2][1]) / 2
            people.setdefault(g.get("data-sample"), []).append(
                {"sex": g.get("data-sex"), "gen": int(g.get("data-gen")), "seq": g.get("data-sequenced"),
                 "shape": shape.tag[len(NS):], "x": cx, "y": cy, "fill": shape.get("fill")})
    descents = []
    for g in root.iter(NS + "g"):
        if g.get("class") == "descent":
            lines = [{k: float(ln.get(k)) for k in ("x1", "y1", "x2", "y2")} for ln in g.iter(NS + "line")]
            descents.append({"father": g.get("data-father"), "mother": g.get("data-mother"),
                             "children": g.get("data-children").split(), "lines": lines})
    couples = [(ln.get("data-a"), ln.get("data-b")) for ln in root.iter(NS + "line") if ln.get("class") == "couple"]
    return people, descents, couples


def people_ok(people):
    bad = [s for s in SEX if len(people.get(s, [])) != 1 or people[s][0]["sex"] != SEX[s]
           or people[s][0]["seq"] != SEQUENCED.get(s, "yes")
           or people[s][0]["shape"] != {"male": "rect", "female": "circle"}[SEX[s]]]
    return not bad and set(people) == set(SEX), bad or sorted(set(people) ^ set(SEX))


def rows_ok(people):
    bad = []
    for child, (f, m) in EXPECTED.items():
        for p in (f, m):
            if p != "." and not people[p][0]["y"] < people[child][0]["y"]:
                bad.append((p, child))
    ys = {}
    for s, v in people.items():
        ys.setdefault(v[0]["gen"], set()).add(v[0]["y"])
    return not bad and all(len(y) == 1 for y in ys.values()), (bad, ys)


def descent_ok(people, descents):
    """Every child of EXPECTED hangs from a descent that names its expected parents, starts at
    their middle and ends on the child's symbol."""
    bad, seen = [], {}
    for d in descents:
        for c in d["children"]:
            seen[c] = (d["father"], d["mother"])
        ps = [people[p][0] for p in (d["father"], d["mother"]) if p != "."]
        ox = sum(p["x"] for p in ps) / len(ps)
        first = d["lines"][0]
        if abs(first["x1"] - ox) > 0.5 or abs(first["x2"] - ox) > 0.5:
            bad.append(("start", d["father"], d["mother"]))
        for c in d["children"]:
            p = people[c][0]
            if not any(abs(ln["x1"] - p["x"]) < 0.5 and abs(ln["y2"] - (p["y"] - 18)) < 0.5 for ln in d["lines"]):
                bad.append(("end", c))
    if seen != EXPECTED:
        bad.append(("parents", seen))
    return not bad, bad


def run_checks(folder):
    order, members, _ = pedigree_svg.load(folder)
    people, descents, couples = parse(pedigree_svg.draw(order, members, "family_a"))
    return [people_ok(people), rows_ok(people), descent_ok(people, descents)], couples


def main():
    work = tempfile.mkdtemp(prefix="pedigree-svg-", dir=os.environ.get("TMPDIR"))
    try:
        results, couples = run_checks(os.path.join(FIX, "family_a"))
        names = ["every person drawn once with the file's sex and sequenced state",
                 "one row per generation, parents above their children",
                 "each line of descent runs from the expected parents to each child"]
        for n, (ok, detail) in zip(names, results):
            check(n, ok, detail)
        check("partners joined: father_a with mother_a, gf_a with gm_a",
              {frozenset(c) for c in couples} == {frozenset(("father_a", "mother_a")), frozenset(("gf_a", "gm_a"))},
              couples)

        # Planted mistake: child_a's parents swapped with father_a's parents.
        bad = os.path.join(work, "swapped")
        os.makedirs(bad)
        with open(os.path.join(FIX, "family_a", "members.tsv")) as f:
            text = f.read()
        lines = text.split("\n")
        for i, ln in enumerate(lines):
            c = ln.split("\t")
            if c[0] == "child_a":
                c[2], c[3] = "gf_a", "gm_a"
            elif c[0] == "father_a":
                c[2], c[3] = "", ""
            lines[i] = "\t".join(c)
        with open(os.path.join(bad, "members.tsv"), "w") as f:
            f.write("\n".join(lines))
        results, _ = run_checks(bad)
        control("child_a drawn under its grandparents", results[2][0] and results[1][0])

        order, members, _ = pedigree_svg.load(os.path.join(FIX, "family_b"))
        people, descents, couples = parse(pedigree_svg.draw(order, members, "family_b"))
        check("family_b from pedigree.tsv: a joined, filled couple and no descent",
              couples == [("partner_x", "partner_y")] and not descents
              and [people[s][0]["seq"] for s in ("partner_x", "partner_y")] == ["yes", "yes"], (couples, people))

        unk = os.path.join(work, "unknown")
        os.makedirs(unk)
        with open(os.path.join(unk, "pedigree.tsv"), "w") as f:
            f.write("sample\tsex\tfather\tmother\tpartner\nkid_u\tunknown\tdad_u\tmum_u\t.\n")
        order, members, _ = pedigree_svg.load(unk)
        people, descents, _ = parse(pedigree_svg.draw(order, members))
        check("a sex not known is a diamond; named parents without a row are drawn not sequenced",
              people["kid_u"][0]["shape"] == "polygon" and people["dad_u"][0]["seq"] == "no"
              and descents and descents[0]["children"] == ["kid_u"], people)

        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "pedigree_svg.py"),
                            os.path.join(FIX, "family_a")], capture_output=True, text=True, timeout=60)
        order, members, _ = pedigree_svg.load(os.path.join(FIX, "family_a"))
        check("the command line writes the same SVG", r.returncode == 0 and r.stdout ==
              pedigree_svg.draw(order, members, "family_a"), r.stderr)
        empty = os.path.join(work, "empty")
        os.makedirs(empty)
        r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "pedigree_svg.py"), empty],
                           capture_output=True, text=True, timeout=60)
        check("a folder without a pedigree exits 2 and says why", r.returncode == 2 and "members.tsv" in r.stderr,
              (r.returncode, r.stderr))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if FAILS:
        print("%d check(s) failed: %s" % (len(FAILS), ", ".join(FAILS)))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
