#!/usr/bin/env python3
"""pedigree_svg.py: draw a family's pedigree as an SVG.

  pedigree_svg.py FAMILY_DIR_OR_FILE [--title TEXT] > pedigree.svg

Reads members.tsv (step 38's table of every person the pedigree names, with
in_pedigree and sequenced) or, when there is none, pedigree.tsv (the rows a
person wrote; parents and partners named without a row are added, as step 38
adds them). A folder argument reads FAMILY_DIR/members.tsv, else
FAMILY_DIR/pedigree.tsv.

The drawing uses the usual pedigree symbols: a square for a male, a circle for
a female, a diamond when the sex is not known. A filled symbol is a sequenced
person (members.tsv's sequenced column; from pedigree.tsv, a person with a row
of their own). A horizontal line joins two partners, or the two parents of a
child; a line drops from it to a bar over their children. Each generation is
one row, parents above their children.

It only draws what the file says. It reads no genotype and checks no
relationship: step 38's relationship check does that.

Every symbol is a <g class="person" data-sample data-sex data-gen
data-sequenced>, every partner line a <line class="couple" data-a data-b>,
and every line of descent a <g class="descent" data-father data-mother
data-children>, so a test can read the drawing back.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SIZE = 36          # a symbol's width and height
STEP_X = 120       # from one person to the next in a row
STEP_Y = 130       # from one generation to the next
MARGIN = 40
TOP = 50           # room for the title
BG, FG, LINE, FILL = "#16181c", "#e6e6e6", "#c0c4cc", "#7aa2f7"
SEXES = ("male", "female")


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def blank(v):
    return v is None or v.strip() in ("", ".", "-9", "0", "NA", "na")


def read_members_tsv(path):
    """(order, members) from members.tsv; a member is {sex, father, mother,
    partner, sequenced}."""
    order, members = [], {}
    with open(path) as f:
        head = f.readline().rstrip("\r\n").lstrip("#").split("\t")
        if "sample" not in head:
            raise ValueError("%s: the first line must be a header with a 'sample' column" % path)
        for line in f:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith("#"):
                continue
            vals = line.split("\t")
            r = dict(zip(head, vals + [""] * (len(head) - len(vals))))
            s = r["sample"].strip()
            if blank(s) or s in members:
                continue
            sex = r.get("sex", "").strip().lower()
            members[s] = {"sex": sex if sex in SEXES else "unknown",
                          "father": None if blank(r.get("father")) else r["father"].strip(),
                          "mother": None if blank(r.get("mother")) else r["mother"].strip(),
                          "partner": None if blank(r.get("partner")) else r["partner"].strip(),
                          "sequenced": r.get("sequenced", "").strip() == "yes"}
            order.append(s)
    return order, members


def read_pedigree_tsv(path):
    """(order, members) from a pedigree file, with step 38's own reader."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import family_report
    rows, errors = family_report.read_pedigree(path)
    members, more = family_report.build_members(rows, {}, path) if not errors else ({}, [])
    if errors or more:
        raise ValueError("; ".join(errors + more))
    out = {}
    for s, m in members.items():
        out[s] = {"sex": m.get("sex") if m.get("sex") in SEXES else "unknown", "father": m.get("father"),
                  "mother": m.get("mother"), "partner": m.get("partner"), "sequenced": bool(m.get("row"))}
    return list(members), out


def load(path):
    """(order, members, source file) of a family folder or a members/pedigree file."""
    if os.path.isdir(path):
        for name in ("members.tsv", "pedigree.tsv"):
            p = os.path.join(path, name)
            if os.path.isfile(p):
                return load(p)
        raise ValueError("%s holds no members.tsv or pedigree.tsv" % path)
    with open(path) as f:
        head = f.readline().rstrip("\r\n").lstrip("#").split("\t")
    if "sequenced" in head:
        order, members = read_members_tsv(path)
    else:
        order, members = read_pedigree_tsv(path)
    return order, members, path


def couples(order, members):
    """Pairs drawn joined: partners, and the two parents of a child, in file order."""
    seen, out = set(), []
    for s in order:
        m = members[s]
        for a, b in ((s, m["partner"]), (m["father"], m["mother"])):
            if a and b and a in members and b in members:
                key = frozenset((a, b))
                if key not in seen:
                    seen.add(key)
                    out.append((a, b))
    return out


def generations(order, members, pairs):
    """Generation of each person: 0 with no parent named, else one below the
    lower parent. A partner with no parent named moves down to the other
    partner's row; two partners who both have parents keep their own rows
    (a line then joins them across rows)."""
    gen = {s: 0 for s in order}

    def rooted(s):
        return any(members[s][k] in gen for k in ("father", "mother"))

    for _ in range(len(order) + 2):
        changed = False
        for s in order:
            ps = [members[s][k] for k in ("father", "mother") if members[s][k] in gen]
            g = max([gen[p] + 1 for p in ps] + [gen[s]])
            if g != gen[s]:
                gen[s], changed = g, True
        for a, b in pairs:
            for x, y in ((a, b), (b, a)):
                if not rooted(x) and gen[x] < gen[y]:
                    gen[x], changed = gen[y], True
        if not changed:
            break
    return gen


def layout(order, members):
    """{person: (x, y)} of each symbol's centre, plus the partner pairs."""
    pairs = couples(order, members)
    gen = generations(order, members, pairs)
    partners = {}
    for a, b in pairs:
        partners.setdefault(a, []).append(b)
        partners.setdefault(b, []).append(a)
    rank = {s: i for i, s in enumerate(order)}
    rows, x_of = {}, {}
    for g in sorted(set(gen.values())):
        people = [s for s in order if gen[s] == g]

        def key(s):
            ps = [x_of[members[s][k]] for k in ("father", "mother") if members[s][k] in x_of]
            return (sum(ps) / len(ps) if ps else None, rank[s])

        # People with a parent drawn go under their parents; one without (a partner who married in,
        # or a whole first generation) goes beside the partner already placed, else in file order.
        anchored = sorted([s for s in people if key(s)[0] is not None], key=key)
        row, placed = [], set()

        def put(s):
            if s in placed:
                return
            placed.add(s)
            row.append(s)
            for p in partners.get(s, []):
                if gen.get(p) == g and key(p)[0] is None:
                    put(p)

        for s in anchored + [s for s in people if key(s)[0] is None]:
            put(s)
        rows[g] = row
        for i, s in enumerate(row):
            x_of[s] = i
    widest = max(len(r) for r in rows.values())
    pos = {}
    for g, row in rows.items():
        shift = (widest - len(row)) * STEP_X / 2
        for i, s in enumerate(row):
            pos[s] = (MARGIN + SIZE / 2 + shift + i * STEP_X, TOP + SIZE / 2 + g * STEP_Y)
    return pos, gen, pairs, widest, max(rows) + 1


def symbol(s, m, x, y, g):
    h = SIZE / 2
    fill = FILL if m["sequenced"] else BG
    if m["sex"] == "male":
        shape = '<rect x="%.1f" y="%.1f" width="%d" height="%d"' % (x - h, y - h, SIZE, SIZE)
    elif m["sex"] == "female":
        shape = '<circle cx="%.1f" cy="%.1f" r="%.1f"' % (x, y, h)
    else:
        shape = '<polygon points="%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f"' % (x, y - h, x + h, y, x, y + h, x - h, y)
    return ('<g class="person" data-sample="%s" data-sex="%s" data-gen="%d" data-sequenced="%s">'
            '%s fill="%s" stroke="%s" stroke-width="2"/>'
            '<text x="%.1f" y="%.1f" fill="%s" text-anchor="middle">%s</text>'
            '<title>%s: %s, %s</title></g>'
            % (esc(s), m["sex"], g, "yes" if m["sequenced"] else "no", shape, fill, LINE,
               x, y + h + 16, FG, esc(s), esc(s), m["sex"], "sequenced" if m["sequenced"] else "not sequenced"))


def draw(order, members, title="Pedigree"):
    if not order:
        raise ValueError("no person to draw")
    pos, gen, pairs, widest, ngen = layout(order, members)
    width = int(2 * MARGIN + SIZE + (widest - 1) * STEP_X)
    legend_y = TOP + SIZE + (ngen - 1) * STEP_Y + 50
    height = legend_y + 30
    L = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
         'font-family="sans-serif" font-size="12">' % (width, height, width, height),
         '<rect width="100%%" height="100%%" fill="%s"/>' % BG,
         '<text x="10" y="22" fill="%s" font-size="14">%s</text>' % (FG, esc(title))]
    h = SIZE / 2
    for a, b in pairs:
        (xa, ya), (xb, yb) = pos[a], pos[b]
        if xa > xb:
            (xa, ya), (xb, yb) = (xb, yb), (xa, ya)
        L.append('<line class="couple" data-a="%s" data-b="%s" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                 'stroke="%s" stroke-width="2"/>' % (esc(a), esc(b), xa + h, ya, xb - h, yb, LINE))
    # Children grouped by their pair of parents, in file order.
    sibships = {}
    for s in order:
        f, mo = members[s]["father"], members[s]["mother"]
        f, mo = (f if f in pos else None), (mo if mo in pos else None)
        if f or mo:
            sibships.setdefault((f, mo), []).append(s)
    for (f, mo), kids in sibships.items():
        if f and mo:
            ox, oy = (pos[f][0] + pos[mo][0]) / 2, (pos[f][1] + pos[mo][1]) / 2
        else:
            ox, oy = pos[f or mo][0], pos[f or mo][1] + h
        top = min(pos[k][1] for k in kids) - h
        bar = top - 22
        xs = [pos[k][0] for k in kids] + [ox]
        parts = ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>' % (ox, oy, ox, bar),
                 '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>' % (min(xs), bar, max(xs), bar)]
        parts += ['<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>' % (pos[k][0], bar, pos[k][0], pos[k][1] - h)
                  for k in kids]
        L.append('<g class="descent" data-father="%s" data-mother="%s" data-children="%s" stroke="%s" '
                 'stroke-width="2">%s</g>' % (esc(f or "."), esc(mo or "."), esc(" ".join(kids)), LINE, "".join(parts)))
    for s in order:
        L.append(symbol(s, members[s], pos[s][0], pos[s][1], gen[s]))
    L.append('<text x="10" y="%d" fill="%s">Square: male. Circle: female. Diamond: sex not known. Filled: '
             'sequenced. A line joins partners; children hang from it.</text>' % (legend_y, FG))
    L.append("</svg>")
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", help="a family folder, members.tsv or pedigree.tsv")
    ap.add_argument("--title", help="default: the folder's name")
    a = ap.parse_args(argv)
    try:
        order, members, src = load(a.path)
        title = a.title or os.path.basename(os.path.dirname(os.path.abspath(src))) or "Pedigree"
        sys.stdout.write(draw(order, members, title))
    except (OSError, ValueError) as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
