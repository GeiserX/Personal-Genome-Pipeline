#!/usr/bin/env python3
"""family_segments.py: the chromosome map of step 38 (scripts/38-family.sh)
and of the Nextflow FAMILY workflow. family_report.py report calls it, so
both write the same files; it also runs on its own over a family folder.
Standard library only.

  family_segments.py --members D/members.tsv --vcf FAMILY.vcf.gz [--par-bed par_grch38.bed] --out-dir D

For each sequenced child with a sequenced parent, the stretches of each
chromosome the child got from that parent. A site tells something when the
parent is heterozygous and the genotypes give the allele the parent passed
on. When one of that parent's own parents is sequenced too, the site also
tells which of the parent's two copies it was (the one from the
grandfather or the one from the grandmother), and a change from one to the
other along the chromosome is a crossover. With only one parent sequenced,
the child's other copy is mapped over the same stretches as coming from the
other side.

Stretches are smoothed: a run of sites that points to the other copy is
kept only with MIN_SITES sites over MIN_BP bases; a shorter run (a genotype
error, a gene conversion) is absorbed and counted as discordant. A stretch
with no site for MAX_GAP bases (a centromere, a gap of the reference) is
split there, since nothing is known inside it. Autosomes, and chrX outside
the pseudoautosomal regions from the mother only: a father passes his one X
whole to a daughter and none to a son.

Writes, per child:
  <child>_segments.tsv     child, side, parent, chrom, start, end, length_bp,
                           sites, discordant, ancestor, basis
  <child>_dnapainter.csv   the same segments in DNA Painter's own import
                           format (Chr,Start,End,cM,SNPs,Match,Confidence,
                           Group,Side,Notes,Colour,Ancestor); GRCh38 positions
  <child>_segments.svg     a map of every chromosome, the father's copy above
                           the mother's
"""
import argparse
import csv
import io
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_report as fr  # noqa: E402  the pedigree and genotype rules of step 38

# A run of sites pointing to the other grandparent's copy needs this much
# support to be a segment of its own (a crossover on each side of it).
MIN_SITES = 20
MIN_BP = 500_000
# No site for this many bases: the stretch is split, nothing is claimed inside.
MAX_GAP = 5_000_000
AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX"]
COLOURS = {"paternal": ("#4f8fd6", "#9cc3ef"), "maternal": ("#d65f8f", "#efa3c0")}
SIDE_ONLY = {"paternal": "#3b6fa8", "maternal": "#a84870"}
OTHER = "#6b6b6b"


def unique(vals):
    vals = set(vals)
    return next(iter(vals)) if len(vals) == 1 else None


class Link:
    """One child and one sequenced parent: the parent's informative sites."""

    def __init__(self, child, side, parent, members):
        self.child, self.side, self.parent = child, side, parent
        p = members[parent]
        self.gf = p["father"]
        self.gm = p["mother"]
        self.gf_seq = bool(self.gf and members[self.gf]["sequenced"])
        self.gm_seq = bool(self.gm and members[self.gm]["sequenced"])
        self.sites = {}   # chrom -> [pos] where the allele passed on is known
        self.hap = {}     # chrom -> [(pos, 0 from the grandfather | 1 from the grandmother)]

    @property
    def grandparents(self):
        return self.gf_seq or self.gm_seq

    def ancestor(self, label):
        name = self.gf if label == 0 else self.gm
        return name or "%s's %s" % (self.parent, "father" if label == 0 else "mother")


def gset(g):
    return set(g) if g is not None else None


def passed_on(child_gt, parent_gt, other_gt, parent_is_father):
    """The allele the parent passed on, when the genotypes give one."""
    f, m = (parent_gt, other_gt) if parent_is_father else (other_gt, parent_gt)
    opts = fr.resolve(child_gt, gset(f), gset(m))
    return unique(o[0] if parent_is_father else o[1] for o in opts) if opts else None


def collect(members, vcf, par):
    """The Link of every sequenced child and sequenced parent, filled from the joint VCF."""
    seq = sorted(s for s, m in members.items() if m["sequenced"])
    links = {}
    for s in seq:
        m = members[s]
        for side, p in (("paternal", m["father"]), ("maternal", m["mother"])):
            if p and members[p]["sequenced"]:
                links.setdefault(s, []).append(Link(s, side, p, members))
    if not links:
        return members, links
    idx = None
    with fr.open_text(vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")[9:]
                idx = {s: 9 + i for i, s in enumerate(cols)}
                lost = [s for s in seq if s not in idx]
                if lost:
                    fr.die(["%s has no column for %s" % (vcf, ", ".join(lost))])
                continue
            c = line.rstrip("\n").split("\t")
            chrom = c[0]
            if chrom not in CHROMS:
                continue
            pos = int(c[1])
            if chrom == "chrX" and fr.in_par(par, chrom, pos):
                continue
            fmt = c[8].split(":")
            if "GT" not in fmt:
                continue
            gi = fmt.index("GT")
            cache = {}

            def g(s):
                if s is None or not members.get(s, {}).get("sequenced"):
                    return None
                if s not in cache:
                    v = c[idx[s]].split(":")
                    cache[s] = fr.parse_gt(v[gi] if gi < len(v) else ".")[0]
                return cache[s]

            for child, ls in links.items():
                cg = g(child)
                if cg is None:
                    continue
                for ln in ls:
                    pg = g(ln.parent)
                    if pg is None or len(pg) != 2 or pg[0] == pg[1]:
                        continue   # only a heterozygous parent tells which copy was passed on
                    if chrom == "chrX":
                        if ln.side != "maternal":
                            continue
                        if members[child]["sex"] == "male":
                            if len(set(cg)) != 1:
                                continue
                            t = cg[0]
                        else:
                            fg = g(members[child]["father"])
                            if fg is None or len(set(fg)) != 1 or len(cg) != 2:
                                continue
                            t = passed_on(cg, pg, (fg[0], fg[0]), False)
                    else:
                        if len(cg) != 2:
                            continue
                        other = members[child]["mother" if ln.side == "paternal" else "father"]
                        t = passed_on(cg, pg, g(other), ln.side == "paternal")
                    if t is None or t not in pg:
                        continue
                    ln.sites.setdefault(chrom, []).append(pos)
                    if not ln.grandparents:
                        continue
                    # which of the parent's copies holds t: the parent's own origin here
                    if chrom == "chrX":
                        # the mother's X: one from her father (his only X), one from her mother
                        gfg, gmg = g(ln.gf), g(ln.gm)
                        if gfg is not None and len(set(gfg)) != 1:
                            continue
                        po = passed_on(pg, gmg, (gfg[0], gfg[0]) if gfg else None, False)
                        if po is None:
                            continue
                        label = 1 if t == po else 0
                    else:
                        opts = fr.resolve(pg, gset(g(ln.gf)), gset(g(ln.gm)))
                        from_gf = unique(o[0] for o in opts) if opts else None
                        if from_gf is None:
                            continue
                        label = 0 if t == from_gf else 1
                    ln.hap.setdefault(chrom, []).append((pos, label))
    return members, links


def blocks(positions):
    """Split sorted positions where no site falls for MAX_GAP bases."""
    out, cur = [], []
    for p in positions:
        if cur and p - cur[-1] > MAX_GAP:
            out.append(cur)
            cur = []
        cur.append(p)
    if cur:
        out.append(cur)
    return out


def smooth(sites):
    """Runs [label, start, end, agree, discordant] of one block of (pos, label)."""
    runs = []
    for pos, label in sites:
        if runs and runs[-1][0] == label:
            runs[-1][2] = pos
            runs[-1][3] += 1
        else:
            runs.append([label, pos, pos, 1, 0])

    def weak(r):
        return r[3] < MIN_SITES or r[2] - r[1] + 1 < MIN_BP

    while len(runs) > 1:
        w = [i for i, r in enumerate(runs) if weak(r)]
        if not w:
            break
        i = min(w, key=lambda k: (runs[k][3], runs[k][2] - runs[k][1], k))
        lo, hi = max(i - 1, 0), min(i + 1, len(runs) - 1)
        keep = runs[lo] if lo != i else runs[hi]
        merged = [keep[0], runs[lo][1], runs[hi][2], 0, 0]
        for r in runs[lo:hi + 1]:
            if r is runs[i]:
                # its sites now point to the other copy; the ones it had absorbed agree again
                merged[3] += r[4]
                merged[4] += r[3]
            else:
                merged[3] += r[3]
                merged[4] += r[4]
        runs[lo:hi + 1] = [merged]
    return [r for r in runs if r[3] >= MIN_SITES]


def segments(members, links):
    """{child: [row dict]}, rows in chromosome order, father's side first."""
    order = {c: i for i, c in enumerate(CHROMS)}
    out = {}
    for child, ls in sorted(links.items()):
        rows = []
        for ln in ls:
            hap = ln.grandparents and sum(len(v) for v in ln.hap.values()) >= MIN_SITES
            for chrom in sorted(ln.sites, key=order.get):
                if hap and chrom in ln.hap:
                    sites = ln.hap[chrom]
                    for blk in blocks([p for p, _ in sites]):
                        lo, hi = blk[0], blk[-1]
                        part = [s for s in sites if lo <= s[0] <= hi]
                        for label, s, e, n, disc in smooth(part):
                            rows.append(row(child, ln.side, ln.parent, chrom, s, e, n, disc,
                                            ln.ancestor(label), "parent and grandparents"))
                else:
                    for blk in blocks(ln.sites[chrom]):
                        if len(blk) >= MIN_SITES:
                            rows.append(row(child, ln.side, ln.parent, chrom, blk[0], blk[-1], len(blk), 0,
                                            ".", "parent"))
        # One parent sequenced: the other copy over the same stretches
        if len(ls) == 1:
            ln = ls[0]
            other_side = "maternal" if ln.side == "paternal" else "paternal"
            other = members[child]["mother" if ln.side == "paternal" else "father"]
            for chrom in sorted(ln.sites, key=order.get):
                if chrom == "chrX":
                    continue   # the X a child got from the father is not seen from the mother's sites
                for blk in blocks(ln.sites[chrom]):
                    if len(blk) >= MIN_SITES:
                        rows.append(row(child, other_side, other or ".", chrom, blk[0], blk[-1], len(blk), 0,
                                        ".", "other copy"))
        rows.sort(key=lambda r: (order[r["chrom"]], r["side"] != "paternal", r["start"]))
        out[child] = rows
    return out


def row(child, side, parent, chrom, start, end, n, disc, ancestor, basis):
    return {"child": child, "side": side, "parent": parent, "chrom": chrom, "start": start, "end": end,
            "length_bp": end - start + 1, "sites": n, "discordant": disc, "ancestor": ancestor, "basis": basis}


TSV = ("child", "side", "parent", "chrom", "start", "end", "length_bp", "sites", "discordant", "ancestor", "basis")
DNAPAINTER = ("Chr", "Start", "End", "cM", "SNPs", "Match", "Confidence", "Group", "Side", "Notes", "Colour",
              "Ancestor")


def contig_lengths(vcf):
    out = {}
    with fr.open_text(vcf) as f:
        for line in f:
            if not line.startswith("##"):
                break
            m = re.match(r"##contig=<ID=([^,>]+),.*?length=(\d+)", line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def crossovers(rows, side):
    """Changes of ancestor between neighbouring segments of one chromosome."""
    n, prev = 0, None
    for r in rows:
        if r["side"] != side or r["basis"] != "parent and grandparents":
            continue
        if prev and prev["chrom"] == r["chrom"] and prev["ancestor"] != r["ancestor"]:
            n += 1
        prev = r
    return n


def write_files(rows_by_child, lengths, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for child, rows in rows_by_child.items():
        with open(os.path.join(out_dir, "%s_segments.tsv" % child), "w") as f:
            f.write("\t".join(TSV) + "\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in TSV) + "\n")
        buf = io.StringIO()
        w = csv.writer(buf, lineterminator="\n")
        w.writerow(DNAPAINTER)
        for r in rows:
            who = r["ancestor"] if r["ancestor"] != "." else r["parent"] if r["parent"] != "." else "%s side" % r["side"]
            note = "GRCh38 positions; %s; %d discordant sites" % (
                {"parent and grandparents": "the parent's copy, told by a sequenced grandparent",
                 "parent": "from the sequenced parent",
                 "other copy": "the copy not from the sequenced parent"}[r["basis"]], r["discordant"])
            w.writerow([r["chrom"][3:], r["start"], r["end"], "", r["sites"], "%s from %s" % (child, who), "",
                        r["parent"], r["side"], note, "", r["ancestor"] if r["ancestor"] != "." else ""])
        with open(os.path.join(out_dir, "%s_dnapainter.csv" % child), "w") as f:
            f.write(buf.getvalue())
        with open(os.path.join(out_dir, "%s_segments.svg" % child), "w") as f:
            f.write(svg(child, rows, lengths))


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def svg(child, rows, lengths):
    chroms = [c for c in CHROMS if c in lengths or any(r["chrom"] == c for r in rows)]
    longest = max([lengths.get(c, 0) for c in chroms] + [r["end"] for r in rows] + [1])
    width, left, bar, row_h = 900, 60, 9, 28
    scale = (width - left - 20) / longest
    colour = {}
    for r in rows:
        if r["basis"] == "other copy":
            colour[(r["side"], r["ancestor"])] = OTHER
        elif r["ancestor"] == ".":
            colour[(r["side"], ".")] = SIDE_ONLY[r["side"]]
        else:
            seen = [k for k in colour if k[0] == r["side"] and k[1] != "."]
            if (r["side"], r["ancestor"]) not in colour:
                colour[(r["side"], r["ancestor"])] = COLOURS[r["side"]][len(seen) % 2]
    height = 50 + row_h * len(chroms) + 20 * (len(colour) + 1)
    L = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">'
         % (width, height),
         '<rect width="100%" height="100%" fill="#16181c"/>',
         '<text x="10" y="20" fill="#e6e6e6" font-size="14">%s: which copy of each chromosome came from whom '
         '(GRCh38)</text>' % esc(child)]
    for i, c in enumerate(chroms):
        y = 40 + i * row_h
        L.append('<text x="10" y="%d" fill="#e6e6e6">%s</text>' % (y + bar + 3, esc(c[3:])))
        n = lengths.get(c, max([r["end"] for r in rows if r["chrom"] == c] + [1]))
        for k, side in enumerate(("paternal", "maternal")):
            L.append('<rect x="%d" y="%d" width="%.1f" height="%d" fill="#2a2d33"/>'
                     % (left, y + k * (bar + 1), n * scale, bar))
        for r in rows:
            if r["chrom"] != c:
                continue
            k = 0 if r["side"] == "paternal" else 1
            L.append('<rect class="seg" x="%.1f" y="%d" width="%.1f" height="%d" fill="%s"><title>%s %s:%d-%d, '
                     '%s, %d sites</title></rect>'
                     % (left + r["start"] * scale, y + k * (bar + 1), max((r["end"] - r["start"]) * scale, 1), bar,
                        colour[(r["side"], r["ancestor"])],
                        esc(r["side"]), esc(r["chrom"]), r["start"], r["end"],
                        esc(r["ancestor"] if r["ancestor"] != "." else r["parent"]), r["sites"]))
    y = 50 + row_h * len(chroms)
    L.append('<text x="10" y="%d" fill="#e6e6e6">Upper bar: the copy from the father. Lower bar: from the mother. '
             'Grey: not mapped.</text>' % y)
    for j, ((side, anc), col) in enumerate(sorted(colour.items())):
        label = {OTHER: "%s side, the copy not from the sequenced parent" % side}.get(col)
        if label is None:
            label = "%s side, %s" % (side, "from the parent" if anc == "." else "from " + anc)
        L.append('<rect x="10" y="%d" width="12" height="10" fill="%s"/><text x="28" y="%d" fill="#e6e6e6">%s</text>'
                 % (y + 10 + 20 * j, col, y + 19 + 20 * j, esc(label)))
    L.append("</svg>")
    return "\n".join(L) + "\n"


def markdown(members, rows_by_child):
    """The "Chromosome map" section of family_report.md."""
    L = ["", "## Chromosome map", ""]
    if not rows_by_child:
        L.append("No sequenced child has a sequenced parent, so there is no map.")
        return L
    L += ["For each child, the stretches of each chromosome that came from each parent, from the sites where "
          "that parent is heterozygous and the genotypes give the allele passed on. When a parent's own parent "
          "is sequenced too, each stretch also says which of the parent's two copies it is (the grandfather's or "
          "the grandmother's), and a change between them is a crossover. `<child>_segments.tsv` lists them, "
          "`<child>_dnapainter.csv` is the same list in DNA Painter's import format (GRCh38 positions), "
          "`<child>_segments.svg` draws it.", ""]
    for child, rows in rows_by_child.items():
        m = members[child]
        for side, who in (("paternal", "father"), ("maternal", "mother")):
            mine = [r for r in rows if r["side"] == side]
            if not mine:
                continue
            bp = sum(r["length_bp"] for r in mine)
            name = m[who] or "not named"
            bases = {r["basis"] for r in mine}
            if "other copy" in bases:
                L.append("- %s: from the %s (%s, not sequenced): %d stretches over %.1f Mb, the copy that did not "
                         "come from the sequenced parent" % (child, who, name, len(mine), bp / 1e6))
            elif "parent and grandparents" in bases:
                x = crossovers(rows, side)
                anc = sorted({r["ancestor"] for r in mine})
                L.append("- %s: from the %s (%s): %d stretches over %.1f Mb, told apart as the copy from %s; "
                         "%d crossover%s, %d discordant sites absorbed"
                         % (child, who, name, len(mine), bp / 1e6, " or ".join(anc), x, "" if x == 1 else "s",
                            sum(r["discordant"] for r in mine)))
            else:
                L.append("- %s: from the %s (%s): %d stretches over %.1f Mb; no grandparent on this side is "
                         "sequenced, so the two copies of the %s are not told apart"
                         % (child, who, name, len(mine), bp / 1e6, who))
    L += ["", "A stretch starts and ends at its first and last informative site: a crossover lies somewhere "
          "between two stretches, not at a base. Short reads give no long-range phase, so without a sequenced "
          "grandparent the map says only which parent a stretch came from."]
    return L


def run(members, vcf, par, out_dir):
    """Map every child, write the files, return the report section's lines."""
    members, links = collect(members, vcf, par)
    rows = segments(members, links)
    write_files(rows, contig_lengths(vcf), out_dir)
    return markdown(members, rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--members", required=True)
    p.add_argument("--vcf", required=True)
    p.add_argument("--par-bed")
    p.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)
    print("\n".join(run(fr.read_members(a.members), a.vcf, fr.read_par(a.par_bed), a.out_dir)))


if __name__ == "__main__":
    main()
