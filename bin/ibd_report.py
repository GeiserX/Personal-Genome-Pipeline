#!/usr/bin/env python3
"""ibd_report.py: the DNA segments two sequenced relatives share by descent,
for step 38 (scripts/38-family.sh) and the Nextflow FAMILY workflow
(modules/local/ibd). Standard library only, so it runs in the plain python
image.

  ibd_report.py --pedigree P.tsv --chroms chroms.txt --maps MAPS --phased PHASED --ibd IBD --out-dir D

Reads, for each chromosome of chroms.txt (one per line): Beagle's phased
family VCF PHASED/<chrom>.vcf.gz, hap-ibd's segments IBD/<chrom>.ibd.gz and
the PLINK genetic map MAPS/<chrom>.map. hap-ibd reports a segment per pair of
haplotypes; here the segments of one pair of people are merged into the
stretches they share (one copy or both), as consumer DNA sites show them.
A phasing error breaks hap-ibd's segments, and a piece shorter than its
2 cM seed between two breaks is lost. Two people who share a copy are never
opposite homozygotes (one 0/0, the other 1/1) bar a genotype error, so two
stretches of a pair are joined when the markers between them hold at most
one opposite homozygote, and a pair's first and last stretch reach the
chromosome's first and last marker on the same test. The Browning lab's
merge-ibd-segments joins segments by the same test in gaps under 0.6 cM;
here a gap has no length limit, since two unrelated people are opposite
homozygotes at several percent of the panel's common SNVs.

Writes:
  D/ibd_segments.tsv  one row per shared stretch: sample_a, sample_b, chrom,
                      start, end, cm, ibd2_cm (the part where both copies
                      are shared)
  D/ibd_pairs.tsv     one row per pair of members: the declared relationship,
                      the stretches, total and longest cM, the share of the
                      analysed genetic map, the relatedness it gives and the
                      degree it points to
  D/ibd_report.md     the same in words

The analysed map is, per chromosome, the cM between the first and last
phased marker. Relatedness is twice the kinship coefficient,
(IBD1 / 4 + IBD2 / 2) / analysed cM * 2, binned as KING bins it (the same
bins as family_report.py, plus a fourth degree).
"""
import argparse
import bisect
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_report as fr  # noqa: E402

PAIR_COLUMNS = ("sample_a", "sample_b", "declared", "expected_relatedness", "segments", "shared_cm", "ibd2_cm",
                "longest_cm", "analysed_cm", "shared_fraction", "relatedness", "degree")
SEGMENT_COLUMNS = ("sample_a", "sample_b", "chrom", "start", "end", "cm", "ibd2_cm")
# hap-ibd's default min-output: shorter segments are not reported.
MIN_OUTPUT_CM = 2.0
# Opposite homozygotes allowed between two stretches that are joined.
MAX_DISCORD = 1
# Share of the first-degree pair's genome in both copies above which the
# pair reads as full siblings (about 25%) rather than parent and child (0%).
SIBLING_MIN_IBD2 = 0.10
DEGREES = ((0.707, "same person or identical twins"), (0.354, "first degree"), (0.177, "second degree"),
           (0.0884, "third degree"), (0.0442, "fourth degree"))


def die(msg, code=2):
    sys.stderr.write("ERROR: %s\n" % msg)
    sys.exit(code)


# --- inputs -----------------------------------------------------------------

class Map:
    """A PLINK genetic map (chrom, id, cM, bp); cM at any position by linear
    interpolation, held flat beyond the first and last row."""

    def __init__(self, path):
        rows = []
        with open(path) as f:
            for n, line in enumerate(f, 1):
                c = line.split()
                if not c:
                    continue
                if len(c) < 4:
                    die("%s line %d: a PLINK map row has 4 columns (chrom, id, cM, bp)" % (path, n))
                try:
                    rows.append((int(c[3]), float(c[2])))
                except ValueError:
                    die("%s line %d: cM and bp must be numbers" % (path, n))
        if len(rows) < 2:
            die("%s: a genetic map needs two rows or more" % path)
        rows.sort()
        self.bp = [r[0] for r in rows]
        self.cm = [r[1] for r in rows]

    def at(self, pos):
        i = bisect.bisect_left(self.bp, pos)
        if i == 0:
            return self.cm[0]
        if i == len(self.bp):
            return self.cm[-1]
        b0, b1 = self.bp[i - 1], self.bp[i]
        if b1 == b0:
            return self.cm[i]
        return self.cm[i - 1] + (self.cm[i] - self.cm[i - 1]) * (pos - b0) / (b1 - b0)


def dosage(gt):
    a = gt.split(":", 1)[0].replace("|", "/").split("/")
    return sum(int(x) for x in a) if len(a) == 2 and all(x in ("0", "1") for x in a) else -1


def read_phased(path):
    """Sample names, marker positions and each sample's ALT count per marker
    (-1 when not 0/1 alleles) of a phased VCF."""
    samples, pos, dos = None, [], None
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                dos = [[] for _ in samples]
            elif not line.startswith("#"):
                if samples is None:
                    die("%s has no #CHROM line before its records" % path)
                c = line.rstrip("\r\n").split("\t")
                pos.append(int(c[1]))
                for i, g in enumerate(c[9:9 + len(samples)]):
                    dos[i].append(dosage(g))
    if samples is None:
        die("%s has no #CHROM line" % path)
    return samples, pos, dict(zip(samples, dos))


def opposite_homozygotes(pos, dx, dy, a, b):
    """Markers strictly between positions a and b where one is 0/0 and the other 1/1."""
    i, j = bisect.bisect_right(pos, a), bisect.bisect_left(pos, b)
    return sum(1 for k in range(i, j) if {dx[k], dy[k]} == {0, 2})


def stitch(iv, pos, dx, dy):
    """Join stretches with at most MAX_DISCORD opposite homozygotes between
    them, and reach the first and last marker on the same test."""
    if not iv:
        return iv
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if opposite_homozygotes(pos, dx, dy, out[-1][1], s) <= MAX_DISCORD:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    if opposite_homozygotes(pos, dx, dy, pos[0] - 1, out[0][0]) <= MAX_DISCORD:
        out[0][0] = min(out[0][0], pos[0])
    if opposite_homozygotes(pos, dx, dy, out[-1][1], pos[-1] + 1) <= MAX_DISCORD:
        out[-1][1] = max(out[-1][1], pos[-1])
    return [tuple(x) for x in out]


def read_ibd(path, chrom):
    """hap-ibd's .ibd.gz: sample1, hap1, sample2, hap2, chrom, start, end, cM."""
    rows = []
    with gzip.open(path, "rt") as f:
        for n, line in enumerate(f, 1):
            c = line.rstrip("\r\n").split("\t")
            if c == [""]:
                continue
            if len(c) < 8:
                die("%s line %d: hap-ibd rows have 8 columns, this one %d" % (path, n, len(c)))
            if c[4] != chrom:
                die("%s line %d: chromosome %s in the file of %s" % (path, n, c[4], chrom))
            try:
                a, ha, b, hb, s, e = c[0], int(c[1]), c[2], int(c[3]), int(c[5]), int(c[6])
            except ValueError:
                die("%s line %d: haplotypes and positions must be integers" % (path, n))
            if ha not in (1, 2) or hb not in (1, 2) or e < s:
                die("%s line %d: a haplotype is 1 or 2 and the end is not before the start" % (path, n))
            if a == b:
                continue
            if b < a:
                a, ha, b, hb = b, hb, a, ha
            rows.append((a, b, ha, hb, s, e))
    return rows


# --- intervals --------------------------------------------------------------

def union(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def intersect(x, y):
    out, i, j = [], 0, 0
    while i < len(x) and j < len(y):
        s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
        if s < e:
            out.append((s, e))
        if x[i][1] < y[j][1]:
            i += 1
        else:
            j += 1
    return out


def ibd2(segs):
    """Where both copies are shared: a's haplotype 1 with one of b's and a's
    haplotype 2 with the other."""
    u = {(ha, hb): union([(s, e) for ha_, hb_, s, e in segs if (ha_, hb_) == (ha, hb)])
         for ha in (1, 2) for hb in (1, 2)}
    return union(intersect(u[1, 1], u[2, 2]) + intersect(u[1, 2], u[2, 1]))


# --- the estimate -----------------------------------------------------------

def degree(rel, ibd2_frac, n):
    if n == 0:
        return "no shared segment"
    for lo, name in DEGREES:
        if rel >= lo:
            if name == "first degree":
                return name + (" (full siblings)" if ibd2_frac >= SIBLING_MIN_IBD2 else " (parent-child)")
            return name
    return "distant or unrelated (short segments only)"


def declared_table(path, samples):
    rows, errors = fr.read_pedigree(path)
    members, more = fr.build_members(rows, {}, path)
    errors += more
    if errors:
        die("; ".join(errors))
    phi = fr.kinship_table(members)

    def declared(a, b):
        if a not in members or b not in members:
            return "not in the pedigree", ""
        e = 2 * phi(a, b)
        return fr.declared_relationship(a, b, members, e), "%.3f" % e
    return declared


def f2(x):
    return "%.2f" % x


def build(a):
    chroms = [c.strip() for c in open(a.chroms) if c.strip()]
    samples, glen, by_pair, maps, geno = None, {}, {}, {}, {}
    for c in chroms:
        maps[c] = Map(os.path.join(a.maps, c + ".map"))
        s, pos, dos = read_phased(os.path.join(a.phased, c + ".vcf.gz"))
        if samples is None:
            samples = s
        elif s != samples:
            die("%s.vcf.gz has the samples %s, %s.vcf.gz %s" % (c, ",".join(s), chroms[0], ",".join(samples)))
        if not pos:
            die("%s.vcf.gz has no marker" % c)
        geno[c] = (pos, dos)
        glen[c] = maps[c].at(pos[-1]) - maps[c].at(pos[0])
        for x, y, ha, hb, s0, e0 in read_ibd(os.path.join(a.ibd, c + ".ibd.gz"), c):
            if x not in s or y not in s:
                die("%s.ibd.gz names %s and %s; the phased VCF holds %s" % (c, x, y, ",".join(s)))
            by_pair.setdefault((x, y), {}).setdefault(c, []).append((ha, hb, s0, e0))
    samples = sorted(samples or [])
    total = sum(glen.values())
    declared = declared_table(a.pedigree, samples)

    pairs, segments = [], []
    for i, x in enumerate(samples):
        for y in samples[i + 1:]:
            shared = both = longest = 0.0
            n = 0
            for c in chroms:
                segs = by_pair.get((x, y), {}).get(c, [])
                if not segs:
                    continue
                m = maps[c]
                pos, dos = geno[c]
                two = ibd2(segs)
                for s0, e0 in stitch(union([(s, e) for _, _, s, e in segs]), pos, dos[x], dos[y]):
                    cm = m.at(e0) - m.at(s0)
                    cm2 = sum(m.at(e) - m.at(s) for s, e in intersect([(s0, e0)], two))
                    segments.append((x, y, c, s0, e0, cm, cm2))
                    shared += cm
                    both += cm2
                    longest = max(longest, cm)
                    n += 1
            frac = shared / total if total > 0 else 0.0
            frac2 = both / total if total > 0 else 0.0
            rel = 2 * ((frac - frac2) / 4 + frac2 / 2)
            d, e = declared(x, y)
            pairs.append({"sample_a": x, "sample_b": y, "declared": d, "expected_relatedness": e, "segments": n,
                          "shared_cm": shared, "ibd2_cm": both, "longest_cm": longest, "analysed_cm": total,
                          "shared_fraction": frac, "relatedness": rel, "degree": degree(rel, frac2, n)})
    return chroms, glen, total, pairs, segments


def write(a, chroms, glen, total, pairs, segments):
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "ibd_segments.tsv"), "w") as f:
        f.write("\t".join(SEGMENT_COLUMNS) + "\n")
        for x, y, c, s, e, cm, cm2 in segments:
            f.write("\t".join([x, y, c, str(s), str(e), f2(cm), f2(cm2)]) + "\n")
    with open(os.path.join(a.out_dir, "ibd_pairs.tsv"), "w") as f:
        f.write("\t".join(PAIR_COLUMNS) + "\n")
        for p in pairs:
            row = dict(p)
            for k in ("shared_cm", "ibd2_cm", "longest_cm", "analysed_cm"):
                row[k] = f2(p[k])
            row["shared_fraction"] = "%.4f" % p["shared_fraction"]
            row["relatedness"] = "%.3f" % p["relatedness"]
            f.write("\t".join(str(row[k]) for k in PAIR_COLUMNS) + "\n")

    L = ["# Shared DNA segments", ""]
    if not chroms:
        L.append("No chromosome was phased, so no shared segment was looked for.")
    else:
        L.append("hap-ibd looked for the stretches of DNA each pair of members shares by descent, on the joint "
                 "call phased by Beagle against the reference panel. Analysed: %d chromosome%s (%s), %.1f cM of "
                 "genetic map. Only stretches of %.0f cM or more count: shorter ones cannot be told apart from "
                 "chance sharing." % (len(chroms), "" if len(chroms) == 1 else "s", ", ".join(chroms), total,
                                      MIN_OUTPUT_CM))
        if len(chroms) < 22:
            L += ["", "Fewer than the 22 autosomes were analysed, so the degree below is a rough estimate."]
        L += ["", "| Pair | Declared | Shared | Of the analysed map | Longest | Stretches | Both copies | Estimate |",
              "|---|---|---|---|---|---|---|---|"]
        for p in pairs:
            L.append("| %s and %s | %s | %.1f cM | %.1f%% | %.1f cM | %d | %.1f cM | %s |"
                     % (p["sample_a"], p["sample_b"], p["declared"], p["shared_cm"], 100 * p["shared_fraction"],
                        p["longest_cm"], p["segments"], p["ibd2_cm"], p["degree"]))
        L += ["", "Every stretch is in ibd_segments.tsv (chromosome, start, end, cM, and the cM where both "
                  "copies are shared). Two stretches a phasing error split are joined when the markers between "
                  "them hold at most one opposite homozygote (one 0/0, the other 1/1), which two people who share "
                  "a copy never are bar a genotype error. A parent and child share one copy along the whole "
                  "genome; full siblings share both copies over about a quarter of it."]
    with open(os.path.join(a.out_dir, "ibd_report.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--pedigree", required=True)
    p.add_argument("--chroms", required=True)
    p.add_argument("--maps", required=True)
    p.add_argument("--phased", required=True)
    p.add_argument("--ibd", required=True)
    p.add_argument("--out-dir", required=True)
    a = p.parse_args(argv)
    write(a, *build(a))
    return 0


if __name__ == "__main__":
    sys.exit(main())
