#!/usr/bin/env python3
"""bin/family_segments.py (step 38's chromosome map) on a synthetic
three-generation family: two grandparents, their son (the father), the
mother, and a grandchild. Genotypes are drawn from known haplotypes, so the
truth is known.

Checks:
  1. trio (father and mother sequenced, the father's parents too): on chr1
     the grandchild got the grandfather's copy up to a planted crossover at
     20 Mb and the grandmother's after it. The map has those two paternal
     segments, each boundary within 50 kb of the crossover; three planted
     genotype errors inside the first segment stay inside it (discordant 3),
     not new segments; chr2 has no crossover, and its 10 Mb gap with no
     sites splits it into two segments from the grandmother; the mother's
     side, with no grandparents sequenced, is mapped as one stretch per
     block; a son's chrX comes from his mother only, outside the PARs;
  2. duo (the mother named, not sequenced): the paternal crossover is still
     found, and the maternal side is written as the copy that did not come
     from the father, over the same stretches;
  3. the DNA Painter CSV has DNA Painter's own header, a row per segment,
     chromosome names without 'chr', and the grandparent as Ancestor;
  4. the SVG is well-formed XML with a bar per segment;
  5. family_report.py report writes the files and a "Chromosome map"
     section in family_report.md that counts the crossover; a second run
     writes the same bytes.

Run: python3 tests/test_family_segments.py
"""
import csv
import os
import random
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(REPO, "bin", "family_report.py")
SEGMENTS = os.path.join(REPO, "bin", "family_segments.py")
PAR = os.path.join(REPO, "assets", "par_grch38.bed")
sys.path.insert(0, os.path.join(REPO, "bin"))

FAILS = []
CROSSOVER = 20_000_000
STEP = 2_000


def check(desc, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + desc + ("" if ok else "  " + str(detail)))
    if not ok:
        FAILS.append(desc)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def members_tsv(path, mum_sequenced=True):
    rows = [("kid", "male", "dad", "mum", ""), ("dad", "male", "gpa", "gma", "mum"),
            ("mum", "female", "", "", "dad"), ("gpa", "male", "", "", "gma"), ("gma", "female", "", "", "gpa")]
    lines = ["sample\tsex\tfather\tmother\tpartner\tin_pedigree\tsequenced"]
    for s, sex, fa, mo, pa in rows:
        seq = "no" if s == "mum" and not mum_sequenced else "yes"
        lines.append("\t".join([s, sex, fa, mo, pa, "yes", seq]))
    return write(path, "\n".join(lines) + "\n")


def gt(a, b):
    return "%d/%d" % tuple(sorted((a, b)))


def family_vcf(path, samples):
    """Sites on chr1 (crossover in the father's meiosis at CROSSOVER), chr2
    (no crossover, a 10 Mb gap with no sites) and chrX (a son's single X from
    his mother). Returns the positions of the planted errors."""
    rnd = random.Random(38)
    recs, errors = [], []
    blocks = {"chr1": [(1_000_000, 41_000_000)], "chr2": [(1_000_000, 11_000_000), (21_000_000, 31_000_000)],
              "chrX": [(5_000_000, 15_000_000)]}
    for chrom, spans in blocks.items():
        for lo, hi in spans:
            for pos in range(lo, hi, STEP):
                h = {k: rnd.randint(0, 1) for k in ("gpa1", "gpa2", "gma1", "gma2", "mum1", "mum2")}
                planted = chrom == "chr1" and pos in (5_000_000, 9_000_000, 13_000_000)
                if planted:
                    # a site that tells everything: the father 0/1 from gpa 0/0 and gma 1/1, the mother 0/0
                    h.update(gpa1=0, gpa2=0, gma1=1, gma2=1, mum1=0, mum2=0)
                # the father carries gpa1 (from the grandfather) and gma1 (from the grandmother)
                dad = (h["gpa1"], h["gma1"])
                if chrom == "chrX":
                    g = {"gpa": "%d" % h["gpa1"], "gma": gt(h["gma1"], h["gma2"]), "dad": "%d" % h["gma1"],
                         "mum": gt(h["mum1"], h["mum2"]), "kid": "%d" % h["mum1"]}
                else:
                    from_dad = h["gpa1"] if chrom == "chr1" and pos < CROSSOVER else h["gma1"]
                    if planted:
                        from_dad = h["gma1"]   # a genotype error that looks like the other copy
                        errors.append(pos)
                    g = {"gpa": gt(h["gpa1"], h["gpa2"]), "gma": gt(h["gma1"], h["gma2"]), "dad": gt(*dad),
                         "mum": gt(h["mum1"], h["mum2"]), "kid": gt(from_dad, h["mum1"])}
                recs.append([chrom, str(pos), ".", "A", "G", "50", "PASS", ".", "GT"] + [g[s] for s in samples])
    head = ["##fileformat=VCFv4.2",
            "##contig=<ID=chr1,length=248956422>", "##contig=<ID=chr2,length=242193529>",
            "##contig=<ID=chrX,length=156040895>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(samples))]
    write(path, "\n".join(head + ["\t".join(r) for r in recs]) + "\n")
    return errors


def read_tsv(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f]


def run(script, *args):
    p = subprocess.run([sys.executable, script] + list(args), capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def near(x, target, tol=50_000):
    return abs(int(x) - target) <= tol


def test_trio(d):
    samples = ("gpa", "gma", "dad", "mum", "kid")
    vcf = os.path.join(d, "family.vcf")
    errors = family_vcf(vcf, samples)
    check("three genotype errors were planted on chr1", len(errors) == 3, errors)
    members = members_tsv(os.path.join(d, "members.tsv"))
    out = os.path.join(d, "trio")
    rc, log = run(SEGMENTS, "--members", members, "--vcf", vcf, "--par-bed", PAR, "--out-dir", out)
    check("family_segments.py exits 0", rc == 0, log)
    path = os.path.join(out, "kid_segments.tsv")
    check("kid_segments.tsv written", os.path.exists(path), log)
    if not os.path.exists(path):
        return
    rows = read_tsv(path)
    pat1 = [r for r in rows if r["side"] == "paternal" and r["chrom"] == "chr1"]
    check("chr1, paternal: two segments", len(pat1) == 2, pat1)
    if len(pat1) == 2:
        a, b = pat1
        check("chr1: first from the grandfather, then from the grandmother",
              (a["ancestor"], b["ancestor"]) == ("gpa", "gma"), (a["ancestor"], b["ancestor"]))
        check("chr1: the first segment ends within 50 kb before the crossover (%s)" % a["end"],
              near(a["end"], CROSSOVER) and int(a["end"]) < CROSSOVER)
        check("chr1: the second starts within 50 kb after it (%s)" % b["start"],
              near(b["start"], CROSSOVER) and int(b["start"]) >= CROSSOVER)
        check("chr1: the segments span the chromosome's sites",
              near(a["start"], 1_000_000) and near(b["end"], 41_000_000), (a["start"], b["end"]))
        check("chr1: the three planted errors are discordant sites of the first segment",
              a["discordant"] == "3" and b["discordant"] == "0", (a["discordant"], b["discordant"]))
        check("chr1: basis is the parent and the grandparents", a["basis"] == "parent and grandparents", a["basis"])
    pat2 = [r for r in rows if r["side"] == "paternal" and r["chrom"] == "chr2"]
    check("chr2, paternal: exactly two segments", len(pat2) == 2, pat2)
    check("chr2, paternal: two segments split by the gap, both from the grandmother",
          [(r["ancestor"], near(r["start"], s), near(r["end"], e)) for r, (s, e) in
           zip(pat2, ((1_000_000, 11_000_000), (21_000_000, 31_000_000)))]
          == [("gma", True, True), ("gma", True, True)], pat2)
    check("chrX: nothing from the father of a son",
          not [r for r in rows if r["side"] == "paternal" and r["chrom"] == "chrX"])
    matx = [r for r in rows if r["side"] == "maternal" and r["chrom"] == "chrX"]
    check("chrX: one maternal segment, from the mother", len(matx) == 1 and matx[0]["parent"] == "mum", matx)
    mat = [r for r in rows if r["side"] == "maternal" and r["chrom"] in ("chr1", "chr2")]
    check("maternal side, no grandparents sequenced: one stretch per block (3), basis parent",
          len(mat) == 3 and all(r["basis"] == "parent" and r["ancestor"] == "." for r in mat), mat)
    check("dad (a child of gpa and gma) has a map of his own",
          os.path.exists(os.path.join(out, "dad_segments.tsv")))
    check("no map for a member without a sequenced parent", not os.path.exists(os.path.join(out, "mum_segments.tsv")))

    with open(os.path.join(out, "kid_dnapainter.csv"), newline="") as f:
        r = list(csv.reader(f))
    check("DNA Painter CSV: DNA Painter's own header",
          r and r[0] == ["Chr", "Start", "End", "cM", "SNPs", "Match", "Confidence", "Group", "Side", "Notes",
                         "Colour", "Ancestor"], r[:1])
    body = r[1:]
    check("DNA Painter CSV: a row per segment", len(body) == len(rows), (len(body), len(rows)))
    check("DNA Painter CSV: chromosomes without 'chr', sides maternal or paternal",
          all(x[0] in ("1", "2", "X") and x[8] in ("maternal", "paternal") for x in body))
    first1 = next((x for x in body if x[0] == "1"), None)
    check("DNA Painter CSV: the first chr1 segment is in the grandparent's group, marked a known ancestor",
          first1 is not None and first1[7] == "gpa" and first1[11] == "1", first1)
    check("DNA Painter CSV: the other grandparent has a group of its own",
          any(x[7] == "gma" and x[11] == "1" for x in body), body)
    check("DNA Painter CSV: Ancestor is only 1 or empty", all(x[11] in ("1", "") for x in body), body)
    try:
        svg = ET.parse(os.path.join(out, "kid_segments.svg")).getroot()
        bars = [e for e in svg.iter() if e.tag.endswith("rect") and e.get("class") == "seg"]
        check("SVG: well-formed, a bar per segment", len(bars) == len(rows), (len(bars), len(rows)))
    except (ET.ParseError, OSError) as e:
        check("SVG: well-formed", False, e)
    return vcf, members


def test_duo(d, vcf):
    members = members_tsv(os.path.join(d, "members_duo.tsv"), mum_sequenced=False)
    out = os.path.join(d, "duo")
    rc, log = run(SEGMENTS, "--members", members, "--vcf", vcf, "--par-bed", PAR, "--out-dir", out)
    check("duo: exits 0", rc == 0, log)
    path = os.path.join(out, "kid_segments.tsv")
    if not os.path.exists(path):
        check("duo: kid_segments.tsv written", False, log)
        return
    rows = read_tsv(path)
    pat1 = [r for r in rows if r["side"] == "paternal" and r["chrom"] == "chr1"]
    check("duo, chr1: the crossover is found (two segments, gpa then gma)",
          [r["ancestor"] for r in pat1] == ["gpa", "gma"], pat1)
    if len(pat1) == 2:
        check("duo, chr1: boundary within 50 kb of the crossover",
              near(pat1[0]["end"], CROSSOVER) and near(pat1[1]["start"], CROSSOVER), pat1)
    mat = [r for r in rows if r["side"] == "maternal"]
    check("duo: the maternal side is the other copy, on chr1 and chr2 only",
          mat and all(r["basis"] == "other copy" and r["parent"] == "mum" for r in mat)
          and {r["chrom"] for r in mat} == {"chr1", "chr2"}, mat)
    want = [("chr1", 1_000_000, 41_000_000), ("chr2", 1_000_000, 11_000_000), ("chr2", 21_000_000, 31_000_000)]
    check("duo: the other copy covers chr1 and the two chr2 stretches, ends within 50 kb",
          len(mat) == 3 and all(r["chrom"] == c and near(r["start"], s) and near(r["end"], e)
                                for r, (c, s, e) in zip(mat, want)), mat)


def test_no_par(d, vcf, members):
    out = os.path.join(d, "nopar")
    rc, log = run(SEGMENTS, "--members", members, "--vcf", vcf, "--out-dir", out)
    check("no --par-bed: exits 0", rc == 0, log)
    path = os.path.join(out, "kid_segments.tsv")
    rows = read_tsv(path) if os.path.exists(path) else []
    check("no --par-bed: chrX is left out, the autosomes are still mapped",
          rows and not [r for r in rows if r["chrom"] == "chrX"], rows)


def test_report(d, vcf, members):
    rel = write(os.path.join(d, "relationships.tsv"),
                "sample_a\tsample_b\tdeclared\texpected_relatedness\trelatedness\tibs0_rate\tsites\tverdict\treason\n")
    outs = []
    for k in (1, 2):
        out = os.path.join(d, "report%d" % k)
        rc, log = run(REPORT, "report", "--members", members, "--relationships", rel, "--vcf", vcf,
                      "--par-bed", PAR, "--out-dir", out)
        check("report run %d exits 0" % k, rc == 0, log[-2000:])
        outs.append(out)
    md = open(os.path.join(outs[0], "family_report.md")).read() if os.path.exists(
        os.path.join(outs[0], "family_report.md")) else ""
    check("family_report.md has a Chromosome map section", "## Chromosome map" in md)
    check("the section counts kid's one paternal crossover",
          "kid: from the father (dad)" in md and "1 crossover" in md, md[md.find("## Chromosome map"):][:1500])
    names = ("kid_segments.tsv", "kid_dnapainter.csv", "kid_segments.svg", "dad_segments.tsv")
    for n in names:
        a, b = (os.path.join(o, n) for o in outs)
        same = os.path.exists(a) and os.path.exists(b) and open(a, "rb").read() == open(b, "rb").read()
        check("report writes %s, the same bytes on a second run" % n, same)


def main():
    with tempfile.TemporaryDirectory() as d:
        res = test_trio(d)
        if res:
            vcf, members = res
            test_duo(d, vcf)
            test_no_par(d, vcf, members)
            test_report(d, vcf, members)
        else:
            check("the trio map ran", False)
    if FAILS:
        print("%d check(s) failed" % len(FAILS))
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
