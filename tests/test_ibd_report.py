#!/usr/bin/env python3
"""bin/ibd_report.py: the DNA segments two sequenced relatives share (step 38).

Each case writes synthetic Beagle VCFs, hap-ibd segment files, PLINK maps and
a pedigree, runs the script, and checks what it writes:

  1. the output schema: the exact header of ibd_pairs.tsv and ibd_segments.tsv,
     one pairs row for every pair of members (a pair with no segment too),
     one segments row per shared stretch
  2. a parent and child: hap-ibd's two segments (the parent's haplotype
     changes at a crossover) are one stretch over the whole analysed map,
     shared fraction 1, nothing in both copies, first degree (parent-child)
  3. full siblings: both copies shared over a quarter of the map is IBD2,
     relatedness from (IBD1 / 4 + IBD2 / 2), first degree (full siblings)
  4. a short segment only, written in hap-ibd's other sample order, with
     opposite homozygotes on both sides: distant;
     a pair with no segment: no shared segment; a member without a pedigree
     row: not in the pedigree
  5. two chromosomes: the analysed map is the sum of first to last marker of
     each, read on the genetic map
  6. the report says it in words, and that one chromosome is a rough estimate
  7. damaged input stops the script (exit 2): a hap-ibd row with 7 columns,
     a row on another chromosome, a sample the phased VCF does not hold, a
     map row without its bp column
  8. no chromosome phased: headers only, and the report says so
  9. stretches broken by phasing: two joined across a gap with no opposite
     homozygote, not across one with two; a stretch reaches the chromosome's
     first and last marker across one opposite homozygote, not across three

Run: python3 tests/test_ibd_report.py
"""
import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "ibd_report.py")
FAILS = []
MB = 1000000

PAIRS_HEAD = ["sample_a", "sample_b", "declared", "expected_relatedness", "segments", "shared_cm", "ibd2_cm",
              "longest_cm", "analysed_cm", "shared_fraction", "relatedness", "degree"]
SEGS_HEAD = ["sample_a", "sample_b", "chrom", "start", "end", "cm", "ibd2_cm"]
PEDIGREE = ("sample\tsex\tfather\tmother\tpartner\n"
            "A\tmale\tB\tX\t.\nC\tfemale\tB\tX\t.\nB\tmale\t.\t.\tX\n")


def check(name, got, want):
    if got == want:
        print(f"PASS {name}")
    else:
        print(f"FAIL {name}: got {got!r}, want {want!r}")
        FAILS.append(name)


GT = {0: "0|0", 1: "0|1", 2: "1|1"}


def write_chrom(d, chrom, samples, first, last, ibd_rows, map_end, dos=lambda s, mb: 1):
    """1 cM per Mb from 0 to map_end Mb; markers every Mb from first to last Mb;
    dos(sample, Mb) is the sample's ALT count there."""
    for sub in ("maps", "phased", "ibd"):
        os.makedirs(os.path.join(d, sub), exist_ok=True)
    with open(os.path.join(d, "maps", chrom + ".map"), "w") as f:
        for i in range(0, map_end + 1):
            f.write(f"{chrom}\t.\t{float(i):.4f}\t{max(1, i * MB)}\n")
    with gzip.open(os.path.join(d, "phased", chrom + ".vcf.gz"), "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        for i in range(first, last + 1):
            f.write(f"{chrom}\t{i * MB}\t.\tA\tG\t.\tPASS\t.\tGT\t" + "\t".join(GT[dos(s, i)] for s in samples) + "\n")
    with gzip.open(os.path.join(d, "ibd", chrom + ".ibd.gz"), "wt") as f:
        for r in ibd_rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def run(d, chroms, pedigree=PEDIGREE):
    with open(os.path.join(d, "chroms.txt"), "w") as f:
        f.write("".join(c + "\n" for c in chroms))
    with open(os.path.join(d, "pedigree.tsv"), "w") as f:
        f.write(pedigree)
    out = os.path.join(d, "out")
    r = subprocess.run([sys.executable, SCRIPT, "--pedigree", os.path.join(d, "pedigree.tsv"),
                        "--chroms", os.path.join(d, "chroms.txt"), "--maps", os.path.join(d, "maps"),
                        "--phased", os.path.join(d, "phased"), "--ibd", os.path.join(d, "ibd"), "--out-dir", out],
                       capture_output=True, text=True)
    return r, out


def table(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return head, [dict(zip(head, l.rstrip("\n").split("\t"))) for l in f if l.strip()]


SAMPLES = ["A", "B", "C", "D"]
# D is 1/1 where A, B and C are 0/0, except over its short segment with B
DOS1 = lambda s, mb: (1 if 10 <= mb <= 13 else 2) if s == "D" else 0  # noqa: E731
CHR1 = [("A", 1, "B", 1, "chr1", 1 * MB, 50 * MB, 49.0),      # parent-child, then the
        ("A", 1, "B", 2, "chr1", 50 * MB, 99 * MB, 49.0),     # parent's other haplotype
        ("A", 1, "C", 1, "chr1", 1 * MB, 99 * MB, 98.0),      # siblings: one copy everywhere,
        ("A", 2, "C", 2, "chr1", 20 * MB, 44500000, 24.5),    # both over a quarter
        ("D", 2, "B", 1, "chr1", 10 * MB, 13 * MB, 3.0)]      # short, in the other order


with tempfile.TemporaryDirectory() as d:
    write_chrom(d, "chr1", SAMPLES, 1, 99, CHR1, 100, DOS1)
    r, out = run(d, ["chr1"])
    check("1 chromosome: exit 0", r.returncode, 0)
    if r.returncode:
        print(r.stderr)
    head, pairs = table(os.path.join(out, "ibd_pairs.tsv"))
    check("1 ibd_pairs.tsv header", head, PAIRS_HEAD)
    check("1 one row for every pair of the four members", [(p["sample_a"], p["sample_b"]) for p in pairs],
          [("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D")])
    shead, segs = table(os.path.join(out, "ibd_segments.tsv"))
    check("1 ibd_segments.tsv header", shead, SEGS_HEAD)
    check("1 one segments row per shared stretch", [(s["sample_a"], s["sample_b"]) for s in segs],
          [("A", "B"), ("A", "C"), ("B", "D")])
    P = {(p["sample_a"], p["sample_b"]): p for p in pairs}

    ab = P["A", "B"]
    check("2 parent-child: declared", ab["declared"], "parent-child")
    check("2 parent-child: expected relatedness from the pedigree", ab["expected_relatedness"], "0.500")
    check("2 parent-child: the two hap-ibd segments are one stretch", ab["segments"], "1")
    check("2 parent-child: shared cM over the whole map", (ab["shared_cm"], ab["analysed_cm"]), ("98.00", "98.00"))
    check("2 parent-child: longest stretch", ab["longest_cm"], "98.00")
    check("2 parent-child: shared fraction and nothing in both copies", (ab["shared_fraction"], ab["ibd2_cm"]),
          ("1.0000", "0.00"))
    check("2 parent-child: relatedness and degree", (ab["relatedness"], ab["degree"]),
          ("0.500", "first degree (parent-child)"))
    check("2 parent-child: the stretch row", [segs[0][k] for k in SEGS_HEAD],
          ["A", "B", "chr1", str(1 * MB), str(99 * MB), "98.00", "0.00"])

    ac = P["A", "C"]
    check("3 siblings: declared", ac["declared"], "full siblings")
    check("3 siblings: both copies over a quarter", ac["ibd2_cm"], "24.50")
    check("3 siblings: relatedness (IBD1 / 4 + IBD2 / 2) * 2", ac["relatedness"], "0.625")
    check("3 siblings: degree", ac["degree"], "first degree (full siblings)")

    bd = P["B", "D"]
    check("4 short segment in the other order: counted for B and D", (bd["segments"], bd["shared_cm"]), ("1", "3.00"))
    check("4 short segment: distant", bd["degree"], "distant or unrelated (short segments only)")
    check("4 short segment: D has no pedigree row", bd["declared"], "not in the pedigree")
    check("4 no segment", (P["C", "D"]["segments"], P["C", "D"]["degree"]), ("0", "no shared segment"))
    md = open(os.path.join(out, "ibd_report.md")).read()
    check("6 the report has the pair's line", "| A and B | parent-child | 98.0 cM | 100.0% | 98.0 cM | 1 | 0.0 cM | "
          "first degree (parent-child) |" in md, True)
    check("6 the report says one chromosome is a rough estimate", "rough estimate" in md, True)
    check("6 the report names what was analysed", "Analysed: 1 chromosome (chr1), 98.0 cM" in md, True)

    write_chrom(d, "chr2", SAMPLES, 2, 49, [("B", 2, "A", 2, "chr2", 2 * MB, 49 * MB, 47.0)], 50, DOS1)
    r, out = run(d, ["chr1", "chr2"])
    check("5 two chromosomes: exit 0", r.returncode, 0)
    _, pairs = table(os.path.join(out, "ibd_pairs.tsv"))
    ab = {(p["sample_a"], p["sample_b"]): p for p in pairs}["A", "B"]
    check("5 two chromosomes: the analysed map is the sum (98 + 47 cM)", ab["analysed_cm"], "145.00")
    check("5 two chromosomes: A and B, two stretches, all of it shared",
          (ab["segments"], ab["shared_cm"], ab["longest_cm"], ab["shared_fraction"]), ("2", "145.00", "98.00", "1.0000"))

# 9. P and Q: 1-30 and 35-60 Mb joined (nothing between), 70-99 Mb not (Q is
# 1/1 at 62 and 66 Mb); P and R: 5-50 Mb, R 1/1 at 1-3 and 51-99 Mb, so no
# reach to the ends; P and S: 3-97 Mb reaches 1 and 99 Mb, S 1/1 at 99 Mb only.
def dos9(s, mb):
    if s == "Q":
        return 2 if mb in (62, 66) else 0
    if s == "R":
        return 2 if mb <= 3 or mb >= 51 else 0
    if s == "S":
        return 2 if mb == 99 else 0
    return 0


with tempfile.TemporaryDirectory() as d:
    write_chrom(d, "chr3", ["P", "Q", "R", "S"], 1, 99,
                [("P", 1, "Q", 1, "chr3", 1 * MB, 30 * MB, 29.0), ("P", 1, "Q", 2, "chr3", 35 * MB, 60 * MB, 25.0),
                 ("P", 1, "Q", 1, "chr3", 70 * MB, 99 * MB, 29.0), ("P", 2, "R", 1, "chr3", 5 * MB, 50 * MB, 45.0),
                 ("P", 2, "S", 2, "chr3", 3 * MB, 97 * MB, 94.0)], 100, dos9)
    r, out = run(d, ["chr3"], "sample\tsex\tfather\tmother\tpartner\nP\tmale\t.\t.\t.\n")
    check("9 stitching: exit 0", r.returncode, 0)
    _, segs = table(os.path.join(out, "ibd_segments.tsv"))
    check("9 stitching: the stretches", [(x["sample_a"], x["sample_b"], x["start"], x["end"]) for x in segs],
          [("P", "Q", str(1 * MB), str(60 * MB)), ("P", "Q", str(70 * MB), str(99 * MB)),
           ("P", "R", str(5 * MB), str(50 * MB)), ("P", "S", str(1 * MB), str(99 * MB))])
    _, pairs = table(os.path.join(out, "ibd_pairs.tsv"))
    pq = {(p["sample_a"], p["sample_b"]): p for p in pairs}["P", "Q"]
    check("9 stitching: P and Q, two stretches, 88 cM, longest 59", (pq["segments"], pq["shared_cm"], pq["longest_cm"]),
          ("2", "88.00", "59.00"))

for name, rows, mapline in (
        ("a hap-ibd row with 7 columns", [("A", 1, "B", 1, "chr1", MB, 5 * MB)], None),
        ("a row on another chromosome", [("A", 1, "B", 1, "chr9", MB, 5 * MB, 4.0)], None),
        ("a sample the phased VCF does not hold", [("A", 1, "Z", 1, "chr1", MB, 5 * MB, 4.0)], None),
        ("a map row without its bp column", CHR1, "chr1\t.\t3.0\n")):
    with tempfile.TemporaryDirectory() as d:
        write_chrom(d, "chr1", SAMPLES, 1, 99, rows, 100)
        if mapline:
            with open(os.path.join(d, "maps", "chr1.map"), "a") as f:
                f.write(mapline)
        r, out = run(d, ["chr1"])
        check(f"7 {name}: exit 2", r.returncode, 2)
        check(f"7 {name}: says why", r.stderr.startswith("ERROR: "), True)

with tempfile.TemporaryDirectory() as d:
    r, out = run(d, [])
    check("8 no chromosome: exit 0", r.returncode, 0)
    check("8 no chromosome: headers only", [open(os.path.join(out, f)).read() for f in ("ibd_pairs.tsv", "ibd_segments.tsv")],
          ["\t".join(PAIRS_HEAD) + "\n", "\t".join(SEGS_HEAD) + "\n"])
    check("8 no chromosome: the report says so", "No chromosome was phased" in open(os.path.join(out, "ibd_report.md")).read(), True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
