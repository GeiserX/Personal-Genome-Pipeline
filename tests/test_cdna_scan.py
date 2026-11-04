#!/usr/bin/env python3
"""bin/cdna_scan.py: the genome-wide count of reads across exon junctions.

A made-up 20 kb contig, chrT, with one gene of three exons (E1, E2, E3).

Part 1, plain Python (reads as small objects, no BAM):
  spliced   reads across E1-E2, E2-E3 and E1-E3 (an isoform without E2), as
            a DNA aligner writes cDNA (clipped at the exon end, the rest as
            an SA alignment at the next exon's start), plus one read with an
            N gap from E1 to E2: one locus, 3 junctions, 9 reads
  genomic   the same gene's genomic reads, running on into the introns:
            nothing listed (the red control of the spliced case)
  one       reads across E1-E2 only: one junction, not listed
  thin      one read per junction: no junction reaches MIN_READS
  two       one read across two junctions counts on each, and once in the
            locus total
  and duplicate, secondary and repeated alignments that do not count, a
  supplementary record of a counted read that does not count it twice, a
  read whose primary lies on another contig that counts from its
  supplementary records alone, SA pieces on the other strand or out of read
  order that do not pair, the MIN_INTRON bound, and a control: with MIN_CLIP
  above the read length the spliced fixture lists nothing, so the listed
  check depends on the clip test.

Part 2, the same fixtures as BAM files written with pysam, through main():
the table, its comment line, the log, --skip, and a missing exon file.
Part 2 needs pysam; without it on the host the test runs part 2 in
PYPGX_IMAGE (versions.env) through docker, as step 06 runs the script.

Run: python3 tests/test_cdna_scan.py
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import cdna_scan as cs  # noqa: E402
import clip_check as cc  # noqa: E402

FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


class Read:
    """The pysam.AlignedSegment attributes cdna_scan reads."""

    def __init__(self, name, start, cigar, sa="", read1=True, flag=0):
        self.query_name, self.reference_start = name, start
        self.cigartuples = cc.parse_cigar(cigar)
        self.is_read1 = read1
        self.is_unmapped = self.is_qcfail = False
        self.is_secondary, self.is_supplementary, self.is_duplicate = bool(flag & 256), bool(flag & 2048), bool(flag & 1024)
        self.is_reverse = bool(flag & 16)
        self._sa = sa

    def has_tag(self, t):
        return t == "SA" and bool(self._sa)

    def get_tag(self, t):
        return self._sa


E1, E2, E3 = (1000, 1200), (3000, 3150), (5000, 5200)
LENGTH = 20000


def exons_file(tmp):
    p = os.path.join(tmp, "exons.bed")
    with open(p, "w") as f:
        for a, b in (E1, E2, E3):
            f.write(f"chrT\t{a}\t{b}\n")
    return p


def across(name, left, right, aligned):
    """A 100 bp read across the junction from exon left to exon right: aligned
    bases on left's end, clipped there, the rest an SA alignment at right's start."""
    clip = 100 - aligned
    return (name, left[1] - aligned, f"{aligned}M{clip}S", f"chrT,{right[0] + 1},+,{aligned}S{clip}M,60,0;", 0)


def fixture(kind):
    """(name, start, cigar, sa, flag) per read."""
    if kind == "genomic":
        return [(f"g{i}", e[1] - 50 + 7 * i, "100M", "", 0) for e in (E1, E2, E3) for i in range(4)] + \
               [(f"s{i}", e[0] - 50 + 7 * i, "100M", "", 0) for e in (E2, E3) for i in range(4)]
    if kind == "one":
        return [across(f"a{i}", E1, E2, 60 + 5 * i) for i in range(4)]
    if kind == "thin":
        return [across("a0", E1, E2, 60), across("b0", E2, E3, 60), across("c0", E1, E3, 60)]
    return ([across(f"a{i}", E1, E2, 60 + 5 * i) for i in range(3)]
            + [("n0", E1[1] - 50, f"50M{E2[0] - E1[1]}N50M", "", 0)]
            + [across(f"b{i}", E2, E3, 60 + 5 * i) for i in range(3)]
            + [across(f"c{i}", E1, E3, 60 + 5 * i) for i in range(2)])


def objects(kind):
    return [Read(n, s, c, sa, flag=f) for n, s, c, sa, f in fixture(kind)]


SPLICED_ROW = ["chrT", str(E1[1]), str(E3[0] + 1), "3", "9",
               f"{E1[1]}-{E2[0] + 1}:4,{E1[1]}-{E3[0] + 1}:2,{E2[1]}-{E3[0] + 1}:3"]


def part1():
    with tempfile.TemporaryDirectory() as tmp:
        ex = cs.Edges(exons_file(tmp))
    expect("merged exon regions", ex.regions["chrT"], [list(E1), list(E2), list(E3)])
    expect(f"an exon end {cc.EDGE_SLACK} bp off is found", cs.Edges.nearest(ex.ends["chrT"], E1[1] + cc.EDGE_SLACK), E1[1])
    expect(f"{cc.EDGE_SLACK + 1} bp off is not", cs.Edges.nearest(ex.ends["chrT"], E1[1] + cc.EDGE_SLACK + 1), None)

    def listed(kind, extra=()):
        return cs.loci(cs.scan([("chrT", objects(kind) + list(extra))], ex))

    expect("spliced: one locus with 3 junctions and 9 reads", listed("spliced"), [SPLICED_ROW])
    expect("genomic only: nothing listed (red control)", listed("genomic"), [])
    expect("one junction alone: not listed", listed("one"), [])
    expect("one read per junction: not listed", listed("thin"), [])
    a = across("a0", E1, E2, 60)
    dup = [Read(a[0], a[1], a[2], a[3]),                      # the same read again (another region fetched it)
           Read("d", a[1], a[2], a[3], flag=1024),            # a duplicate
           Read("s", a[1], a[2], a[3], flag=256),             # a secondary alignment
           Read(a[0], E2[0], "60S40M", f"chrT,{a[1] + 1},+,60M40S,60,0;", flag=2048)]  # its supplementary record
    expect("a repeated read, a duplicate, a secondary alignment and the supplementary record of a counted read "
           "do not count", listed("spliced", dup), [SPLICED_ROW])
    # A 150 bp read: bases 0-50 on a vector contig (the primary), 50-100 at
    # E1's end, 100-150 at E2's start. Only its two supplementary records sit
    # in an exon region of chrT.
    vec = [Read("v", E1[1] - 50, "50S50M50S", f"chrV,1,+,50M100S,60,0;chrT,{E2[0] + 1},+,100S50M,60,0;", flag=2048),
           Read("v", E2[0], "100S50M", f"chrV,1,+,50M100S,60,0;chrT,{E1[1] - 49},+,50S50M50S,60,0;", flag=2048)]
    expect("a read whose primary is on another contig counts from its supplementary records, once",
           cs.scan([("chrT", vec)], ex), {("chrT", E1[1], E2[0]): {("v", True)}})
    p = (E1[1] - 60, cc.parse_cigar("60M40S"))
    expect("SA pieces on the same strand in read order pair",
           cs.junctions("chrT", *p, f"chrT,{E2[0] + 1},+,60S40M,60,0;", ex), {(E1[1], E2[0])})
    expect(f"{cc.EDGE_SLACK} read bases of overlap at the junction still pair",
           cs.junctions("chrT", *p, f"chrT,{E2[0] + 1},+,{60 - cc.EDGE_SLACK}S{40 + cc.EDGE_SLACK}M,60,0;", ex),
           {(E1[1], E2[0])})
    expect("an SA piece on the other strand does not pair",
           cs.junctions("chrT", *p, f"chrT,{E2[0] + 1},-,60S40M,60,0;", ex), set())
    expect("an SA piece earlier in the read than the exon-end piece does not pair",
           cs.junctions("chrT", *p, f"chrT,{E2[0] + 1},+,10S40M50S,60,0;", ex), set())
    expect("both reversed, in read order: they pair",
           cs.junctions("chrT", *p, f"chrT,{E2[0] + 1},-,60S40M,60,0;", ex, True), {(E1[1], E2[0])})
    two = Read("n1", E1[1] - 50, f"50M{E2[0] - E1[1]}N150M{E3[0] - E2[1]}N50M")
    expect("a read across two junctions counts on each junction and once in the locus total",
           [r[4:] for r in listed("spliced", [two])],
           [["10", f"{E1[1]}-{E2[0] + 1}:5,{E1[1]}-{E3[0] + 1}:2,{E2[1]}-{E3[0] + 1}:4"]])
    expect(f"an intron of MIN_INTRON ({cs.MIN_INTRON}) bases counts, one base less does not",
           (cs.gap_ok(1200, 1200 + cs.MIN_INTRON), cs.gap_ok(1200, 1199 + cs.MIN_INTRON)), (True, False))
    keep = cs.MIN_CLIP
    cs.MIN_CLIP = 101
    expect("control: with MIN_CLIP above the read length the spliced reads list only the N-gap read, so nothing",
           listed("spliced"), [])
    cs.MIN_CLIP = keep
    keep = cs.MIN_JUNCTIONS
    cs.MIN_JUNCTIONS = 1
    expect("with MIN_JUNCTIONS 1 the one-junction gene is listed",
           [r[3:5] for r in listed("one")], [["1", "4"]])
    cs.MIN_JUNCTIONS = keep


def write_bam(pysam, path, kind):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chrT", "LN": LENGTH}]}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, start, cigar, sa, flag in fixture(kind):
            a = pysam.AlignedSegment()
            a.query_name, a.flag = name, flag
            a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
            a.cigarstring = cigar
            a.query_sequence = "A" * 100
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.next_reference_id, a.next_reference_start, a.template_length = -1, -1, 0
            if sa:
                a.set_tag("SA", sa)
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)


def run_main(args):
    with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()):
        rc = cs.main(args)
    return rc, out.getvalue()


def part2(pysam):
    with tempfile.TemporaryDirectory() as tmp:
        exons = exons_file(tmp)
        for kind, want in (("spliced", [SPLICED_ROW]), ("genomic", [])):
            bam = os.path.join(tmp, f"{kind}.bam")
            write_bam(pysam, bam, kind)
            out = os.path.join(tmp, f"{kind}.tsv")
            rc, log = run_main(["--bam", bam, "--exons", exons, "--out", out])
            expect(f"{kind} BAM: exit code", rc, 0)
            lines = [line.rstrip("\n").split("\t") for line in open(out)]
            expect(f"{kind} BAM: the comment line", lines[0][0],
                   f"# cdna_scan: read 3 exon regions on 1 contigs; {len(want)} loci listed")
            expect(f"{kind} BAM: header", lines[1], cs.COLUMNS)
            expect(f"{kind} BAM: rows", lines[2:], want)
            expect(f"{kind} BAM: the log line", f"SPLICED READS chrT:{E1[1]}-{E3[0] + 1}: 3 exon junctions" in log,
                   bool(want))
        out = os.path.join(tmp, "skip.tsv")
        rc, _ = run_main(["--skip", "no BAM", "--out", out])
        expect("--skip: the reason and no rows", (rc, open(out).read().splitlines()),
               (0, ["# cdna_scan: not run: no BAM", "\t".join(cs.COLUMNS)]))
        rc, _ = run_main(["--bam", os.path.join(tmp, "spliced.bam"), "--exons", os.path.join(tmp, "none.bed"),
                          "--out", out])
        expect("a missing exon file exits 2", rc, 2)


def in_container():
    """Run part 2 in PYPGX_IMAGE: its exit code."""
    image = ""
    for line in open(os.path.join(REPO, "versions.env")):
        if line.startswith("PYPGX_IMAGE="):
            image = line.split("=", 1)[1].split("#")[0].strip().strip('"')
    if not image:
        print("FAIL no PYPGX_IMAGE line in versions.env")
        return 1
    cmd = ["docker", "run", "--rm", "--network", "none", "-u", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp",
           "-e", "CDNA_SCAN_TEST_PART=2", "-v", f"{REPO}:/repo:ro", image, "python3", "/repo/tests/test_cdna_scan.py"]
    print("+ " + " ".join(cmd), flush=True)
    try:
        return subprocess.run(cmd).returncode
    except FileNotFoundError:
        print("FAIL part 2 needs pysam or docker, and this machine has neither")
        return 1


def main():
    part = os.environ.get("CDNA_SCAN_TEST_PART")
    if part != "2":
        part1()
    try:
        import pysam
    except ImportError:
        if part == "2":
            print("FAIL pysam is not in this image")
            return 1
        if in_container() != 0:
            FAILS.append("part 2 (in PYPGX_IMAGE)")
    else:
        part2(pysam)
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
