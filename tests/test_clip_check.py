#!/usr/bin/env python3
"""bin/clip_check.py: are a ClinVar hit's alt reads spliced cDNA?

Part 1, plain Python (reads as small objects, no BAM):
  - CIGAR and SA parsing, clip lengths, the exon-edge distance (EDGE_SLACK in,
    one more base out, the wrong side out);
  - the allele a read shows for a SNV, a deletion and an insertion;
  - duplicate, secondary and repeated alignments are not counted;
  - the verdict: flagged, clear when the ref reads are clipped too, not
    checked with too few alt reads.

Part 2, a tiny hand-made BAM written with pysam, a made-up 10 kb contig with
two exons, one SNV in the first exon near its end:
  spliced   8 alt reads clipped at the exon end with the rest on the second
            exon (SA tag), 1 hard-clipped supplementary alt alignment, 6
            genomic ref reads running on into the intron: flagged
  genomic   the same SNV on 8 genomic alt reads, the same ref reads: clear
  both      the spliced alt reads, and ref reads clipped there too: clear
  and the output table, --skip, a contig the BAM lacks, an indel hit, and a
  control: with ALT_SHARE_MIN above 1 the spliced fixture is not flagged, so
  the flagged check depends on the threshold.

Part 2 needs pysam. Without it on the host, the test runs part 2 in
PYPGX_IMAGE (versions.env) through docker, as step 06 runs the script.

Run: python3 tests/test_clip_check.py
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import clip_check as cc  # noqa: E402

FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


class Read:
    """The pysam.AlignedSegment attributes clip_check reads."""

    def __init__(self, name, start, cigar, seq, sa="", read1=True, flag_dup=False, secondary=False):
        self.query_name, self.reference_start = name, start
        self.cigartuples, self.query_sequence = cc.parse_cigar(cigar), seq
        self.is_read1, self.is_read2 = read1, not read1
        self.is_unmapped = self.is_qcfail = False
        self.is_duplicate, self.is_secondary = flag_dup, secondary
        self._sa = sa

    def has_tag(self, t):
        return t == "SA" and bool(self._sa)

    def get_tag(self, t):
        return self._sa


def genome(n=10000):
    """A made-up contig: the same bases every run."""
    s, out = 7, []
    for _ in range(n):
        s = (s * 75 + 74) % 65537
        out.append("ACGT"[s % 4])
    return "".join(out)


G = genome()
EXON_A, EXON_B = (1000, 1200), (3000, 3200)
POS0 = 1150
REF = G[POS0]
ALT = "ACGT"[("ACGT".index(REF) + 1) % 4]


def with_alt(seq, start):
    i = POS0 - start
    return seq[:i] + ALT + seq[i + 1:]


def exons_file(tmp):
    p = os.path.join(tmp, "exons.bed")
    with open(p, "w") as f:
        f.write(f"chrT\t{EXON_A[0]}\t{EXON_A[1]}\nchrT\t{EXON_B[0]}\t{EXON_B[1]}\n")
    return p


def fixture_reads(kind):
    """(name, start, cigar, seq, sa, flag) for the three fixtures."""
    reads = []
    ref_clipped = kind == "both"
    for i in range(6):
        start = 1112 + 2 * i
        if ref_clipped:
            aligned = EXON_A[1] - start
            reads.append((f"ref{i}", start, f"{aligned}M{100 - aligned}S",
                          G[start:EXON_A[1]] + G[EXON_B[0]:EXON_B[0] + 100 - aligned],
                          f"chrT,{EXON_B[0] + 1},+,{aligned}S{100 - aligned}M,60,0;", 0))
        else:
            reads.append((f"ref{i}", start, "100M", G[start:start + 100], "", 0))
    for i in range(8):
        start = 1110 + 2 * i
        if kind == "genomic":
            reads.append((f"alt{i}", start, "100M", with_alt(G[start:start + 100], start), "", 0))
        else:
            aligned = EXON_A[1] - start
            seq = with_alt(G[start:EXON_A[1]], start) + G[EXON_B[0]:EXON_B[0] + 100 - aligned]
            reads.append((f"alt{i}", start, f"{aligned}M{100 - aligned}S", seq,
                          f"chrT,{EXON_B[0] + 1},+,{aligned}S{100 - aligned}M,60,0;", 0))
    if kind != "genomic":
        # A supplementary alignment (flag 2048) of a read whose primary lies on
        # the second exon: hard-clipped, as aligners write supplementary ones.
        reads.append(("suppl", 1140, "60M40H", with_alt(G[1140:EXON_A[1]], 1140),
                      f"chrT,{EXON_B[0] + 1},+,60S40M,60,0;", 2048))
    return reads


def part1():
    expect("parse_cigar", cc.parse_cigar("5S70M2I3D20M10H"),
           [(cc.SOFT, 5), (cc.M, 70), (cc.INS, 2), (cc.DEL, 3), (cc.M, 20), (cc.HARD, 10)])
    expect("segment: start, end, left and right clip", cc.segment(100, cc.parse_cigar("5S70M2I3D20M10H")),
           (100, 193, 5, 10))
    expect("sa_segments: 1-based SA position to 0-based", cc.sa_segments("chr2,3001,-,80S20M,60,0;"),
           [("chr2", (3000, 3020, 80, 0))])

    with tempfile.TemporaryDirectory() as tmp:
        ex = cc.Exons(exons_file(tmp))
    end = EXON_A[1]
    expect("a right clip at the exon end counts", ex.edge_clipped("chrT", (1100, end, 0, 20)), True)
    expect(f"a right clip {cc.EDGE_SLACK} bp from the exon end counts",
           ex.edge_clipped("chrT", (1100, end + cc.EDGE_SLACK, 0, 20)), True)
    expect(f"a right clip {cc.EDGE_SLACK + 1} bp from the exon end does not",
           ex.edge_clipped("chrT", (1100, end + cc.EDGE_SLACK + 1, 0, 20)), False)
    expect(f"a clip of {cc.MIN_CLIP - 1} bases does not", ex.edge_clipped("chrT", (1100, end, 0, cc.MIN_CLIP - 1)),
           False)
    expect("a left clip at an exon end (the wrong side) does not", ex.edge_clipped("chrT", (end, 1300, 20, 0)), False)
    expect("a left clip at an exon start counts", ex.edge_clipped("chrT", (EXON_B[0], 3100, 20, 0)), True)
    expect("an SA alignment clipped at an exon start counts",
           ex.read_clipped("chrT", (1100, 1190, 0, 5), f"chrT,{EXON_B[0] + 1},+,80S20M,60,0;"), True)
    expect("another contig has no edges", ex.edge_clipped("chrQ", (1100, end, 0, 20)), False)

    seq = G[1100:1200]
    b = cc.walk(1100, cc.parse_cigar("100M"), seq)
    expect("SNV: the ref base", cc.allele(*b, POS0, REF, ALT), "ref")
    b = cc.walk(1100, cc.parse_cigar("100M"), with_alt(seq, 1100))
    expect("SNV: the alt base", cc.allele(*b, POS0, REF, ALT), "alt")
    b = cc.walk(1100, cc.parse_cigar("40M"), seq[:40])
    expect("SNV: a read that ends before the site", cc.allele(*b, POS0, REF, ALT), None)
    dref = G[POS0:POS0 + 4]
    b = cc.walk(1100, cc.parse_cigar("51M3D46M"), G[1100:POS0 + 1] + G[POS0 + 4:1200])
    expect("deletion: a read with the 3 bp deletion", cc.allele(*b, POS0, dref, dref[0]), "alt")
    b = cc.walk(1100, cc.parse_cigar("100M"), seq)
    expect("deletion: a read without it", cc.allele(*b, POS0, dref, dref[0]), "ref")
    b = cc.walk(1100, cc.parse_cigar("51M2I49M"), G[1100:POS0 + 1] + "TT" + G[POS0 + 1:1149 + 51])
    expect("insertion: a read with the inserted TT", cc.allele(*b, POS0, REF, REF + "TT"), "alt")
    b = cc.walk(1100, cc.parse_cigar("100M"), seq)
    expect("insertion: a read without it", cc.allele(*b, POS0, REF, REF + "TT"), "ref")
    expect("a complex allele is not checked", cc.complex_allele("AC", "GTT"), True)
    expect("a SNV is checked", cc.complex_allele("A", "G"), False)

    with tempfile.TemporaryDirectory() as tmp:
        ex = cc.Exons(exons_file(tmp))
    reads = [Read(n, s, c, q, sa, secondary=bool(f & 256)) for n, s, c, q, sa, f in fixture_reads("spliced")]
    expect("spliced, as objects: alt reads, clipped, ref reads, clipped",
           cc.count_reads(reads, ex, "chrT", POS0, REF, ALT), [9, 9, 6, 0])
    extra = [Read("alt0", 1110, "100M", with_alt(G[1110:1210], 1110)),
             Read("dup", 1110, "100M", with_alt(G[1110:1210], 1110), flag_dup=True),
             Read("sec", 1110, "100M", with_alt(G[1110:1210], 1110), secondary=True)]
    expect("a repeated name, a duplicate and a secondary alignment are not counted",
           cc.count_reads(reads + extra, ex, "chrT", POS0, REF, ALT), [9, 9, 6, 0])
    expect("verdict: most alt reads clipped, no ref read", cc.verdict(9, 9, 6, 0)[0], "flagged")
    expect("verdict: the ref reads are clipped too", cc.verdict(9, 9, 6, 6)[0], "clear")
    expect("verdict: few alt reads clipped", cc.verdict(9, 2, 6, 0)[0], "clear")
    expect("verdict: one alt read", cc.verdict(1, 1, 6, 0)[0], "not_checked")
    expect("verdict: the flag text", cc.verdict(9, 9, 0, 0)[1],
           "possible cDNA/transgene contamination: check the reads")


def write_bam(pysam, path, kind):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chrT", "LN": len(G)}]}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, start, cigar, seq, sa, flag in fixture_reads(kind):
            a = pysam.AlignedSegment()
            a.query_name, a.flag = name, flag
            a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
            a.cigarstring = cigar
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.next_reference_id, a.next_reference_start, a.template_length = -1, -1, 0
            if sa:
                a.set_tag("SA", sa)
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)


def run_main(args):
    with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()):
        rc = cc.main(args)
    return rc, out.getvalue()


def part2(pysam):
    with tempfile.TemporaryDirectory() as tmp:
        exons = exons_file(tmp)
        hits = os.path.join(tmp, "hits.tsv")
        with open(hits, "w") as f:
            f.write("chrom\tpos\tref\talt\tgenotype\tclinvar_id\tgeneinfo\tclnsig\tclnrevstat\n")
            f.write(f"chrT\t{POS0 + 1}\t{REF}\t{ALT}\t0/1\t1\tGENET:1\tPathogenic\tx\n")
            f.write("chrQ\t500\tA\tG\t0/1\t2\tGENEQ:2\tPathogenic\tx\n")
            f.write(f"chrT\t{POS0 + 1}\t{REF}AC\tGT\t0/1\t3\tGENET:1\tPathogenic\tx\n")
        rows = {}
        for kind in ("spliced", "genomic", "both"):
            bam = os.path.join(tmp, f"{kind}.bam")
            write_bam(pysam, bam, kind)
            out = os.path.join(tmp, f"{kind}.tsv")
            rc, log = run_main(["--hits", hits, "--bam", bam, "--exons", exons, "--out", out])
            expect(f"{kind}: exit code", rc, 0)
            lines = [line.rstrip("\n").split("\t") for line in open(out)]
            expect(f"{kind}: header", lines[0], cc.COLUMNS)
            expect(f"{kind}: one row per hit", len(lines), 4)
            rows[kind] = lines[1]
            if kind == "spliced":
                expect("spliced: the other two hits are not checked, with the reason",
                       [(r[8], r[9]) for r in lines[2:]],
                       [("not_checked", "the BAM has no contig chrQ"),
                        ("not_checked", "not a SNV, MNV or simple indel")])
                expect("spliced: the log names the flagged hit", f"FLAGGED chrT:{POS0 + 1}" in log, True)
        expect("spliced fixture: flagged, 9 of 9 alt and 0 of 6 ref reads clipped",
               rows["spliced"][4:9], ["9", "9", "6", "0", "flagged"])
        expect("spliced fixture: the flag text", rows["spliced"][9], cc.FLAG_NOTE)
        expect("genomic fixture with the same SNV: clear", rows["genomic"][4:9], ["8", "0", "6", "0", "clear"])
        expect("ref reads clipped at the exon edge too: clear", rows["both"][4:9], ["9", "9", "6", "6", "clear"])

        # Control: the flagged check above depends on the alt-share threshold.
        keep = cc.ALT_SHARE_MIN
        cc.ALT_SHARE_MIN = 1.01
        out = os.path.join(tmp, "control.tsv")
        run_main(["--hits", hits, "--bam", os.path.join(tmp, "spliced.bam"), "--exons", exons, "--out", out])
        cc.ALT_SHARE_MIN = keep
        expect("control: with ALT_SHARE_MIN 1.01 the spliced fixture is not flagged",
               open(out).readlines()[1].split("\t")[8], "clear")

        out = os.path.join(tmp, "skip.tsv")
        rc, _ = run_main(["--hits", hits, "--skip", "no BAM", "--out", out])
        expect("--skip: every hit not checked, with the reason",
               [line.rstrip("\n").split("\t")[8:] for line in open(out)][1:], [["not_checked", "no BAM"]] * 3)
        rc, _ = run_main(["--hits", hits, "--bam", os.path.join(tmp, "spliced.bam"),
                          "--exons", os.path.join(tmp, "none.bed"), "--out", out])
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
           "-e", "CLIP_CHECK_TEST_PART=2", "-v", f"{REPO}:/repo:ro", image, "python3", "/repo/tests/test_clip_check.py"]
    print("+ " + " ".join(cmd), flush=True)
    try:
        return subprocess.run(cmd).returncode
    except FileNotFoundError:
        print("FAIL part 2 needs pysam or docker, and this machine has neither")
        return 1


def main():
    part = os.environ.get("CLIP_CHECK_TEST_PART")
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
