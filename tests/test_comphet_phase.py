#!/usr/bin/env python3
"""bin/comphet_phase.py: is a compound-het pair in cis or in trans?

The fixture is tests/fixtures/comphet/make_fixture.py: read pairs on a
made-up contig with three pairs, GENEC (an SNV and a deletion 3 bp apart,
both on one copy), GENET (two SNVs 400 bp apart, one on each copy, seen only
across mates) and GENEF (two SNVs 50 kb apart).

Part 1, plain Python: the SAM text is read into small read objects, so the
counting runs without pysam:
  - GENEC reads cis, GENET trans, GENEF unknown, with the read counts;
  - control: the fixture with GENEC's deletion moved to the other copy
    (--swap) reads trans, so the cis verdict comes from the alleles;
  - control: with JOIN_MATES off, GENET falls to unknown, so its trans verdict
    rests on joining mates;
  - the verdict thresholds, duplicates and low-MAPQ reads left out, a pair of
    mates that disagree, the per-variant PHASE column and the slivar pairs.
Part 2, the real route: the SAM becomes a sorted, indexed BAM with pysam and
main() runs on it, --skip, and the PHASE column. Without pysam on the host
part 2 runs in PYPGX_IMAGE (versions.env) through docker, as step 31 runs the
script.

Run: python3 tests/test_comphet_phase.py
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "fixtures", "comphet"))
import comphet_phase as cp  # noqa: E402
import make_fixture  # noqa: E402
from clip_check import parse_cigar  # noqa: E402

FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


class Read:
    """The pysam.AlignedSegment attributes comphet_phase reads, from a SAM line."""

    def __init__(self, line):
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        self.query_name, self.reference_name = f[0], f[2]
        self.reference_start, self.mapping_quality = int(f[3]) - 1, int(f[4])
        self.cigartuples, self.query_sequence = parse_cigar(f[5]), f[9]
        self.is_unmapped, self.is_read1 = bool(flag & 4), bool(flag & 64)
        self.is_secondary, self.is_qcfail = bool(flag & 256), bool(flag & 512)
        self.is_duplicate, self.is_supplementary = bool(flag & 1024), bool(flag & 2048)
        self.reference_end = self.reference_start + sum(n for op, n in self.cigartuples if op in (0, 2, 3, 7, 8))


class Bam:
    """pysam.AlignmentFile's references and fetch, over SAM text."""

    def __init__(self, text):
        self.reads = [Read(x) for x in text.splitlines() if x and not x.startswith("@")]
        self.references = sorted({r.reference_name for r in self.reads})

    def fetch(self, contig, start, end):
        return [r for r in self.reads if r.reference_name == contig
                and r.reference_start < end and r.reference_end > start]


def fixture(swap=False):
    seq = make_fixture.reference()
    vs = make_fixture.variants(seq, swap)
    return make_fixture.sam(seq, vs), make_fixture.comphet_vcf(vs)


def run(sam_text, vcf_text, tmp):
    vcf = os.path.join(tmp, "comphet.vcf")
    with open(vcf, "w") as fh:
        fh.write(vcf_text)
    return cp.phase(cp.read_pairs(vcf), cp.Reads(Bam(sam_text)), None)


def by_gene(rows):
    return {r[0]: (r[9], r[10], r[11]) for r in rows}


def part1():
    with tempfile.TemporaryDirectory() as tmp:
        sam, vcf = fixture()
        rows = run(sam, vcf, tmp)
        expect("three pairs from slivar's entries, each named once",
               [(r[0], r[2], r[5]) for r in rows],
               [("GENEC", "1000", "1003"), ("GENEF", "20000", "70000"), ("GENET", "5000", "5400")])
        g = by_gene(rows)
        expect("GENEC (SNV and deletion 3 bp apart on one copy): cis on 22 reads", g["GENEC"], ("cis", "22", "0"))
        expect("GENET (400 bp apart, one on each copy): trans on 20 read pairs", g["GENET"], ("trans", "0", "20"))
        expect("GENEF (50 kb apart): unknown, no read pair covers both", g["GENEF"], ("unknown", "0", "0"))
        expect("GENEF's note", [r[12] for r in rows if r[0] == "GENEF"], [
            "no read or read pair covers both positions"])
        expect("the distance column", [r[8] for r in rows], ["3", "50000", "400"])

        # Control 1: the same positions with the deletion on the other copy.
        sam_s, vcf_s = fixture(swap=True)
        expect("control: GENEC with the deletion moved to the other copy reads trans",
               by_gene(run(sam_s, vcf_s, tmp))["GENEC"], ("trans", "0", "22"))

        # Control 2: without joining mates no single read covers GENET's pair.
        keep, cp.JOIN_MATES = cp.JOIN_MATES, False
        try:
            g2 = by_gene(run(sam, vcf, tmp))
        finally:
            cp.JOIN_MATES = keep
        expect("control: with JOIN_MATES off GENET falls to unknown", g2["GENET"], ("unknown", "0", "0"))
        expect("control: GENEC still cis on single reads", g2["GENEC"][0], "cis")

        # Reads that do not count.
        lines = sam.splitlines()
        c_reads = [x for x in lines if x.startswith("c0_")]
        dup = [x.replace("\t99\t", "\t1123\t").replace("\t147\t", "\t1171\t") for x in c_reads]
        lowq = [x.replace("\t60\t", "\t10\t", 1) for x in c_reads]
        rest = "\n".join(x for x in lines if not x.startswith("c0_"))
        expect("duplicates do not count: GENEC on copy 2's 11 reads alone is cis",
               by_gene(run(rest + "\n" + "\n".join(dup) + "\n", vcf, tmp))["GENEC"], ("cis", "11", "0"))
        expect("MAPQ below 20 does not count",
               by_gene(run(rest + "\n" + "\n".join(lowq) + "\n", vcf, tmp))["GENEC"], ("cis", "11", "0"))

    # Verdicts.
    expect("verdict: none", cp.verdict(0, 0)[0], "unknown")
    expect("verdict: one read is not enough", cp.verdict(0, 1), ("unknown", "only 1 read or read pair covers both positions"))
    expect("verdict: 8 of 10 trans is trans", cp.verdict(2, 8)[0], "trans")
    expect("verdict: 7 of 10 either way is unknown", cp.verdict(3, 7), ("unknown", "the reads disagree: 3 cis, 7 trans"))
    expect("verdict: 9 of 10 cis is cis", cp.verdict(9, 1)[0], "cis")

    # Mates that disagree at one position are left out.
    class R:
        def __init__(self, name, read1, seq):
            self.query_name, self.is_read1 = name, read1
            self.reference_start, self.cigartuples, self.query_sequence = 99, [(0, 3)], seq
            self.mapping_quality = 60
            self.is_unmapped = self.is_secondary = self.is_supplementary = False
            self.is_duplicate = self.is_qcfail = False
    got = cp.site_alleles([R("a", True, "AGA"), R("a", False, "ACA"), R("b", True, "AGA")], ("chrX", 101, "C", "G"))
    expect("mates showing alt and ref at one site: the pair is left out", got, {"b": "alt"})
    got = cp.site_alleles([R("a", True, None), R("b", True, "AGA")], ("chrX", 101, "C", "G"))
    expect("a read with no stored sequence (SEQ *) is skipped, the rest still counted", got, {"b": "alt"})

    # Per-variant PHASE: trans over unknown over cis.
    rows = [["G", "c", "1", "A", "C", "5", "G", "T", ".", "cis", "3", "0", "."],
            ["G", "c", "1", "A", "C", "9", "T", "A", ".", "trans", "0", "3", "."],
            ["G", "c", "5", "G", "T", "9", "T", "A", ".", "unknown", "0", "0", "."]]
    expect("per-variant PHASE: trans with one partner wins, unknown over cis",
           cp.variant_phase(rows), {("c", "1", "A", "C"): "trans", ("c", "5", "G", "T"): "unknown",
                                    ("c", "9", "T", "A"): "trans"})


def write_bam(pysam, sam_text, path):
    with tempfile.TemporaryDirectory() as tmp:
        sam = os.path.join(tmp, "reads.sam")
        with open(sam, "w") as fh:
            fh.write(sam_text)
        unsorted = os.path.join(tmp, "unsorted.bam")
        with pysam.AlignmentFile(sam, "r") as src, pysam.AlignmentFile(unsorted, "wb", template=src) as out:
            for r in src:
                out.write(r)
        pysam.sort("-o", path, unsorted)
        pysam.index(path)


def run_main(args):
    with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()):
        rc = cp.main(args)
    return rc, out.getvalue()


def tsv(path):
    return [line.rstrip("\n").split("\t") for line in open(path)]


def part2(pysam):
    with tempfile.TemporaryDirectory() as tmp:
        make_fixture.write(tmp)
        bam = os.path.join(tmp, "reads.bam")
        write_bam(pysam, open(os.path.join(tmp, "reads.sam")).read(), bam)
        vcf = os.path.join(tmp, "comphet.vcf")
        variants = os.path.join(tmp, "variants.tsv")
        with open(variants, "w") as fh:
            fh.write("GENE\tCHROM\tPOS\tREF\tALT\tIMPACT\tConsequence\tGT\tMAX_AF\n")
            for line in open(vcf):
                if not line.startswith("#"):
                    f = line.split("\t")
                    fh.write(f"G\t{f[0]}\t{f[1]}\t{f[3]}\t{f[4]}\tHIGH\tstop_gained\t0/1\t0.0001\n")
        out, vout = os.path.join(tmp, "pairs.tsv"), os.path.join(tmp, "variants.out.tsv")
        rc, log = run_main(["--vcf", vcf, "--bam", bam, "--out", out, "--variants", variants, "--variants-out", vout])
        expect("BAM: exit code", rc, 0)
        rows = tsv(out)
        expect("BAM: header", rows[0], cp.COLUMNS)
        expect("BAM: GENEC cis, GENEF unknown, GENET trans, with the counts",
               [(r[0], r[9], r[10], r[11]) for r in rows[1:]],
               [("GENEC", "cis", "22", "0"), ("GENEF", "unknown", "0", "0"), ("GENET", "trans", "0", "20")])
        expect("BAM: the log counts them", "3 compound-het pairs from the reads: 1 trans, 1 cis, 1 unknown" in log, True)
        v = tsv(vout)
        expect("variants: PHASE goes right after GT", v[0][7:10], ["GT", "PHASE", "MAX_AF"])
        expect("variants: each variant's PHASE", [(r[2], r[8]) for r in v[1:]],
               [("1000", "cis"), ("1003", "cis"), ("5000", "trans"), ("5400", "trans"),
                ("20000", "unknown"), ("70000", "unknown")])

        rc, log = run_main(["--vcf", vcf, "--skip", "no BAM", "--out", out])
        expect("--skip: every pair unknown with the reason",
               {(r[9], r[10], r[12]) for r in tsv(out)[1:]}, {("unknown", ".", "no BAM")})
        rc, _ = run_main(["--vcf", os.path.join(tmp, "none.vcf"), "--bam", bam, "--out", out])
        expect("a missing VCF exits 2", rc, 2)

        # The swapped fixture through the BAM too.
        sw = os.path.join(tmp, "swap")
        make_fixture.write(sw, swap=True)
        write_bam(pysam, open(os.path.join(sw, "reads.sam")).read(), os.path.join(sw, "reads.bam"))
        run_main(["--vcf", os.path.join(sw, "comphet.vcf"), "--bam", os.path.join(sw, "reads.bam"), "--out", out])
        expect("BAM control: the swapped GENEC reads trans", [r[9] for r in tsv(out)[1:] if r[0] == "GENEC"], ["trans"])


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
           "-e", "COMPHET_TEST_PART=2", "-v", f"{REPO}:/repo:ro", image, "python3", "/repo/tests/test_comphet_phase.py"]
    print("+ " + " ".join(cmd), flush=True)
    try:
        return subprocess.run(cmd).returncode
    except FileNotFoundError:
        print("FAIL part 2 needs pysam or docker, and this machine has neither")
        return 1


def main():
    part = os.environ.get("COMPHET_TEST_PART")
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
