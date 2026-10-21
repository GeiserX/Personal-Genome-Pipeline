#!/usr/bin/env python3
"""bin/rfc1_motif.py: which motif do the reads at RFC1 carry?

Every case is a set of made-up reads at the GRCh38 RFC1 repeat
(chr4:39348425-39348479, the reference (AAAAG)11) with random flanks:
  benign      reads across the repeat with an expanded AAAAG allele, and
              reads whose AAAAG part is soft-clipped: benign_motif
  pathogenic  AAGGG reads of four kinds: soft-clipped next to the repeat,
              unmapped and placed at their mate's position, mates mapped to
              another contig (mate rescue), and reverse-strand CCCTT reads
              across the repeat; three reads of an expanded AAAAG allele;
              and eight reads across an unexpanded (AAAAG)11 allele, which
              are not counted: pathogenic_motif, with AAGGG and AAAAG shares
  mixed       AAAAG and AAAGG reads only: benign_motif, both shares listed
  acagg       ACAGG reads: pathogenic_motif
  unknown     AAGAG reads (not a known benign or pathogenic motif):
              uncertain_motif
  few         two repeat reads: not_resolved
and duplicates, secondary and QC-failed alignments are not counted.

Controls, each a check that fails when the rescue it names is off:
  - the pathogenic case counts its 6 soft-clipped reads; with the soft-clip
    rescue off (--no-softclip) it counts fewer;
  - it counts its 3 rescued mates; with --no-mates fewer;
  - a homopolymer and a too-short run are not counted.

Part 1 runs the counting in plain Python on read objects. Part 2 writes the
same reads to a BAM with pysam and runs the script's command line on it,
as step 9b does. Without pysam on the host, part 2 runs in PYPGX_IMAGE
(versions.env) through docker.

Run: python3 tests/test_rfc1_motif.py
"""
import contextlib
import io
import os
import random
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import rfc1_motif as rm  # noqa: E402

FAILS = []
CHROM, START, END = "chr4", 39348424, 39348479  # 0-based start, end exclusive: 55 bases, (AAAAG)11
OTHER = "chr1"
RNG = random.Random(7)
LEFT = "".join(RNG.choice("ACGT") for _ in range(300))
RIGHT = "".join(RNG.choice("ACGT") for _ in range(300))
VCF = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
       f"{CHROM}\t{START + 1}\t.\tA\t<STR100>\t.\tPASS\tEND={END};REF=11;RL=55;RU=AARRG;VARID=RFC1;REPID=RFC1;"
       "STR_STATUS=full_mutation\tGT:SO:REPCN\t" + ":".join(["/".join("11"), "INREPEAT/INREPEAT", "/".join(["100", "150"])]) + "\n")


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def expect_true(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


class Read:
    """The pysam.AlignedSegment attributes rfc1_motif reads, and what part 2 writes."""

    def __init__(self, name, contig, start, cigar, seq, read1=True, unmapped=False, mate=None, reverse=False,
                 dup=False, secondary=False, qcfail=False):
        self.query_name, self.contig, self.reference_start = name, contig, start
        self.reference_id = 0 if contig == CHROM else 1
        self.cigar_text = cigar
        self.cigartuples = rm.parse_cigar(cigar) if cigar else None
        self.query_sequence = seq
        self.is_read1, self.is_read2 = read1, not read1
        self.is_unmapped, self.is_reverse = unmapped, reverse
        self.is_duplicate, self.is_secondary, self.is_supplementary, self.is_qcfail = dup, secondary, False, qcfail
        self.is_paired = True
        mate = mate or (contig, start)
        self.mate_contig, self.next_reference_start = mate
        self.next_reference_id = 0 if mate[0] == CHROM else 1
        self.mate_is_unmapped = False


def spanning(name, unit, units, offset=40, read1=True, **kw):
    """A 150-base read aligned across the repeat: offset flank bases, the
    repeat (55 reference bases plus an insertion of the rest), flank."""
    rep = unit * units
    right = 150 - offset - len(rep)
    assert right > 0 and len(rep) >= END - START
    ins = len(rep) - (END - START)
    cig = f"{offset}M{END - START}M" + (f"{ins}I" if ins else "") + f"{right}M"
    return Read(name, CHROM, START - offset, cig, LEFT[-offset:] + rep + RIGHT[:right], read1=read1, **kw)


def clipped(name, unit, read1=True):
    """A read whose 60 aligned bases end at the repeat start; the 90 repeat bases are soft-clipped."""
    seq = LEFT[-60:] + (unit * 18)
    return Read(name, CHROM, START - 60, "60M90S", seq, read1=read1)


def placed_unmapped(name, unit):
    """An unmapped read placed at its mapped mate's position, all repeat."""
    return Read(name, CHROM, START - 100, "", (unit * 30), read1=False, unmapped=True)


def rescued(name, unit):
    """A flank read at the locus and its repeat mate mapped to another contig."""
    a = Read(name, CHROM, START - 120, "100M", LEFT[-120:-20], read1=True, mate=(OTHER, 5000))
    b = Read(name, OTHER, 5000, "150M", unit * 30, read1=False, mate=(CHROM, START - 120))
    return a, b


def revcomp(s):
    return rm.revcomp(s)


def case(kind):
    """(reads at the locus, mates elsewhere by name) of one case."""
    at, elsewhere = [], {}
    if kind == "benign":
        at += [spanning(f"b{i}", "AAAAG", 20) for i in range(6)]
        at += [clipped(f"bc{i}", "AAAAG") for i in range(4)]
    elif kind == "pathogenic":
        at += [clipped(f"pc{i}", "AAGGG") for i in range(6)]
        at += [placed_unmapped(f"pu{i}", "AAGGG") for i in range(2)]
        for i in range(3):
            a, b = rescued(f"pm{i}", "AAGGG")
            at.append(a)
            elsewhere[a.query_name] = b
        at += [Read(f"pr{i}", CHROM, START - 40, "150M", revcomp(LEFT[-40:] + "AAGGG" * 22), reverse=True)
               for i in range(2)]
        at += [spanning(f"aaaag{i}", "AAAAG", 20) for i in range(3)]
        # Not counted: reads across an unexpanded (AAAAG)11 allele, and one
        # across 12 AAGGG units would be (more than the reference's 11).
        at += [spanning(f"ref{i}", "AAAAG", 11) for i in range(8)]
        # Not counted: a duplicate, a secondary and a QC-failed alignment.
        at.append(clipped("pc0", "AAGGG", read1=True))  # the same read again
        at.append(Read("dup", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, dup=True))
        at.append(Read("sec", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, secondary=True))
        at.append(Read("qc", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 18, qcfail=True))
        # Not counted either: a homopolymer and a run of three units.
        at.append(Read("poly", CHROM, START - 60, "60M90S", LEFT[-60:] + "A" * 90))
        at.append(Read("short", CHROM, START - 60, "60M90S", LEFT[-60:] + "AAGGG" * 3 + RIGHT[:75]))
    elif kind == "mixed":
        at += [spanning(f"m{i}", "AAAAG", 20) for i in range(6)]
        at += [spanning(f"n{i}", "AAAGG", 20) for i in range(4)]
    elif kind == "acagg":
        at += [clipped(f"a{i}", "ACAGG") for i in range(5)]
    elif kind == "unknown":
        at += [clipped(f"u{i}", "AAGAG") for i in range(5)]
    elif kind == "few":
        at += [clipped(f"f{i}", "AAGGG") for i in range(2)]
    return at, elsewhere


def run_plain(kind, softclip=True, mates=True):
    at, elsewhere = case(kind)
    lo, hi = START - rm.FLANK, END + rm.FLANK
    mate_of = (lambda r: elsewhere.get(r.query_name)) if mates else None
    units = rm.min_units(START, END)
    counts = rm.count_motifs(rm.read_segments(at, lo, hi, mate_of, softclip), units)
    return counts, rm.verdict(counts, units)


def part1():
    expect("longest run: (AAAAG)11 in flanks", rm.longest_run(LEFT[-30:] + "AAAAG" * 11 + RIGHT[:30])[1] >= 11, True)
    expect("motif name: GGGAA and CCCTT are AAGGG", (rm.motif_name("GGGAA"), rm.motif_name("CCCTT")), ("AAGGG", "AAGGG"))
    expect("motif name: an unknown motif by its smallest form", rm.motif_name("GAGAA"), "AAGAG")
    expect("a homopolymer is not a pentamer run", rm.longest_run("A" * 60), (None, 0))
    r = clipped("x", "AAGGG")
    lo, hi = START - rm.FLANK, END + rm.FLANK
    expect("soft clip on: the clipped repeat bases are kept", rm.repeat_bases(r, lo, hi, True).endswith("AAGGG" * 18), True)
    expect("soft clip off: only the aligned flank bases", rm.repeat_bases(r, lo, hi, False), LEFT[-rm.FLANK:])

    counts, v = run_plain("benign")
    expect("benign: 10 AAAAG reads", counts, {"AAAAG": 10})
    expect("benign: benign_motif", v[4], "benign_motif")
    expect_true("benign: the note says the flag is not supported", "the flag is not supported" in v[5], v[5])

    counts, v = run_plain("pathogenic")
    expect("pathogenic: 13 AAGGG reads (6 clipped, 2 placed unmapped, 3 rescued mates, 2 reverse) and 3 expanded "
           "AAAAG; the 8 reference-allele reads are not counted", counts, {"AAGGG": 13, "AAAAG": 3})
    at, _ = case("pathogenic")
    lo, hi = START - rm.FLANK, END + rm.FLANK
    expect("control: with a 4-unit floor the 8 reference-allele reads would count",
           rm.count_motifs(rm.read_segments(at, lo, hi, None, True), rm.MIN_UNITS).get("AAAAG"), 3 + 8)
    expect("the floor is one unit more than the reference allele", rm.min_units(START, END), 12)
    expect("pathogenic: pathogenic_motif", v[4], "pathogenic_motif")
    expect("pathogenic: the shares of both motifs", v[3], "AAGGG=0.81,AAAAG=0.19")
    expect_true("pathogenic: the note asks for repeat-primed PCR", "repeat-primed PCR" in v[5], v[5])

    counts, _ = run_plain("pathogenic", softclip=False)
    expect_true("control: with the soft-clip rescue off the 6 clipped reads are lost",
                counts.get("AAGGG", 0) == 13 - 6, counts)
    counts, _ = run_plain("pathogenic", mates=False)
    expect_true("control: with the mate rescue off the 3 rescued mates are lost", counts.get("AAGGG", 0) == 13 - 3, counts)

    counts, v = run_plain("mixed")
    expect("mixed: AAAAG and AAAGG counted apart", counts, {"AAAAG": 6, "AAAGG": 4})
    expect("mixed: benign_motif with both shares", (v[4], v[3]), ("benign_motif", "AAAAG=0.60,AAAGG=0.40"))
    expect("acagg: pathogenic_motif", run_plain("acagg")[1][4], "pathogenic_motif")
    expect("unknown: uncertain_motif", run_plain("unknown")[1][:2] + (run_plain("unknown")[1][4],),
           ("AAGAG", "unknown", "uncertain_motif"))
    expect("few: not_resolved", run_plain("few")[1][4], "not_resolved")


def write_bam(pysam, path, kind):
    at, elsewhere = case(kind)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": 190214555}, {"SN": OTHER, "LN": 248956422}]}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for r in at + list(elsewhere.values()):
            a = pysam.AlignedSegment()
            a.query_name = r.query_name
            a.query_sequence = r.query_sequence
            a.query_qualities = pysam.qualitystring_to_array("I" * len(r.query_sequence))
            a.flag = (1 | (64 if r.is_read1 else 128) | (4 if r.is_unmapped else 0) | (16 if r.is_reverse else 0)
                      | (256 if r.is_secondary else 0) | (512 if r.is_qcfail else 0) | (1024 if r.is_duplicate else 0))
            a.reference_id = r.reference_id
            a.reference_start = r.reference_start
            a.next_reference_id = r.next_reference_id
            a.next_reference_start = r.next_reference_start
            if not r.is_unmapped:
                a.cigarstring = r.cigar_text
                a.mapping_quality = 60
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.remove(unsorted)


def run_main(args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = rm.main(args)
    return rc, buf.getvalue()


def table(path):
    lines = [line.rstrip("\n").split("\t") for line in open(path)]
    return dict(zip(lines[0], lines[1]))


def part2(pysam):
    with tempfile.TemporaryDirectory() as tmp:
        vcf = os.path.join(tmp, "S_eh_stranger.vcf")
        with open(vcf, "w") as fh:
            fh.write(VCF)
        out = os.path.join(tmp, "out.tsv")
        for kind, verdict, reads in (("benign", "benign_motif", "10"), ("pathogenic", "pathogenic_motif", "16"),
                                     ("mixed", "benign_motif", "10"), ("few", "not_resolved", "2")):
            bam = os.path.join(tmp, f"{kind}.bam")
            write_bam(pysam, bam, kind)
            rc, _ = run_main(["--bam", bam, "--vcf", vcf, "--out", out])
            t = table(out)
            expect(f"BAM {kind}: exit 0, {verdict}, {reads} repeat reads",
                   (rc, t["verdict"], t["repeat_reads"]), (0, verdict, reads))
        bam = os.path.join(tmp, "pathogenic.bam")
        run_main(["--bam", bam, "--vcf", vcf, "--out", out])
        t = table(out)
        expect("BAM pathogenic: the locus from the VCF, AAGGG dominant, both shares",
               (t["chrom"], t["start"], t["end"], t["dominant_motif"], t["motif_shares"]),
               (CHROM, str(START + 1), str(END), "AAGGG", "AAGGG=0.81,AAAAG=0.19"))
        run_main(["--bam", bam, "--vcf", vcf, "--out", out, "--no-softclip"])
        expect_true("BAM control: --no-softclip loses the 6 clipped reads", table(out)["repeat_reads"] == "10",
                    table(out))
        run_main(["--bam", bam, "--vcf", vcf, "--out", out, "--no-mates"])
        expect_true("BAM control: --no-mates loses the 3 rescued mates", table(out)["repeat_reads"] == "13", table(out))
        rc, _ = run_main(["--vcf", vcf, "--out", out, "--skip", "no BAM"])
        expect("--skip: not_checked with the reason", (rc, table(out)["verdict"], table(out)["note"]),
               (0, "not_checked", "no BAM"))
        novcf = os.path.join(tmp, "none.vcf")
        with open(novcf, "w") as fh:
            fh.write(VCF.replace("REPID=RFC1", "REPID=HTT").replace("VARID=RFC1", "VARID=HTT"))
        run_main(["--bam", bam, "--vcf", novcf, "--out", out])
        expect("a VCF without an RFC1 record: not_checked", table(out)["verdict"], "not_checked")


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
           "-e", "RFC1_MOTIF_TEST_PART=2", "-v", f"{REPO}:/repo:ro", image, "python3", "/repo/tests/test_rfc1_motif.py"]
    print("+ " + " ".join(cmd), flush=True)
    try:
        return subprocess.run(cmd).returncode
    except FileNotFoundError:
        print("FAIL part 2 needs pysam or docker, and this machine has neither")
        return 1


def main():
    part = os.environ.get("RFC1_MOTIF_TEST_PART")
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
