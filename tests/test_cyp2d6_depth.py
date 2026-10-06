#!/usr/bin/env python3
"""bin/cyp2d6_depth_check.py: is the depth at CYP2D6 fit for a copy-number call?

Five depth tables (tests/fixtures/pgx/depth/<case>.q0.bed and .q1.bed, the
mosdepth region files of all reads and of MAPQ >= 1):

  clean_no_alt          the fixture's CYP2D slice mapped to the no-ALT
                        reference (ALT depth A/B run 37177862776): ok
  multimapped_with_alt  the same reads mapped to the Broad hg38 FASTA with
                        ALT contigs (same run): unreliable, multi-mapped
  deleted_both          CYP2D6 deleted on both copies (*5/*5): almost no reads
                        at the gene, normal flanks: ok, the deletion is real
  deleted_one           one copy deleted: half the flank depth: ok
  no_flank_reads        no reads at all: unreliable

Also: the check file has the status, the message and the four depths, the
`bed` command prints CYP2D6 and both flanks, and an input without a CYP2D6
row, or with regions other than those (a zero-length one included), or an
output that cannot be written, exits 2.

Run: python3 tests/test_cyp2d6_depth.py
"""
import contextlib
import io
import os
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import cyp2d6_depth_check as dc  # noqa: E402

DEPTH = os.path.join(REPO, "tests", "fixtures", "pgx", "depth")
FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def run_check(case, tmp):
    out = os.path.join(tmp, f"{case}.tsv")
    with contextlib.redirect_stdout(io.StringIO()):
        rc = dc.main(["check", "--all", os.path.join(DEPTH, f"{case}.q0.bed"),
                      "--mapq1", os.path.join(DEPTH, f"{case}.q1.bed"), "--out", out])
    rows = dict(line.rstrip("\n").split("\t", 1) for line in open(out)) if rc == 0 else {}
    return rc, rows


CASES = [
    ("clean_no_alt", "ok", "CYP2D6 depth is fit for a copy-number call"),
    ("multimapped_with_alt", "unreliable", dc.UNRELIABLE),
    ("deleted_both", "ok", "CYP2D6 depth is fit for a copy-number call"),
    ("deleted_one", "ok", "CYP2D6 depth is fit for a copy-number call"),
    ("no_flank_reads", "unreliable", dc.NO_DEPTH),
]

with tempfile.TemporaryDirectory() as tmp:
    for case, status, message in CASES:
        rc, rows = run_check(case, tmp)
        expect(f"{case}: exit code", rc, 0)
        expect(f"{case}: status", rows.get("status"), status)
        expect(f"{case}: message", rows.get("message"), message)

    # The bead's wording, word for word, is in the message a multi-mapped BAM gets.
    expect("the multi-mapped message names the cause", dc.UNRELIABLE,
           "CYP2D6 copy number unreliable: reads are multi-mapped "
           "(was this BAM aligned to a reference with ALT contigs?)")

    rc, rows = run_check("clean_no_alt", tmp)
    expect("check file: header and the six rows",
           list(rows), ["metric", "status", "message", "gene_depth_all", "gene_depth_mapq1",
                        "flank_depth_all", "flank_depth_mapq1"])
    expect("check file: depths, length-weighted",
           [rows.get(k) for k in ("gene_depth_all", "gene_depth_mapq1", "flank_depth_all", "flank_depth_mapq1")],
           ["18.53", "16.79", "29.68", "29.67"])

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        dc.main(["bed"])
    names = [line.split("\t")[3] for line in buf.getvalue().splitlines()]
    expect("bed: CYP2D6 between two flanks", names, ["flank", "CYP2D6", "flank"])

    bad = os.path.join(tmp, "bad.bed")
    with open(bad, "w") as f:
        f.write("chr22\t42050000\t42100000\tflank\t30\n")
    with contextlib.redirect_stderr(io.StringIO()):
        rc = dc.main(["check", "--all", bad, "--mapq1", bad, "--out", os.path.join(tmp, "x.tsv")])
    expect("an input without a CYP2D6 row exits 2", rc, 2)

    moved = os.path.join(tmp, "moved.bed")
    with open(os.path.join(DEPTH, "clean_no_alt.q0.bed")) as f, open(moved, "w") as out:
        out.write(f.read().replace("42123192", "42123000"))
    with contextlib.redirect_stderr(io.StringIO()):
        rc = dc.main(["check", "--all", moved, "--mapq1", moved, "--out", os.path.join(tmp, "y.tsv")])
    expect("regions other than the ones `bed` prints exit 2", rc, 2)

    clean = os.path.join(DEPTH, "clean_no_alt.q0.bed")
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        rc = dc.main(["check", "--all", clean, "--mapq1", clean, "--out", os.path.join(tmp, "no", "such", "dir.tsv")])
    expect("an output that cannot be written exits 2", rc, 2)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
