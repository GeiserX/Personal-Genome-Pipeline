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

Also: the check file has the status, the message, the four depths and the
MAPQ shares, the `bed` command prints CYP2D6 and both flanks, and an input
without a CYP2D6 row, or with regions other than those (a zero-length one
included), or an output that cannot be written, exits 2.

The MAPQ shares beside the verdict (`measure`, then `check --measure`): on
synthetic reads at CYP2D6, reads in REP6 at MAPQ 0 lower the share over
Cyrius's region but not over the gene body; a MAPQ 10 read counts at MAPQ >= 1
and not at >= 20; a duplicate counts nowhere; without the measure file every
share it gives reads 'not measured'; the shares never change the status.

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
    expect("check file: header and the eleven rows",
           list(rows), ["metric", "status", "message", "gene_depth_all", "gene_depth_mapq1",
                        "flank_depth_all", "flank_depth_mapq1", "gene_mapq1_share", "gene_mapq20_share",
                        "body_depth_all", "body_mapq1_share", "body_mapq20_share"])
    expect("check file: the verdict's MAPQ >= 1 share over Cyrius's region (16.79 / 18.53)",
           rows.get("gene_mapq1_share"), "0.91")
    expect("check file: without `measure` the other shares are not measured",
           [rows.get(k) for k in ("gene_mapq20_share", "body_depth_all", "body_mapq1_share", "body_mapq20_share")],
           ["not measured"] * 4)
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

# --- the MAPQ shares: `measure` on synthetic reads, then `check --measure` ------------
GENE = next((s, e) for c, s, e, n in dc.REGIONS if n == "CYP2D6")
B0, B1 = dc.BODY[1:]
L = 100
# Cyrius's region tiled once with reads (depth 1 everywhere): in the body at
# MAPQ 60, its first ten at MAPQ 10; outside the body (REP6 and the rest) at
# MAPQ 0; plus one duplicate read on the body.
def tile(s, e, mapq):
    return [(0, mapq, [(x, min(x + L, e))]) for x in range(s, e, L)]


reads = tile(GENE[0], B0, 0) + tile(B0, B0 + 10 * L, 10) + tile(B0 + 10 * L, B1, 60) + tile(B1, GENE[1], 0)
reads.append((0x400, 60, [(B0, B0 + L)]))
m = {(n, q): d for n, q, d in dc.measure_rows(reads)}
expect("measure: depth 1 over the body and the region, the duplicate counts nowhere",
       (round(m[("body", 0)], 6), round(m[("gene", 0)], 6)), (1.0, 1.0))
expect("measure: the MAPQ 0 reads in REP6 leave the body's MAPQ >= 1 share at 1",
       m[("body", 1)], m[("body", 0)])
expect("measure: a MAPQ 10 read counts at MAPQ >= 1, not at >= 20",
       round(m[("body", 0)] - m[("body", 20)], 6), round(10 * L / (B1 - B0), 6))
expect("measure: over Cyrius's region the MAPQ 0 reads lower the MAPQ >= 1 share to the body's part",
       round(m[("gene", 1)], 6), round((B1 - B0) / (GENE[1] - GENE[0]), 6))

with tempfile.TemporaryDirectory() as tmp:
    meas = os.path.join(tmp, "shares.tsv")
    with open(meas, "w") as f:
        f.write("region\tmin_mapq\tmean_depth\n")
        f.writelines(f"{n}\t{q}\t{d:.4f}\n" for (n, q), d in m.items())
    out = os.path.join(tmp, "check.tsv")
    with contextlib.redirect_stdout(io.StringIO()):
        rc = dc.main(["check", "--all", os.path.join(DEPTH, "clean_no_alt.q0.bed"),
                      "--mapq1", os.path.join(DEPTH, "clean_no_alt.q1.bed"), "--measure", meas, "--out", out])
    rows = dict(line.rstrip("\n").split("\t", 1) for line in open(out)) if rc == 0 else {}
    expect("check --measure: exit code", rc, 0)
    expect("check --measure: the shares from the measure file",
           [rows.get(k) for k in ("gene_mapq20_share", "body_mapq1_share", "body_mapq20_share")],
           [f"{m[('gene', 20)] / m[('gene', 0)]:.2f}", "1.00", f"{m[('body', 20)] / m[('body', 0)]:.2f}"])
    expect("check --measure: the body depth", rows.get("body_depth_all"), f"{m[('body', 0)]:.2f}")
    expect("check --measure: the shares do not change the status", rows.get("status"), "ok")
    with contextlib.redirect_stdout(io.StringIO()):
        dc.main(["check", "--all", os.path.join(DEPTH, "multimapped_with_alt.q0.bed"),
                 "--mapq1", os.path.join(DEPTH, "multimapped_with_alt.q1.bed"), "--measure", meas, "--out", out])
    expect("check --measure: good shares do not rescue a multi-mapped verdict",
           dict(line.rstrip("\n").split("\t", 1) for line in open(out)).get("status"), "unreliable")
    with open(meas, "w") as f:
        f.write("chrom\tstart\n")
    with contextlib.redirect_stderr(io.StringIO()):
        rc = dc.main(["check", "--all", os.path.join(DEPTH, "clean_no_alt.q0.bed"),
                      "--mapq1", os.path.join(DEPTH, "clean_no_alt.q1.bed"), "--measure", meas, "--out", out])
    expect("check --measure: a file that is not measure's exits 2", rc, 2)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
