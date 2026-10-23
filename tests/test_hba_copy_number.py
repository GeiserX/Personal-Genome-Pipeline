#!/usr/bin/env python3
"""bin/hba_copy_number.py: step 35's HBA1/HBA2 count from short-read depth.

On `samtools bedcov` tables made up for each case (summed depth = mean depth
x length), at 30x in the flanks, with a quarter of the gene-body depth from
reads of MAPQ 0 (as BWA-MEM and minimap2 give on HG002) unless a case says
otherwise:
  1. no deletion: 30x over HBA2 and HBA1 is four alpha-globin genes;
  2. a planted -a3.7-like loss (HBA2 at half depth) is three, a silent carrier,
     marked unconfirmed (the count is checked on one real genome only);
  3. both genes at half depth (--SEA-like) is two, with the warning that
     depth cannot tell --/aa from -a/-a;
  4. a triplication (HBA2 at 1.5x) is five;
  5. an estimate between two whole numbers, flanks under 8x, and flanks that
     disagree are each 'not assessed', with the reason;
  6. no read of MAPQ 0 over the genes (GIAB's novoalign BAM of HG002: HBA2
     15.5x, HBA1 25.7x, flanks 30x and 31x, which main counted as three or
     left between two and three) is 'not assessed', naming the cause; 1%
     is too, 3% is called; no gene at all (no reads to judge by) is still 0;
  7. the CLI: `bed` prints the four regions, `check` writes the metric table
     and exits 0; a table with other regions, or a short row, in either
     input exits 2, and so does a missing --bedcov-mapq1.

Run: python3 tests/test_hba_copy_number.py
"""
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "hba_copy_number.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import hba_copy_number as hba  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {detail}"))
    if not ok:
        FAILS += 1


def bedcov(path, depth):
    """A bedcov table with the given mean depth per region name."""
    with open(path, "w") as f:
        for c, s, e, n in hba.REGIONS:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{int(depth[n] * (e - s))}\n")


def mapq1(depth, share=0.25):
    """The depths of reads of MAPQ 1 or more: the flanks whole, the genes
    without the given share."""
    return {n: v * (1 - share) if n in ("HBA2", "HBA1") else v for n, v in depth.items()}


def run(depth, share=0.25):
    with tempfile.TemporaryDirectory() as t:
        p, q = os.path.join(t, "cov.tsv"), os.path.join(t, "cov_q1.tsv")
        bedcov(p, depth)
        bedcov(q, mapq1(depth, share))
        return hba.assess(hba.read_bedcov(p), hba.read_bedcov(q))


def depths(hba2, hba1, left=30.0, right=30.0):
    return {"flank_left": left, "HBA2": hba2, "HBA1": hba1, "flank_right": right}


def main():
    r = run(depths(30, 30))
    check("no deletion: four alpha-globin genes, called",
          (r["status"], r["alpha_globin_genes"], r["estimate"]) == ("called", "4", "4.00")
          and "no deletion seen" in r["message"], r)
    r = run(depths(15, 30))
    check("planted one-gene loss (HBA2 at half depth): three genes, a silent carrier",
          (r["status"], r["alpha_globin_genes"]) == ("called", "3") and "silent carrier" in r["message"], r)
    check("a count under four says it is unconfirmed; four does not",
          "unconfirmed" in r["message"] and "needs a clinical test" in r["message"]
          and "unconfirmed" not in run(depths(30, 30))["message"]
          and "unconfirmed" not in run(depths(45, 30))["message"], r)
    r = run(depths(15, 15))
    check("both genes at half depth: two, and depth cannot tell --/aa from -a/-a",
          (r["status"], r["alpha_globin_genes"]) == ("called", "2") and "cannot tell these apart" in r["message"], r)
    r = run(depths(45, 30))
    check("HBA2 at 1.5x: five genes, a triplication",
          r["alpha_globin_genes"] == "5" and "triplication" in r["message"], r)
    r = run(depths(22.5, 30))
    check("an estimate of 3.50 is not assessed, with both whole numbers named",
          r["status"] == "not assessed" and r["estimate"] == "3.50" and "between 3 and 4" in r["message"], r)
    r = run(depths(5, 5, 6, 6))
    check("flanks under 8x: not assessed, too low",
          r["status"] == "not assessed" and r["alpha_globin_genes"] == "." and "too low" in r["message"], r)
    r = run(depths(30, 30, 30, 15))
    check("flanks 30x and 15x: not assessed, they disagree",
          r["status"] == "not assessed" and "disagree" in r["message"], r)
    novo = depths(15.5, 25.68, 30.11, 31.25)
    r = run(novo, share=0)
    check("no read of MAPQ 0 over the genes (novoalign's HG002): not assessed, the aligner named",
          r["status"] == "not assessed" and r["alpha_globin_genes"] == "." and r["mapq0_share"] == "0.000"
          and "mapping quality 0" in r["message"] and "novoalign" in r["message"], r)
    r = run(depths(15, 30), share=0)
    check("no read of MAPQ 0 with an estimate of 3.00: not assessed either, never three",
          r["status"] == "not assessed" and r["estimate"] == "3.00" and "mapping quality 0" in r["message"], r)
    r = run(depths(30, 30), share=0.01)
    check("1% of the gene depth at MAPQ 0: not assessed", r["status"] == "not assessed"
          and "1.0%" in r["message"], r)
    r = run(depths(30, 30), share=0.03)
    check("3% at MAPQ 0: called, four", (r["status"], r["alpha_globin_genes"], r["mapq0_share"]) == ("called", "4", "0.030"), r)
    r = run(depths(0, 0), share=0)
    check("no gene at all: still called 0, with no reads to judge the aligner by",
          (r["status"], r["alpha_globin_genes"]) == ("called", "0") and "Bart" in r["message"], r)

    with tempfile.TemporaryDirectory() as t:
        out = subprocess.run([sys.executable, SCRIPT, "bed"], capture_output=True, text=True)
        rows = [line.split("\t") for line in out.stdout.splitlines()]
        check("bed: four regions, the two genes between the flanks",
              out.returncode == 0 and [x[3] for x in rows] == ["flank_left", "HBA2", "HBA1", "flank_right"], out)
        cov, covq, res = os.path.join(t, "cov.tsv"), os.path.join(t, "cov_q1.tsv"), os.path.join(t, "hba.tsv")
        bedcov(cov, depths(15, 30))
        bedcov(covq, mapq1(depths(15, 30)))
        cli = [sys.executable, SCRIPT, "check", "--bedcov", cov, "--bedcov-mapq1", covq, "--out", res]
        out = subprocess.run(cli, capture_output=True, text=True)
        table = dict(line.split("\t", 1) for line in open(res).read().splitlines()[1:]) if os.path.exists(res) else {}
        check("check: exit 0, the table has the status, the count and the message",
              out.returncode == 0 and table.get("status") == "called" and table.get("alpha_globin_genes") == "3"
              and table.get("mapq0_share") == "0.250"
              and list(table) == list(hba.ORDER) and "HBA1/HBA2: three" in out.stdout, (out, table))
        out = subprocess.run(cli[:4] + cli[6:], capture_output=True, text=True)
        check("check: without --bedcov-mapq1 exits 2", out.returncode == 2 and "--bedcov-mapq1" in out.stderr, out)
        for name, path in (("all reads", cov), ("MAPQ 1 or more", covq)):
            bedcov(cov, depths(15, 30))
            bedcov(covq, mapq1(depths(15, 30)))
            with open(path, "a") as f:
                f.write("chr16\t1\t100\tother\t500\n")
            out = subprocess.run(cli, capture_output=True, text=True)
            check(f"check: a table ({name}) with other regions exits 2", out.returncode == 2 and "differ" in out.stderr, out)
            with open(path, "w") as f:
                f.write("chr16\t100000\t140000\n")
            out = subprocess.run(cli, capture_output=True, text=True)
            check(f"check: a short row ({name}) exits 2", out.returncode == 2 and "expected chrom" in out.stderr, out)

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
