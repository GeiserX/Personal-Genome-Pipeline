#!/usr/bin/env python3
"""bin/mosaic_sex_loss.py on synthetic mosdepth windows and VCFs.

Checks:
  1. a man with one Y in every cell reads a depth ratio of 0.5 and a loss
     near zero, with zero inside its 95% range, not detected; the chrY
     windows outside the single-copy part (ampliconic, centromere, gaps) do
     not move it, as they would move a whole-chromosome mean;
  2. 30% of the Y reads gone reads a loss of about 0.30, detected; an extra
     Y (47,XYY) reads as a gain, not a loss;
  3. a woman: no loss reads near zero from depth and from allele fractions;
     one X lost in 30% of the cells (chrX depth at 0.85, heterozygous chrX
     sites at 0.5 +- 0.088) reads about 0.30 from both, detected; shifted
     allele fractions with no depth change are detected and called
     copy-neutral;
  4. a noisy library (allele fractions spread more than binomial on every
     chromosome) is not read as a shift, because the spread comes from
     chr1-22;
  5. nothing to read: no sex, no mosdepth windows, a sliced BAM, a VCF that
     is not one: not_estimated with the reason, and the tool never fails;
  6. the sex from the reads wins over the declared one; contigs named
     without 'chr' work;
  7. the command adds its rows to the table once: a rerun replaces them and
     keeps the other rows in their order.

Run: python3 tests/test_mosaic_sex_loss.py
"""
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(REPO, "bin", "mosaic_sex_loss.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import mosaic_sex_loss as m  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def windows(f, chrom, start, end, depth, rnd, prefix="chr"):
    for s in range(start - start % 500, end, 500):
        d = max(0.0, rnd.gauss(depth, depth * 0.1)) if depth else 0.0
        f.write(f"{prefix}{chrom}\t{s}\t{s + 500}\t{d:.2f}\n")


def regions(path, sex, x_depth=None, y_depth=None, seed=1, prefix="chr", slice_only=False):
    """mosdepth --by 500 windows of a 30x genome."""
    rnd = random.Random(seed)
    with gzip.open(path, "wt") as f:
        windows(f, "1", 1_000_000, 4_000_000, 30, rnd, prefix)
        windows(f, "1", 4_000_000, 4_100_000, 0, rnd, prefix)          # a gap
        windows(f, "2", 5_000_000, 6_000_000, 30, rnd, prefix)
        windows(f, "2", 6_000_000, 6_050_000, 400, rnd, prefix)         # a satellite
        if slice_only:
            windows(f, "Y", 2_790_000, 2_800_000, 15, rnd, prefix)
            return path
        xd = x_depth if x_depth is not None else (15 if sex == "male" else 30)
        windows(f, "X", 1_000_000, 2_000_000, 30, rnd, prefix)          # PAR1
        windows(f, "X", 73_200_000, 74_500_000, xd, rnd, prefix)
        windows(f, "X", 89_500_000, 90_500_000, 30 if sex == "male" else 0, rnd, prefix)  # XTR
        if sex == "male":
            yd = y_depth if y_depth is not None else 15
            windows(f, "Y", 0, 2_781_479, 0, rnd, prefix)               # PAR, masked
            windows(f, "Y", 2_781_479, 3_050_044, yd, rnd, prefix)      # single copy
            windows(f, "Y", 3_050_044, 6_235_111, 15, rnd, prefix)      # XTR
            windows(f, "Y", 6_748_713, 7_574_481, yd, rnd, prefix)      # single copy
            windows(f, "Y", 7_574_481, 9_046_914, 60, rnd, prefix)      # ampliconic
            windows(f, "Y", 10_316_944, 10_544_039, 900, rnd, prefix)   # centromere
            windows(f, "Y", 26_673_214, 56_673_214, 0, rnd, prefix)     # heterochromatin, N
    return path


def binom(rnd, n, p):
    return sum(rnd.random() < p for _ in range(n))


def vcf(path, d=0.0, spread=0.0, seed=2, n_auto=3000, n_x=600):
    """Heterozygous SNVs: chr1-22 at 0.5, chrX at 0.5 +- d; spread adds a
    per-site spread of the true fraction (a noisy library), on every
    chromosome."""
    rnd = random.Random(seed)
    lines = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"]

    def site(chrom, pos, p):
        n = rnd.randint(20, 40)
        p = min(0.95, max(0.05, p + rnd.gauss(0, spread))) if spread else p
        a = binom(rnd, n, p)
        lines.append(f"{chrom}\t{pos}\t.\tA\tG\t40\tPASS\t.\tGT:AD:DP\t0/1:{n - a},{a}:{n}\n")
    for i in range(n_auto):
        site("chr1", 1_000_000 + 1000 * i, 0.5)
    for i in range(n_x):
        site("chrX", 73_200_000 + 2000 * i, 0.5 + (d if rnd.random() < 0.5 else -d))
    # Records the estimate must skip: hom-ALT, not PASS, an indel, PAR1, depth below 10.
    lines += ["chrX\t73100000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t1/1:0,30\n",
              "chrX\t73100100\t.\tA\tG\t40\tRefCall\t.\tGT:AD\t0/1:5,25\n",
              "chrX\t73100200\t.\tAT\tA\t40\tPASS\t.\tGT:AD\t0/1:5,25\n",
              "chrX\t1500000\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:2,28\n",
              "chrX\t73100300\t.\tA\tG\t40\tPASS\t.\tGT:AD\t0/1:1,5\n"]
    with gzip.open(path, "wt") as f:
        f.writelines(lines)
    return path


def run(table, regions_path=None, vcf_path=None):
    return dict(m.estimate(table, regions_path, vcf_path))


def num(row, k):
    try:
        return float(row[k])
    except (KeyError, ValueError):
        return None


def near(row, k, want, tol):
    v = num(row, k)
    return v is not None and abs(v - want) <= tol


MALE = [("inferred_sex", "male"), ("declared_sex", "male")]
FEMALE = [("inferred_sex", "female"), ("declared_sex", "female")]


def main():
    w = tempfile.mkdtemp()
    try:
        # 1. a man with one Y in every cell
        r = run(MALE, regions(f"{w}/m.bed.gz", "male"))
        print(r)
        check("male: the normalised chrY ratio reads 0.5", near(r, "mosaic_depth_ratio", 0.5, 0.02), r)
        check("male: the loss is near zero", near(r, "mosaic_loss_depth", 0.0, 0.03), r)
        check("male: zero lies inside the 95% range",
              (num(r, "mosaic_loss_depth_low") or 1) <= 0 <= (num(r, "mosaic_loss_depth_high") or -1), r)
        check("male: not detected, the reason names the 10% floor",
              r["mosaic_call"] == "not_detected" and "below that nothing can be told at 30x" in r["mosaic_reason"], r)
        check("male: the rows say chrY and the floor", r["mosaic_chrom"] == "chrY" and r["mosaic_floor"] == "0.100"
              and r["mosaic_depth_ratio_expected"] == "0.5" and int(r["mosaic_depth_blocks"]) > 100, r)
        check("male: no allele-fraction estimate", r["mosaic_loss_af"] == "" and r["mosaic_x_het_sites"] == "", r)

        # 2. 30% of the Y reads gone; 47,XYY
        r = run(MALE, regions(f"{w}/m30.bed.gz", "male", y_depth=10.5))
        print(r)
        check("male, 30% of Y gone: the loss reads about 0.30", near(r, "mosaic_loss_depth", 0.30, 0.03), r)
        check("male, 30% of Y gone: detected, the low end above 0.10",
              r["mosaic_call"] == "detected" and num(r, "mosaic_loss_depth_low") > 0.10, r)
        r = run(MALE, regions(f"{w}/xyy.bed.gz", "male", y_depth=30))
        check("47,XYY: a gain, not a loss", r["mosaic_call"] == "gain" and near(r, "mosaic_loss_depth", -1.0, 0.05)
              and "47,XYY" in r["mosaic_reason"], r)

        # 3. a woman
        reg_f = regions(f"{w}/f.bed.gz", "female")
        r = run(FEMALE, reg_f, vcf(f"{w}/f.vcf.gz"))
        print(r)
        check("female, no loss: chrX ratio 1.0, the depth loss near zero",
              near(r, "mosaic_depth_ratio", 1.0, 0.03) and near(r, "mosaic_loss_depth", 0.0, 0.05), r)
        check("female, no loss: the allele-fraction loss is below the floor, its low end 0",
              num(r, "mosaic_loss_af") < 0.10 and num(r, "mosaic_loss_af_low") == 0.0, r)
        check("female, no loss: the 600 chrX sites count, the skipped records do not",
              r["mosaic_x_het_sites"] == "600", r)
        check("female, no loss: not detected", r["mosaic_call"] == "not_detected", r)

        d30 = 0.3 / (2 * (2 - 0.3))
        r = run(FEMALE, regions(f"{w}/f30.bed.gz", "female", x_depth=30 * 0.85), vcf(f"{w}/f30.vcf.gz", d=d30))
        print(r)
        check("female, one X lost in 30%: the depth loss reads about 0.30", near(r, "mosaic_loss_depth", 0.30, 0.04), r)
        check("female, one X lost in 30%: the allele-fraction loss reads about 0.30",
              near(r, "mosaic_loss_af", 0.30, 0.08), r)
        check("female, one X lost in 30%: detected", r["mosaic_call"] == "detected"
              and "lost one X" in r["mosaic_reason"], r)
        r = run(FEMALE, reg_f, vcf(f"{w}/fcn.vcf.gz", d=d30))
        check("female, shifted allele fractions, depth unchanged: detected, called copy-neutral",
              r["mosaic_call"] == "detected" and "copy-neutral" in r["mosaic_reason"]
              and near(r, "mosaic_loss_depth", 0.0, 0.05), r)

        # 4. a noisy library
        r = run(FEMALE, reg_f, vcf(f"{w}/noisy.vcf.gz", spread=0.06))
        print(r)
        check("noisy library, no loss: not read as a shift (spread from chr1-22)",
              r["mosaic_call"] == "not_detected" and num(r, "mosaic_loss_af") < 0.10, r)
        an, au, xn, xu = m.het_sites(f"{w}/noisy.vcf.gz")
        naive = (sum(xu) / len(xu) - 1) / (4 * (sum(xn) / len(xn) - 1))
        check("control: the same sites read against binomial spread alone would show a loss above 0.10",
              naive > 0 and 4 * naive ** 0.5 / (1 + 2 * naive ** 0.5) > 0.10, naive)

        # 5. nothing to read
        r = run([("inferred_sex", "unknown"), ("declared_sex", "")], reg_f)
        check("no sex: not_estimated", r["mosaic_call"] == "not_estimated" and "no sex" in r["mosaic_reason"], r)
        r = run(MALE)
        check("no mosdepth windows: not_estimated, says step 16b",
              r["mosaic_call"] == "not_estimated" and "step 16b" in r["mosaic_reason"], r)
        r = run(MALE, regions(f"{w}/slice.bed.gz", "male", slice_only=True))
        check("a sliced BAM: not_estimated, too few blocks",
              r["mosaic_call"] == "not_estimated" and "too few 10 kb blocks" in r["mosaic_reason"], r)
        with open(f"{w}/placeholder.vcf.gz", "w") as f:
            f.write("placeholder\n")
        r = run(FEMALE, None, f"{w}/placeholder.vcf.gz")
        check("a VCF that is not one: not_estimated with the reason, no exception",
              r["mosaic_call"] == "not_estimated" and "allele fractions not estimated" in r["mosaic_reason"], r)
        r = run(FEMALE, reg_f, None)
        check("female without a VCF: the depth still counts, the reason says no VCF",
              r["mosaic_call"] == "not_detected" and "no VCF" in r["mosaic_reason"], r)

        # 6. the sex from the reads; contigs without chr
        r = run([("inferred_sex", "male"), ("declared_sex", "female")], regions(f"{w}/m2.bed.gz", "male"))
        check("the sex from the reads wins", r["mosaic_sex"] == "male" and r["mosaic_sex_from"] == "reads", r)
        r = run([("inferred_sex", "unknown"), ("declared_sex", "male")], regions(f"{w}/m2.bed.gz", "male"))
        check("unknown from the reads: the declared sex", r["mosaic_sex"] == "male"
              and r["mosaic_sex_from"] == "declared", r)
        r = run(MALE, regions(f"{w}/nochr.bed.gz", "male", y_depth=10.5, prefix=""))
        check("contigs named 1, X, Y: the same estimate", near(r, "mosaic_loss_depth", 0.30, 0.03), r)

        # 7. the command on a table
        t = f"{w}/S_sample_qc.tsv"
        with open(t, "w") as f:
            f.write("key\tvalue\nsample\tS\ninferred_sex\tmale\ndeclared_sex\tmale\nsame_person_as\t\n")
        args = [sys.executable, TOOL, "--table", t, "--regions", f"{w}/m30.bed.gz"]
        p1 = subprocess.run(args, capture_output=True, text=True)
        p2 = subprocess.run(args, capture_output=True, text=True)
        keys = [l.split("\t", 1)[0] for l in open(t).read().splitlines()]
        check("the command exits 0 twice", p1.returncode == 0 and p2.returncode == 0, p1.stderr + p2.stderr)
        check("the table keeps its rows first, then each mosaic row once",
              keys == ["key", "sample", "inferred_sex", "declared_sex", "same_person_as"] + m.KEYS, keys)
        check("the command prints the call", "mosaic_call: detected" in p2.stdout, p2.stdout)
    finally:
        shutil.rmtree(w)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
