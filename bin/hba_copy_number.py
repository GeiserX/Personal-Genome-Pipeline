#!/usr/bin/env python3
"""hba_copy_number.py: HBA1 and HBA2 copy number from short-read depth (step 35).

Parascopy's pinned models have none for the alpha-globin genes, so step 35
counts them from read depth instead. HBA1 and HBA2 are near-identical, so
short reads cannot be placed on one of them; but every read from either
lands on one of the two, so the depth over both gene bodies together, against
single-copy sequence on each side of the cluster, counts the alpha-globin
genes: four on a genome without a deletion (aa/aa), three with one gene lost
on one chromosome (-a3.7 or -a4.2), two with both lost on one chromosome
(--SEA and other large deletions) or one on each (-a/-a).

  hba_copy_number.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  hba_copy_number.py check --bedcov BEDCOV.tsv --out HBA.tsv
      read `samtools bedcov` over those regions (chrom, start, end, name,
      summed depth; all MAPQ, duplicates left out) and write HBA.tsv
      (metric<TAB>value: status, alpha_globin_genes, estimate, the depths and
      the message). Exit 0 whatever the status; 2 on unreadable input.

status is 'called' or 'not assessed'. The estimate is
2 x (HBA2 depth + HBA1 depth) / flank depth. It is called only when it lies
within 0.3 of a whole number, the flanks reach MIN_FLANK_DEPTH and agree
with each other; otherwise the message says why not. Standard library only:
it runs in the Parascopy image and under tests/test_hba_copy_number.py.
"""
import argparse
import sys

# GRCh38, 0-based. The gene bodies are the MANE Select transcripts (HBA2
# NM_000517.6, chr16:172,876-173,710; HBA1 NM_000558.5, chr16:176,680-177,522).
# The flanks are 40 kb of single-copy sequence on each side of the cluster,
# outside the --SEA deletion (about chr16:165,400-184,700) and away from the
# HBZ/HBZP1 and HBM paralogs.
REGIONS = [
    ("chr16", 100000, 140000, "flank_left"),
    ("chr16", 172875, 173710, "HBA2"),
    ("chr16", 176679, 177522, "HBA1"),
    ("chr16", 190000, 230000, "flank_right"),
]
# Below this mean flank depth the ratio is noise: not assessed.
MIN_FLANK_DEPTH = 8.0
# The two flanks must agree within this factor: a deletion that reaches one
# flank, or a dropout, would otherwise skew the ratio.
MAX_FLANK_RATIO = 1.5
# An estimate further than this from a whole number is not called.
MAX_OFF = 0.3

READINGS = {
    0: "no alpha-globin gene: the pattern of Hb Bart's hydrops fetalis (--/--)",
    1: "one alpha-globin gene: the pattern of HbH disease (--/-a)",
    2: ("two alpha-globin genes: alpha-thalassaemia trait, either both genes lost on one chromosome "
        "(--/aa, as --SEA; a child with a partner who carries the same can have Hb Bart's) or one on "
        "each (-a/-a); depth cannot tell these apart"),
    3: "three alpha-globin genes: one gene lost on one chromosome (-a/aa, as -a3.7 or -a4.2), a silent carrier",
    4: "four alpha-globin genes (aa/aa): no deletion seen",
}


# On HG002, whose HiFi reads show four genes (Paraphase), this check on GIAB's
# Illumina reads at about 30x counts three: HBA2 at about half the flank depth
# (the e2e case paralog-checks-and-doc-fixes-1-hba). So a count under four is
# reported with this note, never as a finding on its own.
UNDER_FOUR = ("; unconfirmed: on HG002, which has four genes by long reads, this short-read check counts "
              "three, so a count under four needs a clinical test before it means anything")


def reading(n):
    if n in READINGS:
        return READINGS[n] + (UNDER_FOUR if n < 4 else "")
    return f"{n} alpha-globin genes: more than four, a triplication (as aaa/aa)"


def read_bedcov(path):
    """{'flank_left': mean, 'HBA2': mean, ...} from `samtools bedcov` output
    over the regions of REGIONS, which must all be there and nothing else."""
    means, seen = {}, set()
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                raise ValueError(f"{path}: expected chrom, start, end, name, depth sum; got {line.strip()!r}")
            chrom, start, end, name = c[0], int(c[1]), int(c[2]), c[3]
            seen.add((chrom, start, end, name))
            means[name] = float(c[4]) / (end - start)
    # The regions must be the ones `bed` prints, so a drift in a copy of
    # them cannot compare other stretches unnoticed.
    if seen != set(REGIONS):
        raise ValueError(f"{path}: its regions differ from the ones `hba_copy_number.py bed` prints")
    return means


def assess(m):
    """The result dict for the mean depths of read_bedcov."""
    left, right = m["flank_left"], m["flank_right"]
    flank = (left + right) / 2
    out = {"hba2_depth": round(m["HBA2"], 2), "hba1_depth": round(m["HBA1"], 2),
           "flank_left_depth": round(left, 2), "flank_right_depth": round(right, 2),
           "estimate": ".", "alpha_globin_genes": ".", "status": "not assessed"}
    if min(left, right) < MIN_FLANK_DEPTH:
        out["message"] = (f"not assessed: flank depth under {MIN_FLANK_DEPTH:g}x "
                          f"({left:.1f}x and {right:.1f}x), too low to compare")
        return out
    if max(left, right) / min(left, right) > MAX_FLANK_RATIO:
        out["message"] = (f"not assessed: the two flanks disagree ({left:.1f}x and {right:.1f}x); "
                          "a deletion larger than --SEA, or uneven coverage, would skew the count")
        return out
    est = 2 * (m["HBA2"] + m["HBA1"]) / flank
    n = round(est)
    out["estimate"] = f"{est:.2f}"
    if abs(est - n) > MAX_OFF:
        out["message"] = (f"not assessed: the estimate {est:.2f} lies between {int(est)} and {int(est) + 1} "
                          "alpha-globin genes")
        return out
    out.update(status="called", alpha_globin_genes=str(n), message=reading(n))
    return out


ORDER = ("status", "alpha_globin_genes", "estimate", "hba2_depth", "hba1_depth",
         "flank_left_depth", "flank_right_depth", "message")


def check(args):
    try:
        res = assess(read_bedcov(args.bedcov))
    except (OSError, ValueError, KeyError, ZeroDivisionError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    try:
        with open(args.out, "w") as f:
            f.write("metric\tvalue\n")
            for k in ORDER:
                f.write(f"{k}\t{res[k]}\n")
    except OSError as e:
        print(f"ERROR: cannot write {args.out}: {e}", file=sys.stderr)
        return 2
    print(f"HBA1/HBA2: {res['message']} (estimate {res['estimate']}; HBA2 {res['hba2_depth']}x, "
          f"HBA1 {res['hba1_depth']}x, flanks {res['flank_left_depth']}x and {res['flank_right_depth']}x)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hba_copy_number.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the regions")
    c = sub.add_parser("check", help="judge samtools bedcov over the regions")
    c.add_argument("--bedcov", required=True)
    c.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "bed":
        for r in REGIONS:
            print("\t".join(map(str, r)))
        return 0
    return check(args)


if __name__ == "__main__":
    sys.exit(main())
