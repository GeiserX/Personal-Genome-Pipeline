#!/usr/bin/env python3
"""cyp2d6_depth_check.py: is the read depth at CYP2D6 fit for a copy-number call?

pypgx (step 32) and Cyrius (step 21) call CYP2D6 copy number from read depth.
Before either calls, mosdepth measures the mean depth over CYP2D6 and two
50 kb flanks outside the CYP2D6-CYP2D8 cluster, once for all reads and once
for reads with MAPQ >= 1 (callers ignore MAPQ 0). This script reads the two
mosdepth region files and says whether the depth can be trusted.

  cyp2d6_depth_check.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  cyp2d6_depth_check.py check --all Q0.regions.bed[.gz] --mapq1 Q1.regions.bed[.gz] --out CHECK.tsv
      write CHECK.tsv (metric<TAB>value: status, message and the depths) and
      print the message. Exit 0 whatever the status; 2 on unreadable input.

status is 'ok' or 'unreliable'. Standard library only: it runs in the plain
python image, in the pypgx image and under tests/test_cyp2d6_depth.py.
"""
import argparse
import gzip
import sys

# GRCh38, 0-based. CYP2D6 as Cyrius 1.1.1 defines it (data/CYP2D6_region_38.bed,
# with REP6); the flanks are the ones scripts/ci/alt-depth-ab.sh measures.
REGIONS = [
    ("chr22", 42050000, 42100000, "flank"),
    ("chr22", 42123192, 42132032, "CYP2D6"),
    ("chr22", 42200000, 42250000, "flank"),
]
CHROM = "chr22"

UNRELIABLE = ("CYP2D6 copy number unreliable: reads are multi-mapped "
              "(was this BAM aligned to a reference with ALT contigs?)")
NO_DEPTH = ("CYP2D6 copy number unreliable: no reads in the flanks of CYP2D6 "
            "(chr22:42.05-42.25 Mb), so its depth cannot be compared")
MAX_RATIO = 0.6


def read_regions(path):
    """{'gene': mean, 'flank': mean} from a mosdepth regions file (chrom,
    start, end, name, mean), length-weighted over the rows of each name."""
    total, length = {}, {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}: expected chrom, start, end, name, mean; got {line.strip()!r}")
            name = "gene" if fields[3] == "CYP2D6" else "flank" if "flank" in fields[3].lower() else None
            if name is None:
                continue
            n = int(fields[2]) - int(fields[1])
            total[name] = total.get(name, 0.0) + float(fields[4]) * n
            length[name] = length.get(name, 0) + n
    missing = {"gene", "flank"} - set(length)
    if missing:
        raise ValueError(f"{path}: no {' or '.join(sorted(missing))} region (name CYP2D6, or one containing 'flank')")
    return {k: total[k] / length[k] for k in length}


def assess(gene_all, gene_q1, flank_all, flank_q1):
    """(status, message) for the four mean depths."""
    if flank_all <= 0:
        return "unreliable", NO_DEPTH
    ratio_all = gene_all / flank_all
    ratio_q1 = gene_q1 / flank_q1 if flank_q1 > 0 else 0.0
    if ratio_q1 < MAX_RATIO and ratio_all >= MAX_RATIO:
        return "unreliable", UNRELIABLE
    return "ok", "CYP2D6 depth is fit for a copy-number call"


def check(args):
    try:
        d_all = read_regions(args.all)
        d_q1 = read_regions(args.mapq1)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    status, message = assess(d_all["gene"], d_q1["gene"], d_all["flank"], d_q1["flank"])
    rows = [
        ("status", status),
        ("message", message),
        ("gene_depth_all", f"{d_all['gene']:.2f}"),
        ("gene_depth_mapq1", f"{d_q1['gene']:.2f}"),
        ("flank_depth_all", f"{d_all['flank']:.2f}"),
        ("flank_depth_mapq1", f"{d_q1['flank']:.2f}"),
    ]
    with open(args.out, "w") as f:
        f.write("metric\tvalue\n")
        f.writelines(f"{k}\t{v}\n" for k, v in rows)
    print(f"CYP2D6 depth check: {status}: {message}")
    print("  mean depth, all reads: CYP2D6 {:.2f}, flanks {:.2f}; MAPQ >= 1: CYP2D6 {:.2f}, flanks {:.2f}".format(
        d_all["gene"], d_all["flank"], d_q1["gene"], d_q1["flank"]))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the CYP2D6 and flank regions as BED")
    c = sub.add_parser("check", help="assess two mosdepth region files")
    c.add_argument("--all", required=True, help="mosdepth regions file, all reads (-Q 0)")
    c.add_argument("--mapq1", required=True, help="mosdepth regions file, MAPQ >= 1 (-Q 1)")
    c.add_argument("--out", required=True, help="where to write the check (TSV)")
    args = p.parse_args(argv)
    if args.cmd == "bed":
        for chrom, start, end, name in REGIONS:
            print(f"{chrom}\t{start}\t{end}\t{name}")
        return 0
    return check(args)


if __name__ == "__main__":
    sys.exit(main())
