#!/usr/bin/env python3
"""cyp2d6_depth_check.py: is the read depth at CYP2D6 fit for a copy-number call?

pypgx (step 32) and Cyrius (step 21) call CYP2D6 copy number from read depth.
Before either calls, mosdepth measures the mean depth over CYP2D6 and two
50 kb flanks outside the CYP2D6-CYP2D8 cluster, once for all reads and once
for reads with MAPQ >= 1 (callers ignore MAPQ 0). This script reads the two
mosdepth region files and says whether the depth can be trusted.

  cyp2d6_depth_check.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  cyp2d6_depth_check.py measure --bam BAM --out SHARES.tsv
      count the depth over Cyrius's CYP2D6 region (with REP6) and over the
      gene body alone (without REP6), for all reads, MAPQ >= 1 and
      MAPQ >= 20, from the BAM itself (pysam: run it in the pypgx image)
  cyp2d6_depth_check.py check --all Q0.regions.bed[.gz] --mapq1 Q1.regions.bed[.gz]
                              [--measure SHARES.tsv] --out CHECK.tsv
      write CHECK.tsv (metric<TAB>value: status, message, the depths and the
      MAPQ shares) and print the message. Exit 0 whatever the status; 2 on
      unreadable input. Without --measure the shares it gives read 'not
      measured'.

status is 'ok' or 'unreliable'; the shares do not change it. Standard
library only, except `measure`, which imports pysam: the rest runs in the
plain python image, in the pypgx image and under tests/test_cyp2d6_depth.py.
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
# The gene body without REP6: CYP2D6's MANE Select transcript
# (ENST00000645361, NM_000106.6), GRCh38 chr22:42,126,499-42,130,810, 0-based
# here. Cyrius's region above adds REP6 and the 3' flank, where much of the
# depth is MAPQ 0 even on a clean BAM; this one says how the gene itself maps.
BODY = (CHROM, 42126498, 42130810)
# The MAPQ floors `measure` counts at: all reads, what the callers keep, and
# a floor that leaves out the reads an aligner placed with little confidence.
MEASURE_MAPQ = (0, 1, 20)
# Reads mosdepth leaves out by default (its -F 1796): unmapped, secondary,
# QC-failed and duplicate. `measure` leaves out the same.
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
NOT_MEASURED = "not measured"

UNRELIABLE = ("CYP2D6 copy number unreliable: reads are multi-mapped "
              "(was this BAM aligned to a reference with ALT contigs?)")
NO_DEPTH = ("CYP2D6 copy number unreliable: no reads in the flanks of CYP2D6 "
            "(chr22:42.05-42.25 Mb), so its depth cannot be compared")
# A depth fit for a copy-number call has at least this share of its reads at
# MAPQ >= 1, at CYP2D6 and in the flanks.
MIN_MAPQ1_SHARE = 0.6
# Below this ratio of CYP2D6 to flank depth (all reads) the gene has too few
# reads for its MAPQ share to mean anything: both copies deleted (*5/*5).
MIN_GENE_RATIO = 0.15


def read_regions(path):
    """{'gene': mean, 'flank': mean} from a mosdepth regions file (chrom,
    start, end, name, mean), length-weighted over the rows of each name."""
    total, length, seen = {}, {}, set()
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
            seen.add((fields[0], int(fields[1]), int(fields[2]), name))
            n = int(fields[2]) - int(fields[1])
            total[name] = total.get(name, 0.0) + float(fields[4]) * n
            length[name] = length.get(name, 0) + n
    missing = {"gene", "flank"} - set(length)
    if missing:
        raise ValueError(f"{path}: no {' or '.join(sorted(missing))} region (name CYP2D6, or one containing 'flank')")
    # The regions must be the ones `bed` prints: the Nextflow module writes
    # its own copy of them, and a drift there would compare other stretches.
    want = {(c, s, e, "gene" if n == "CYP2D6" else "flank") for c, s, e, n in REGIONS}
    if seen != want:
        raise ValueError(f"{path}: its regions differ from the ones `cyp2d6_depth_check.py bed` prints")
    return {k: total[k] / length[k] for k in length}


def assess(gene_all, gene_q1, flank_all, flank_q1):
    """(status, message) for the four mean depths.

    Unreliable when more than 40% of the reads at CYP2D6 or in its flanks have
    MAPQ 0: the aligner placed them on more than one sequence. A deletion
    leaves few reads at CYP2D6, but the ones it leaves map uniquely and the
    flanks are untouched, so it passes.

    Why not compare CYP2D6 with the flanks at MAPQ >= 1 against 0.6, as first
    planned: on the fixture's CYP2D slice (ALT depth A/B run 37177862776) the
    no-ALT BAM has CYP2D6/flanks 0.62 for all reads and 0.57 at MAPQ >= 1, the
    normal value of the locus, which that rule flags; the with-ALT BAM places
    a third of the primary alignments on chr22_KI270879v1_alt, so all-reads
    depth falls too (0.28) and that rule passes it. The MAPQ >= 1 share
    separates the two: 0.91 at CYP2D6 without ALT contigs, 0.005 with them.
    """
    if flank_all <= 0:
        return "unreliable", NO_DEPTH
    ratio_all = gene_all / flank_all
    share_gene = gene_q1 / gene_all if gene_all > 0 else 1.0
    share_flank = flank_q1 / flank_all
    if share_flank < MIN_MAPQ1_SHARE or (ratio_all >= MIN_GENE_RATIO and share_gene < MIN_MAPQ1_SHARE):
        return "unreliable", UNRELIABLE
    return "ok", "CYP2D6 depth is fit for a copy-number call"


def mean_depth(reads, start, end, min_mapq):
    """Mean depth over [start, end) of reads given as (flag, mapq, aligned
    blocks), the blocks 0-based half-open as pysam's get_blocks() gives them.
    Only aligned bases count, as in mosdepth; a read mosdepth skips by
    default (SKIP_FLAGS) or below min_mapq counts nowhere."""
    total = 0
    for flag, mapq, blocks in reads:
        if flag & SKIP_FLAGS or mapq < min_mapq:
            continue
        for s, e in blocks:
            total += max(0, min(e, end) - max(s, start))
    return total / (end - start)


def measure_rows(reads):
    """(region, min_mapq, mean depth) for Cyrius's region and the gene body,
    at each floor of MEASURE_MAPQ."""
    gene = next((s, e) for c, s, e, n in REGIONS if n == "CYP2D6")
    return [(name, q, mean_depth(reads, s, e, q))
            for name, (s, e) in (("gene", gene), ("body", BODY[1:])) for q in MEASURE_MAPQ]


def measure(args):
    try:
        import pysam
    except ImportError:
        print("ERROR: measure needs pysam (run it in the pypgx image)", file=sys.stderr)
        return 2
    gene = next((s, e) for c, s, e, n in REGIONS if n == "CYP2D6")
    try:
        with pysam.AlignmentFile(args.bam) as bam:
            chrom = CHROM if CHROM in bam.references else CHROM[3:] if CHROM[3:] in bam.references else None
            if chrom is None:
                raise ValueError(f"{args.bam} has no {CHROM}")
            reads = [(r.flag, r.mapping_quality, r.get_blocks()) for r in bam.fetch(chrom, gene[0], gene[1])]
        rows = measure_rows(reads)
        with open(args.out, "w") as f:
            f.write("region\tmin_mapq\tmean_depth\n")
            f.writelines(f"{n}\t{q}\t{d:.4f}\n" for n, q, d in rows)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    return 0


def read_measure(path):
    """{(region, min_mapq): mean depth} from `measure`'s file."""
    out = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head != ["region", "min_mapq", "mean_depth"]:
            raise ValueError(f"{path}: not a file of `cyp2d6_depth_check.py measure`")
        for line in f:
            if line.strip():
                n, q, d = line.rstrip("\n").split("\t")
                out[(n, int(q))] = float(d)
    return out


def share(num, den):
    """num/den as two decimals, 'no reads' when den is 0, NOT_MEASURED
    when either is missing."""
    if num is None or den is None:
        return NOT_MEASURED
    return f"{num / den:.2f}" if den > 0 else "no reads"


def check(args):
    try:
        d_all = read_regions(args.all)
        d_q1 = read_regions(args.mapq1)
        m = read_measure(args.measure) if args.measure else {}
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    status, message = assess(d_all["gene"], d_q1["gene"], d_all["flank"], d_q1["flank"])
    body_all = m.get(("body", 0))
    rows = [
        ("status", status),
        ("message", message),
        ("gene_depth_all", f"{d_all['gene']:.2f}"),
        ("gene_depth_mapq1", f"{d_q1['gene']:.2f}"),
        ("flank_depth_all", f"{d_all['flank']:.2f}"),
        ("flank_depth_mapq1", f"{d_q1['flank']:.2f}"),
        # The share the verdict reads (mosdepth, Cyrius's region with REP6),
        # then what `measure` counted: the same region at MAPQ >= 20, and the
        # gene body without REP6. These do not change the status.
        ("gene_mapq1_share", share(d_q1["gene"], d_all["gene"])),
        ("gene_mapq20_share", share(m.get(("gene", 20)), m.get(("gene", 0)))),
        ("body_depth_all", NOT_MEASURED if body_all is None else f"{body_all:.2f}"),
        ("body_mapq1_share", share(m.get(("body", 1)), body_all)),
        ("body_mapq20_share", share(m.get(("body", 20)), body_all)),
    ]
    try:
        with open(args.out, "w") as f:
            f.write("metric\tvalue\n")
            f.writelines(f"{k}\t{v}\n" for k, v in rows)
    except OSError as e:
        print(f"ERROR: cannot write {args.out}: {e}", file=sys.stderr)
        return 2
    print(f"CYP2D6 depth check: {status}: {message}")
    print("  mean depth, all reads: CYP2D6 {:.2f}, flanks {:.2f}; MAPQ >= 1: CYP2D6 {:.2f}, flanks {:.2f}".format(
        d_all["gene"], d_all["flank"], d_q1["gene"], d_q1["flank"]))
    kv = dict(rows)
    print("  share of the depth at MAPQ >= 1 / >= 20: Cyrius's region (with REP6) {} / {}; gene body without REP6 {} / {}".format(
        kv["gene_mapq1_share"], kv["gene_mapq20_share"], kv["body_mapq1_share"], kv["body_mapq20_share"]))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the CYP2D6 and flank regions as BED")
    c = sub.add_parser("check", help="assess two mosdepth region files")
    c.add_argument("--all", required=True, help="mosdepth regions file, all reads (-Q 0)")
    c.add_argument("--mapq1", required=True, help="mosdepth regions file, MAPQ >= 1 (-Q 1)")
    c.add_argument("--measure", help="the file of `measure` (the MAPQ >= 20 and gene-body shares)")
    c.add_argument("--out", required=True, help="where to write the check (TSV)")
    m = sub.add_parser("measure", help="count the MAPQ >= 1 and >= 20 depth of the region and the gene body (pysam)")
    m.add_argument("--bam", required=True)
    m.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.cmd == "bed":
        for chrom, start, end, name in REGIONS:
            print(f"{chrom}\t{start}\t{end}\t{name}")
        return 0
    if args.cmd == "measure":
        return measure(args)
    return check(args)


if __name__ == "__main__":
    sys.exit(main())
