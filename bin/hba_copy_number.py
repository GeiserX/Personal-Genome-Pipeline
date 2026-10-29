#!/usr/bin/env python3
"""hba_copy_number.py: HBA1 and HBA2 copy number from short-read depth (step 35).

Parascopy's pinned models have none for the alpha-globin genes, so step 35
counts them from read depth instead. HBA1 and HBA2 are near-identical, so
short reads cannot be placed on one of them; but an aligner that reports
every read puts each read from either on one of the two (with mapping
quality 0 where it fits both equally), so the depth over both gene bodies
together, against
single-copy sequence on each side of the cluster, counts the alpha-globin
genes: four on a genome without a deletion (aa/aa), three with one gene lost
on one chromosome (-a3.7 or -a4.2), two with both lost on one chromosome
(--SEA and other large deletions) or one on each (-a/-a).

The method and the regions are the HBA row of bin/cn_loci.tsv, run by
bin/copy_number.py (steps 35 and the PARASCOPY module call that engine);
this file lays out the step's table and keeps its own command line.

  hba_copy_number.py bed
      print the regions (GRCh38, 0-based BED, name in column 4) to stdout
  hba_copy_number.py check --bedcov BEDCOV.tsv --bedcov-mapq1 BEDCOV_Q1.tsv --out HBA.tsv
      read `samtools bedcov` over those regions (chrom, start, end, name,
      summed depth; all MAPQ, duplicates left out) and `samtools bedcov -Q 1`
      over the same regions (reads of MAPQ 1 or more only), and write HBA.tsv
      (metric<TAB>value: status, alpha_globin_genes, estimate, the depths,
      mapq0_share and the message). Exit 0 whatever the status; 2 on
      unreadable input.

status is 'called' or 'not assessed'. The estimate is
2 x (HBA2 depth + HBA1 depth) / flank depth. It is called only when it lies
within 0.3 of a whole number, the flanks reach MIN_FLANK_DEPTH and agree
with each other, and (unless the estimate rounds to 0) reads of MAPQ 0 make
up at least MIN_MAPQ0_SHARE of the depth over the two genes; otherwise the
message says why not. Standard library only:
it runs in the Parascopy image and under tests/test_hba_copy_number.py.
"""
import argparse
import sys

import copy_number as cn

LOCUS = "HBA"
_ROW = []


def row():
    """The HBA row of the catalog, read once."""
    if not _ROW:
        _ROW.append(cn.load_catalog()[LOCUS])
    return _ROW[0]


def __getattr__(name):
    # The regions and thresholds come from the catalog row, read on first use.
    r = row()
    vals = {"REGIONS": cn.regions(r), "MIN_FLANK_DEPTH": r["_min_control"], "MAX_FLANK_RATIO": r["_max_ratio"],
            "MAX_OFF": r["_call"][1][0], "MIN_MAPQ0_SHARE": r["_mapq0_min"]}
    if name in vals:
        return vals[name]
    raise AttributeError(name)


def read_bedcov(path):
    """{'flank_left': mean, 'HBA2': mean, ...} from `samtools bedcov` output
    over the regions, which must all be there and nothing else."""
    m = cn.read_bedcov(path, row())
    if not m:
        raise ValueError(f"{path}: its regions differ from the ones `hba_copy_number.py bed` prints (empty)")
    return m


def table(b):
    """([(metric, value)], summary line) of the engine's HBA block."""
    r, m = row(), b["depths"]
    if b["reason"] == "no_data":
        out = {k: "." for k in ORDER}
        out.update(status="not assessed", message="not assessed: the BAM has no chr16 reads, so HBA1 and HBA2 have none")
        return [(k, out[k]) for k in ORDER], f"HBA1/HBA2: {out['message']}"
    left, right = m.get("flank_left"), m.get("flank_right")
    share = b["mapq0_share"]
    out = {"hba2_depth": round(m["HBA2"], 2), "hba1_depth": round(m["HBA1"], 2),
           "flank_left_depth": round(left, 2), "flank_right_depth": round(right, 2),
           "mapq0_share": f"{share:.3f}",
           "estimate": ".", "alpha_globin_genes": ".", "status": "not assessed"}
    reason, est = b["reason"], b["estimate"]
    if est is not None:
        out["estimate"] = f"{est:.2f}"
    if reason == "low_control_depth":
        out["message"] = (f"not assessed: flank depth under {r['_min_control']:g}x "
                          f"({left:.1f}x and {right:.1f}x), too low to compare")
    elif reason == "controls_disagree":
        out["message"] = (f"not assessed: the two flanks disagree ({left:.1f}x and {right:.1f}x); "
                          "a deletion larger than --SEA, or uneven coverage, would skew the count")
    elif reason == "low_mapq0_share":
        out["message"] = (f"not assessed: reads of mapping quality 0 make up {100 * share:.1f}% of the depth "
                          "over HBA1 and HBA2, so the aligner left out the reads that fit both genes (as "
                          "novoalign does by default) and the count would come out low; realign with "
                          "minimap2 (step 02) or BWA-MEM")
    elif reason == "between_counts":
        out["message"] = (f"not assessed: the estimate {est:.2f} lies between {int(est)} and {int(est) + 1} "
                          "alpha-globin genes")
    else:
        out.update(status="called", alpha_globin_genes=str(b["copies"]), message=b["reading"])
    line = (f"HBA1/HBA2: {out['message']} (estimate {out['estimate']}; HBA2 {out['hba2_depth']}x, "
            f"HBA1 {out['hba1_depth']}x, flanks {out['flank_left_depth']}x and {out['flank_right_depth']}x; "
            f"MAPQ 0 share over the genes {out['mapq0_share']})")
    return [(k, out[k]) for k in ORDER], line


def assess(m, q1):
    """The result dict for the mean depths of read_bedcov: m over every
    read, q1 over reads of MAPQ 1 or more."""
    return dict(table(cn.assess(row(), m, q1))[0])


ORDER = ("status", "alpha_globin_genes", "estimate", "hba2_depth", "hba1_depth",
         "flank_left_depth", "flank_right_depth", "mapq0_share", "message")


def check(args):
    try:
        rows, line = table(cn.assess(row(), read_bedcov(args.bedcov), read_bedcov(args.bedcov_mapq1)))
    except cn.BrokenControl as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3
    except (OSError, ValueError, KeyError, ZeroDivisionError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    try:
        with open(args.out, "w") as f:
            f.write("metric\tvalue\n")
            for k, v in rows:
                f.write(f"{k}\t{v}\n")
    except OSError as e:
        print(f"ERROR: cannot write {args.out}: {e}", file=sys.stderr)
        return 2
    print(line)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="hba_copy_number.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("bed", help="print the regions")
    c = sub.add_parser("check", help="judge samtools bedcov over the regions")
    c.add_argument("--bedcov", required=True)
    c.add_argument("--bedcov-mapq1", required=True, help="samtools bedcov -Q 1 over the same regions")
    c.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "bed":
        for r in cn.regions(row()):
            print("\t".join(map(str, r)))
        return 0
    return check(args)


if __name__ == "__main__":
    sys.exit(main())
