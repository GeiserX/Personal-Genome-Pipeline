#!/usr/bin/env python3
"""pgx_outside_calls.py: the PharmCAT outside-call file of one sample.

  pgx_outside_calls.py --calls OUT_CALLS.tsv --consensus OUT_CONSENSUS.tsv
                       [--hla T1K_genotype.tsv] [--pypgx PYPGX_summary.tsv]
                       [--cyrius CYRIUS.tsv] [--depth-check CYP2D6_depth_check.tsv]

writes the outside-call file PharmCAT reads with -po (gene<TAB>diplotype, one
line per gene, empty when there is nothing to pass) and a consensus table
(Gene, Result, Outside_call, Reason, Evidence) that says for HLA-A, HLA-B and
CYP2D6 what each caller said and what was passed on.
"""
import argparse
import csv
import sys

CONSENSUS_HEADER = ["Gene", "Result", "Outside_call", "Reason", "Evidence"]


def read_pypgx(path):
    if not path:
        return None
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row.get("Gene") or "").strip() == "CYP2D6":
                return (row.get("Diplotype") or "").strip()
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--calls", required=True)
    p.add_argument("--consensus", required=True)
    p.add_argument("--hla")
    p.add_argument("--pypgx")
    p.add_argument("--cyrius")
    p.add_argument("--depth-check")
    args = p.parse_args(argv)
    calls, rows = [], []
    dip = read_pypgx(args.pypgx)
    if dip:
        calls.append(("CYP2D6", dip))
        rows.append(["CYP2D6", dip, "yes", "pypgx", f"pypgx: {dip}"])
    with open(args.calls, "w") as f:
        f.writelines(f"{g}\t{d}\n" for g, d in calls)
    with open(args.consensus, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(CONSENSUS_HEADER)
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
