#!/usr/bin/env python3
"""locityper_report.py: step 53's table, from what Locityper genotype wrote.

Usage:
  locityper_report.py --gt-dir DIR --loci BED --exit-code N --out TSV

DIR is the output folder of `locityper genotype` (DIR/loci/<locus>/res.json.gz,
with a `success` file beside each finished locus). BED is the locus list the
database was built from (assets/locityper_loci.bed: chrom, start, end, name;
'#' lines are comments). One row per BED locus, in the BED's order:

  locus, region          the name and the GRCh38 region (1-based)
  status                 called, or not called (with the reason)
  genotype               the two pangenome haplotypes Locityper picked,
  haplotype_1/2          as <sample>.<1|2> of the panel (the HPRC person and
                         which of their two copies), "." when not called
  quality                Locityper's genotype quality, Phred: the chance that
                         another genotype (often a very similar one) is right
  total_reads            reads Locityper used at the locus
  unexplained_reads      reads that fit some haplotype but not the two picked
  weighted_distance      the distance from the call to the other likely
                         genotypes, weighted by their probability (edit
                         distance in bp with a database augmented by
                         `locityper augment`, as setup.sh builds it)
  warnings               Locityper's own warnings, "." when none
  filter                 pass or fail: the filter the Locityper authors used
                         for their HLA and KIR benchmark (Prodanov et al.
                         2025, Zenodo 15879519 README): no warning, a weighted
                         distance of 30 or less, and fewer than 1000 or fewer
                         than 20% unexplained reads. "." when not called
  reason                 why a locus was not called, or why the filter failed

A locus Locityper did not finish is "not called", never a genotype. Exit code
N is Locityper's: a run that stopped calls nothing it did not finish. Exit 1
only when this script cannot read its inputs.
"""
import argparse
import csv
import gzip
import json
import math
import os
import sys
import zlib

COLUMNS = ["locus", "region", "status", "genotype", "haplotype_1", "haplotype_2", "quality", "total_reads",
           "unexplained_reads", "weighted_distance", "warnings", "filter", "reason"]

# The authors' filter of their HLA/KIR benchmark (Zenodo 15879519 README,
# "Benchmarking HLA/KIR calling"): qual 1 when there are no warnings, the
# weighted distance is at most 30, and there are fewer than 1000 or fewer
# than 20% unexplained reads.
MAX_WEIGHT_DIST = 30.0
MAX_UNEXPLAINED = 1000
MAX_UNEXPLAINED_SHARE = 0.2


def read_loci(path):
    loci = []
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or not c[1].isdigit() or not c[2].isdigit():
                raise ValueError(f"{path}:{n}: want chrom, start, end, name (tab separated)")
            loci.append((c[3], f"{c[0]}:{int(c[1]) + 1}-{c[2]}"))
    if not loci:
        raise ValueError(f"{path} lists no locus")
    return loci


def fmt(x, digits):
    if x is None:
        return "."
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "."
    if math.isnan(v):
        return "."
    return f"{v:.{digits}f}"


def warnings_of(res):
    w = res.get("warnings")
    if w in (None, "", "*", []):
        return []
    if isinstance(w, str):
        return [x for x in w.replace(";", ",").split(",") if x.strip()]
    return [str(x) for x in w]


def judge(res, warns):
    """pass or fail by the authors' filter, with the reason it failed."""
    why = []
    if warns:
        why.append("Locityper warned: " + ", ".join(warns))
    wd = res.get("weight_dist")
    try:
        wd = float(wd)
    except (TypeError, ValueError):
        wd = float("nan")
    if math.isnan(wd):
        why.append("no weighted distance (the database has no haplotype alignments: setup.sh runs locityper augment)")
    elif wd > MAX_WEIGHT_DIST:
        why.append(f"weighted distance {wd:.2f} is above {MAX_WEIGHT_DIST:g}")
    total = int(res.get("total_reads") or 0)
    unexpl = int(res.get("unexpl_reads") or 0)
    if not (unexpl < MAX_UNEXPLAINED or (total > 0 and unexpl < MAX_UNEXPLAINED_SHARE * total)):
        why.append(f"{unexpl} of {total} reads unexplained (1000 or more, and 20% or more)")
    return ("fail" if why else "pass"), "; ".join(why)


def row_for(name, region, gt_dir, exit_code):
    d = os.path.join(gt_dir, "loci", name)
    res_path = os.path.join(d, "res.json.gz")
    blank = {c: "." for c in COLUMNS}
    blank.update(locus=name, region=region, status="not called")
    if not os.path.isfile(os.path.join(d, "success")) or not os.path.isfile(res_path):
        blank["reason"] = ("Locityper stopped before this locus (exit code %d)" % exit_code if exit_code
                           else "Locityper wrote no result for this locus")
        return blank
    try:
        with gzip.open(res_path, "rt") as f:
            res = json.load(f)
    except (OSError, EOFError, ValueError, zlib.error) as e:  # EOFError: a truncated gzip; zlib.error: a corrupt body
        blank["reason"] = f"res.json.gz cannot be read: {e or type(e).__name__}"
        return blank
    gt = (res.get("genotype") or "").strip()
    haps = [h for h in gt.split(",") if h]
    if len(haps) != 2:
        blank["reason"] = f"Locityper gave no diploid genotype ({gt or 'empty'})"
        return blank
    warns = warnings_of(res)
    filt, why = judge(res, warns)
    return {
        "locus": name, "region": region, "status": "called", "genotype": ",".join(haps),
        "haplotype_1": haps[0], "haplotype_2": haps[1], "quality": fmt(res.get("quality"), 1),
        "total_reads": str(int(res.get("total_reads") or 0)), "unexplained_reads": str(int(res.get("unexpl_reads") or 0)),
        "weighted_distance": fmt(res.get("weight_dist"), 3), "warnings": ",".join(warns) or ".",
        "filter": filt, "reason": why or ".",
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--gt-dir", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--exit-code", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        loci = read_loci(a.loci)
    except (OSError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    rows = [row_for(n, r, a.gt_dir, a.exit_code) for n, r in loci]
    tmp = a.out + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, a.out)
    called = sum(r["status"] == "called" for r in rows)
    print(f"{called} of {len(rows)} loci called, {sum(r['filter'] == 'pass' for r in rows)} pass the authors' filter")
    return 0


if __name__ == "__main__":
    sys.exit(main())
