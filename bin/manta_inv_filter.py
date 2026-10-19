#!/usr/bin/env python
"""Keep a whole-span inversion only where Manta saw both of its junctions.

Manta writes each inversion junction as a pair of breakend (BND) records,
flagged INV3 (open 3' of the reported position) or INV5. A balanced
inversion has two junctions, one of each. Manta's convertInversion.py turns
every same-chromosome pair into one SVTYPE=INV record from POS to END, so a
lone junction whose two breakends lie 150 Mb apart becomes a 150 Mb
inversion, and AnnotSV then annotates it over hundreds of genes.

This script reads Manta's own file and the converted one and writes the
converted file again, with three changes:

  - an INV record whose opposite junction is present (an INV5 record for an
    INV3 one, or the reverse, with POS and END each within --pair-bp, or in
    the same Manta EVENT) gets INV_JUNCTIONS=2;
  - an INV record without one gets INV_JUNCTIONS=1. When it spans
    --max-single-bp or more, it is replaced by its two breakend records from
    Manta's file (each with INV_JUNCTIONS=1), so nothing downstream reads it
    as a whole-span inversion;
  - a single-junction INV record shorter than --short-bp also gets
    INVDUP_CANDIDATE: two breakends that close are more likely an inverted
    insertion or duplication than an inversion.

The records come out sorted by the ##contig order of the header, then POS.
One summary line goes to stderr.

It runs inside the Manta image (Python 2.7) in step 04 and in the Nextflow
MANTA process, and under Python 3 in the unit tests, so it uses only what
both have.

Usage: manta_inv_filter.py RAW_VCF CONVERTED_VCF [--pair-bp N]
       [--max-single-bp N] [--short-bp N] > OUT.vcf
"""
from __future__ import print_function

import argparse
import gzip
import sys


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        if sys.version_info[0] >= 3:
            return gzip.open(path, "rt")
        return gzip.open(path, "rb")
    return open(path)


def info_dict(info):
    d = {}
    if info in ("", "."):
        return d
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def add_info(fields, item):
    if fields[7] in ("", "."):
        fields[7] = item
    else:
        fields[7] = fields[7] + ";" + item


def read_vcf(path):
    header, records = [], []
    fh = open_text(path)
    try:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header.append(line)
            else:
                records.append(line.split("\t"))
    finally:
        if fh is not sys.stdin:
            fh.close()
    return header, records


def pair_key(record_id):
    """MantaBND:12:0:1:0:0:0:0 -> :12:0:1:0:0:0, the part convertInversion.py keeps."""
    if not record_id.startswith("MantaBND"):
        return None
    rest = record_id[len("MantaBND"):]
    i = rest.rfind(":")
    return rest[:i] if i > 0 else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("raw", help="Manta's own VCF (diploidSV.raw.vcf.gz)")
    ap.add_argument("converted", help="convertInversion.py's output for it")
    ap.add_argument("--pair-bp", type=int, default=1000,
                    help="largest POS and END difference of two reciprocal junctions (default 1000)")
    ap.add_argument("--max-single-bp", type=int, default=1000000,
                    help="a single-junction inversion this long or longer stays a breakend pair (default 1000000)")
    ap.add_argument("--short-bp", type=int, default=1000,
                    help="a single-junction inversion shorter than this gets INVDUP_CANDIDATE (default 1000)")
    args = ap.parse_args(argv)

    _, raw_records = read_vcf(args.raw)
    raw_by_pair = {}
    for f in raw_records:
        k = pair_key(f[2])
        if k is not None:
            raw_by_pair.setdefault(k, []).append(f)

    header, records = read_vcf(args.converted)

    # The inversion records convertInversion.py wrote, by orientation.
    invs = []
    for n, f in enumerate(records):
        d = info_dict(f[7])
        if d.get("SVTYPE") != "INV" or not f[2].startswith("MantaINV"):
            continue
        if "INV3" in d:
            side = 3
        elif "INV5" in d:
            side = 5
        else:
            continue
        try:
            end = int(d["END"])
        except (KeyError, ValueError):
            continue
        invs.append((n, f[0], int(f[1]), end, side, d.get("EVENT")))

    paired = set()
    for a in invs:
        for b in invs:
            if a[4] != 3 or b[4] != 5 or a[1] != b[1]:
                continue
            near = abs(a[2] - b[2]) <= args.pair_bp and abs(a[3] - b[3]) <= args.pair_bp
            same_event = a[5] is not None and a[5] == b[5]
            if near or same_event:
                paired.add(a[0])
                paired.add(b[0])

    n_paired = n_single = n_short = n_long = 0
    out = []
    replaced = set()
    for n, chrom, pos, end, side, _ in invs:
        f = records[n]
        if n in paired:
            add_info(f, "INV_JUNCTIONS=2")
            n_paired += 1
            continue
        span = end - pos
        mates = raw_by_pair.get(f[2][len("MantaINV"):], [])
        if span >= args.max_single_bp and len(mates) == 2:
            for m in mates:
                m = list(m)
                add_info(m, "INV_JUNCTIONS=1")
                out.append(m)
            replaced.add(n)
            n_long += 1
            continue
        add_info(f, "INV_JUNCTIONS=1")
        n_single += 1
        if span < args.short_bp:
            add_info(f, "INVDUP_CANDIDATE")
            n_short += 1

    out.extend(f for n, f in enumerate(records) if n not in replaced)

    contigs = []
    for h in header:
        if h.startswith("##contig=<ID="):
            contigs.append(h[len("##contig=<ID="):].split(",")[0].rstrip(">"))
    order = dict((c, i) for i, c in enumerate(contigs))
    for f in out:
        if f[0] not in order:
            order[f[0]] = len(order)
    out.sort(key=lambda f: (order[f[0]], int(f[1])))

    new_info = [
        '##INFO=<ID=INV_JUNCTIONS,Number=1,Type=Integer,Description="Junctions Manta found for this inversion: '
        '2 when the reciprocal INV3 and INV5 junctions are both present (POS and END within {0} bp, or one '
        'Manta EVENT), 1 when only one is; a single junction of {1} bp or more stays a breakend pair">'
        .format(args.pair_bp, args.max_single_bp),
        '##INFO=<ID=INVDUP_CANDIDATE,Number=0,Type=Flag,Description="Single-junction inversion shorter than '
        '{0} bp: more likely an inverted insertion or duplication than an inversion">'.format(args.short_bp),
    ]
    w = sys.stdout
    for h in header:
        if h.startswith("#CHROM"):
            for line in new_info:
                w.write(line + "\n")
        w.write(h + "\n")
    for f in out:
        w.write("\t".join(f) + "\n")

    sys.stderr.write(
        "Inversion check: {0} INV records with both junctions, {1} single-junction INV records kept "
        "({2} shorter than {3} bp marked INVDUP_CANDIDATE), {4} single-junction inversions of {5} bp "
        "or more left as their breakend pairs\n".format(
            n_paired, n_single, n_short, args.short_bp, n_long, args.max_single_bp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
