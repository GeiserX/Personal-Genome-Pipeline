#!/usr/bin/env python3
"""Add the CNVpytor and Delly calls Manta does not have to AnnotSV's input.

AnnotSV classifies only the records of the one VCF it is given, and step 05
gives it Manta's calls. A large deletion that CNVpytor (read depth) or Delly
calls, but Manta calls as something else or not at all, never gets an ACMG
class. This script writes the base VCF (Manta's calls, after sv_svlen_fix.py)
again with, after it, the deletions and duplications of the other callers
that no call already in the file covers:

  - Delly (step 19): records with FILTER PASS and SVTYPE DEL or DUP. A call
    longer than 1,000 bp also needs Delly's read-depth copy number (FORMAT
    RDCN) to agree: below the copies a normal region has for a deletion,
    above them for a duplication (1 on a male's chrX or chrY outside the
    PARs given by --par-bed, 2 elsewhere; with no --sex, either on those
    stretches). A missing or negative RDCN is no evidence and the call is
    set aside; a VCF without the RDCN field is read as before. Short-read
    Delly can mark PASS a long SV whose read depth shows no change, and
    one such call reaches the SV report as a class 4 or 5 gene loss;
  - CNVpytor (step 18): DEL and DUP records with e-val (pytorP1) below
    --max-eval (default 0.01) and a fraction of reads with mapping quality 0
    (pytorQ0) below --max-q0 (default 0.5). Copy-neutral LOH is left out.

A call is covered when a record of the same type (DEL, or DUP and DUP:TANDEM)
on the same chromosome overlaps it reciprocally by at least --min-overlap
(default 0.5) of the longer of the two. Delly's calls are added before
CNVpytor's, so a deletion both of them call and Manta misses is added once.

Every record says where it came from in INFO SVCALLER (the base's records
get --base-name, default manta). An added record keeps the caller's ID, END,
CIPOS, CIEND and IMPRECISE, and its genotype as FORMAT GT in the first sample
column. It drops SVLEN: CNVpytor writes it positive for a deletion too, and
AnnotSV, which keeps a given SVLEN, derives it from END with the right sign. The output is sorted by the base header's ##contig order, then
POS. One summary line, with the count of Delly calls read depth set
aside, goes to stderr and into the header as ##svcallers_added.

Usage: cnvpytor_to_vcf.py --base BASE.vcf[.gz]|- [--base-name NAME]
       [--delly DELLY.vcf.gz] [--sex male|female] [--par-bed PAR.bed]
       [--cnvpytor CNVPYTOR.vcf.gz]
       [--max-eval E] [--max-q0 Q] [--min-overlap F] > OUT.vcf
"""
import argparse
import bisect
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_report  # noqa: E402  the PAR intervals, as step 38 reads them

STANDARD_INFO = {
    "END": '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    "SVTYPE": '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    "SVLEN": '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    "IMPRECISE": '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
    "CIPOS": '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">',
    "CIEND": '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">',
    "SVCALLER": '##INFO=<ID=SVCALLER,Number=1,Type=String,Description="The SV caller of this record: the input '
                'of step 05 (manta), or delly or cnvpytor for a call no record of that input covers">',
}
KEPT = ("END", "SVTYPE", "CIPOS", "CIEND")
# A Delly DEL or DUP longer than this needs read depth to agree: its FORMAT
# RDCN (Delly's read-depth copy number) must be below (DEL) or above (DUP) the
# copies a normal region has. 1,000 bp is Delly's own line between small and
# large SVs in its read counts (src/delly.h, indelsize). Short-read Delly can
# mark PASS a long SV whose RDCN shows no change; a short one rests on split
# reads alone. bin/couple_report.py (step 39) uses the same rule.
DEPTH_MIN_SIZE = 1000


def delly_depth_contradicts(kind, chrom, pos, end, fmt, vals, sex, in_par):
    """True when a Delly DEL or DUP longer than DEPTH_MIN_SIZE has a FORMAT
    RDCN that does not show the change: not below the normal copies for a
    DEL, not above them for a DUP. The normal is 1 on a male's chrX or chrY
    outside the PARs and 2 elsewhere; with no sex, either 1 or 2 there, so
    only an RDCN at or past both is held against the call. A missing or
    negative RDCN (no reads in the flanks) is no evidence of a change. A
    record without RDCN (not Delly 2.x) is never contradicted."""
    if end - pos <= DEPTH_MIN_SIZE or "RDCN" not in fmt:
        return False
    i = fmt.index("RDCN")
    cn = vals[i] if i < len(vals) else "."
    if not cn.isdigit():
        return True
    sex_chrom = chrom in ("chrX", "chrY") and not in_par
    if sex_chrom and sex == "male":
        low = high = 1
    elif sex_chrom and sex != "female":
        low, high = 1, 2
    else:
        low = high = 2
    if kind == "DEL":
        return int(cn) >= high
    return int(cn) <= low


def open_text(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_vcf(path):
    header, records = [], []
    fh = open_text(path)
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            header.append(line)
        else:
            records.append(line.split("\t"))
    if fh is not sys.stdin:
        fh.close()
    return header, records


def info_dict(info):
    d = {}
    if info in ("", "."):
        return d
    for item in info.split(";"):
        k, _, v = item.partition("=")
        d[k] = v if _ else True
    return d


def sv_type(fields, info):
    t = str(info.get("SVTYPE", ""))
    if t.startswith("DUP") or fields[4].startswith("<DUP"):
        return "DUP"
    if t == "DEL" or fields[4] == "<DEL>":
        return "DEL"
    return t


def span(fields, info):
    try:
        return int(fields[1]), int(info["END"])
    except (KeyError, ValueError):
        return None


def first_gt(fields):
    if len(fields) < 10:
        return "./."
    keys = fields[8].split(":")
    if "GT" not in keys:
        return "./."
    values = fields[9].split(":")
    i = keys.index("GT")
    return values[i] if i < len(values) and values[i] else "./."


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


class Intervals:
    """Calls by (chromosome, type), for the reciprocal-overlap check."""

    def __init__(self):
        self.groups = {}

    def add(self, chrom, kind, start, end):
        g = self.groups.setdefault((chrom, kind), [[], [], 0])
        i = bisect.bisect(g[0], start)
        g[0].insert(i, start)
        g[1].insert(i, end)
        g[2] = max(g[2], end - start)

    def covers(self, chrom, kind, start, end, frac):
        g = self.groups.get((chrom, kind))
        if not g:
            return False
        length = max(end - start, 1)
        lo = bisect.bisect_left(g[0], start - g[2])
        hi = bisect.bisect_left(g[0], end)
        for s, e in zip(g[0][lo:hi], g[1][lo:hi]):
            overlap = min(e, end) - max(s, start)
            if overlap > 0 and overlap >= frac * max(e - s, length):
                return True
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--base", required=True, help="the VCF AnnotSV annotates (Manta's calls), or - for stdin")
    ap.add_argument("--base-name", default="manta", help="SVCALLER of the base's records (default manta)")
    ap.add_argument("--delly", help="Delly's VCF (step 19)")
    ap.add_argument("--sex", choices=("male", "female", ""), default="",
                    help="the sample's sex, for the normal copies on chrX and chrY (default: not given)")
    ap.add_argument("--par-bed", help="the pseudoautosomal regions (assets/par_grch38.bed)")
    ap.add_argument("--cnvpytor", help="CNVpytor's VCF (step 18)")
    ap.add_argument("--max-eval", type=float, default=0.01, help="CNVpytor e-val (pytorP1) below this (default 0.01)")
    ap.add_argument("--max-q0", type=float, default=0.5, help="CNVpytor q0 (pytorQ0) below this (default 0.5)")
    ap.add_argument("--min-overlap", type=float, default=0.5,
                    help="reciprocal overlap that makes two calls the same (default 0.5)")
    args = ap.parse_args(argv)

    header, base = read_vcf(args.base)
    n_samples = 0
    for h in header:
        if h.startswith("#CHROM"):
            n_samples = max(len(h.split("\t")) - 9, 0)

    known = Intervals()
    for f in base:
        info = info_dict(f[7])
        kind, sp = sv_type(f, info), span(f, info)
        if kind in ("DEL", "DUP") and sp:
            known.add(f[0], kind, sp[0], sp[1])
        f[7] = ("" if f[7] in ("", ".") else f[7] + ";") + "SVCALLER=" + args.base_name

    added = []
    summary = []

    def take(caller, records, keep):
        n_kept = n_added = 0
        for f in records:
            info = info_dict(f[7])
            kind, sp = sv_type(f, info), span(f, info)
            if kind not in ("DEL", "DUP") or not sp or not keep(f, info):
                continue
            n_kept += 1
            if known.covers(f[0], kind, sp[0], sp[1], args.min_overlap):
                continue
            known.add(f[0], kind, sp[0], sp[1])
            items = ["{0}={1}".format(k, info[k]) for k in KEPT if k in info and info[k] is not True]
            if "IMPRECISE" in info:
                items.append("IMPRECISE")
            items.append("SVCALLER=" + caller)
            row = [f[0], f[1], f[2], "N", "<{0}>".format(kind), f[5], "PASS", ";".join(items)]
            if n_samples:
                row += ["GT", first_gt(f)] + ["./."] * (n_samples - 1)
            added.append(row)
            n_added += 1
        return n_kept, n_added

    if args.delly:
        _, recs = read_vcf(args.delly)
        par = family_report.read_par(args.par_bed)
        no_depth = [0]

        def delly_keep(f, info):
            if f[6] != "PASS":
                return False
            sp = span(f, info)
            fmt = f[8].split(":") if len(f) > 9 else []
            vals = f[9].split(":") if len(f) > 9 else []
            if delly_depth_contradicts(sv_type(f, info), f[0], sp[0], sp[1], fmt, vals, args.sex,
                                       family_report.in_par(par, f[0], sp[0])):
                no_depth[0] += 1
                return False
            return True

        kept, n = take("delly", recs, delly_keep)
        summary.append("delly: {0} PASS DEL/DUP calls, {1} added, and {2} PASS calls longer than {3} bp set "
                       "aside because read depth (RDCN) does not show the change".format(kept, n, no_depth[0], DEPTH_MIN_SIZE))
    if args.cnvpytor:
        _, recs = read_vcf(args.cnvpytor)
        kept, n = take("cnvpytor", recs, lambda f, i: as_float(i.get("pytorP1")) < args.max_eval
                       and as_float(i.get("pytorQ0")) < args.max_q0)
        summary.append("cnvpytor: {0} DEL/DUP calls with e-val < {1:g} and q0 < {2:g}, {3} added".format(
            kept, args.max_eval, args.max_q0, n))

    contigs = []
    for h in header:
        if h.startswith("##contig=<ID="):
            contigs.append(h[len("##contig=<ID="):].split(",")[0].rstrip(">"))
    order = {c: i for i, c in enumerate(contigs)}
    records = base + added
    for f in records:
        if f[0] not in order:
            order[f[0]] = len(order)
    records.sort(key=lambda f: (order[f[0]], int(f[1])))

    present = {h.split(",")[0][len("##INFO=<ID="):] for h in header if h.startswith("##INFO=<ID=")}
    note = "Callers beside {0}: {1}".format(args.base_name, "; ".join(summary) if summary else "none given")
    out = sys.stdout
    for h in header:
        if h.startswith("##svcallers_added="):
            continue
        if h.startswith("#CHROM"):
            out.write('##svcallers_added="{0}"\n'.format(note))
            for key, line in STANDARD_INFO.items():
                if key not in present:
                    out.write(line + "\n")
            if added and n_samples and not any(x.startswith("##FORMAT=<ID=GT,") for x in header):
                out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(h + "\n")
    for f in records:
        out.write("\t".join(f) + "\n")

    sys.stderr.write(note + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
