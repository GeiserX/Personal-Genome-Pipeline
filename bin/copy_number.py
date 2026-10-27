#!/usr/bin/env python3
"""copy_number.py: copy number at any locus of bin/cn_loci.tsv from short-read depth.

One method, one catalog row per locus: the depth over the locus's target
segments against the depth over single-copy control segments, as
`samtools bedcov` sums it. With T the target depth and C the control depth
(each combined as the row's target_depth and control_depth say),
    estimate = scale x T / C - offset
then guarded by a minimum control depth, flank agreement and the share of the
target depth from reads of mapping quality 0, and called as the row's `call`
column says. docs/cn-loci.md describes the columns and how to add a locus.

  copy_number.py bed [--catalog TSV] [--locus L ...]
      print the segments of the loci (all bedcov loci by default) as a BED
      (chrom, 0-based start, end, name), sorted by position within a locus
  copy_number.py call [--catalog TSV] --locus L --bedcov T [--bedcov-mapq1 T1]
                      [--json OUT.json] [--table OUT.tsv]
      read `samtools bedcov` over those segments (chrom, start, end, name,
      summed depth; the row's first bedcov pass) and, for a row with a MAPQ 0
      guard, the MAPQ >= 1 pass; write the locus's JSON block and, with
      --table, the step's metric<TAB>value table (bin/hba_copy_number.py and
      bin/lpa_kiv2.py lay those out). Exit 0 whatever the status; 2 on
      unreadable input; 3 (BROKEN_CONTROL) when a control segment has no
      read while the other controls have the depth the row needs: the row's
      coordinates are wrong, or the sample lacks that stretch, and no number
      is given.

Two empty tables mean the BAM has no reads on that contig at all (bedcov is
skipped there): status 'not assessed'. Standard library only: it runs in the
Python and Parascopy images and under tests/test_copy_number.py.
"""
import argparse
import csv
import importlib
import json
import math
import os
import sys

CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cn_loci.tsv")
COLUMNS = ["locus", "kind", "chrom", "targets", "controls", "details", "target_depth", "control_depth",
           "scale", "offset", "bedcov", "mapq0_min", "mapq0_max", "min_control_depth", "max_control_ratio",
           "call", "interval", "readings", "input", "owner", "cross_check", "pmids"]
KINDS = ("gene", "interval", "organelle")
INPUTS = ("bedcov", "cnvpytor-level")
Z95 = 1.96
# The step table of a locus, by the module that lays it out (same folder).
TABLES = {"HBA": "hba_copy_number", "LPA_KIV2": "lpa_kiv2"}


class CatalogError(ValueError):
    pass


class BrokenControl(ValueError):
    """A control segment with no read while the other controls are covered."""


def _num(row, col, default=None):
    v = row[col]
    if v == ".":
        return default
    try:
        return float(v)
    except ValueError:
        raise CatalogError(f"{row['locus']}: {col} is {v!r}, not a number or '.'")


def _segments(row, col):
    """[(name, start, end, label)] of a targets, controls or details cell."""
    if row[col] == ".":
        return []
    out = []
    for part in row[col].split(";"):
        bits = part.split(":", 2)
        try:
            name, (s, e) = bits[0], bits[1].split("-")
            start, end = int(s), int(e)
        except (IndexError, ValueError):
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} is not NAME:start-end[:label]")
        if not name or end <= start or start < 0:
            raise CatalogError(f"{row['locus']}: {col} segment {part!r} has no name or an empty span")
        out.append((name, start, end, bits[2] if len(bits) > 2 else ""))
    return out


def _readings(row):
    out = {}
    if row["readings"] == ".":
        return out
    for part in row["readings"].split("|"):
        k, sep, text = part.partition("=")
        if not sep or not k or not text:
            raise CatalogError(f"{row['locus']}: reading {part!r} is not KEY=text")
        out[k] = text
    return out


def _call(row):
    mode, sep, args = row["call"].partition(":")
    try:
        vals = [float(x) for x in args.split(",")]
    except ValueError:
        vals = []
    if not sep or mode not in ("round", "cut", "class") or not vals \
            or (mode == "round" and len(vals) != 1) or (mode == "class" and len(vals) != 2) \
            or vals != sorted(vals):
        raise CatalogError(f"{row['locus']}: call {row['call']!r} is not round:D, cut:A,B... or class:LOW,HIGH")
    return mode, vals


def parse_row(row):
    """The catalog row checked, with its parsed fields under '_' keys."""
    loc = row["locus"]
    if row["kind"] not in KINDS:
        raise CatalogError(f"{loc}: kind {row['kind']!r} (one of {', '.join(KINDS)})")
    if row["input"] not in INPUTS:
        raise CatalogError(f"{loc}: input {row['input']!r} (one of {', '.join(INPUTS)})")
    if not row["chrom"].startswith("chr"):
        raise CatalogError(f"{loc}: chrom {row['chrom']!r} (GRCh38 names, chr1 to chrM)")
    p = dict(row)
    p["_targets"], p["_controls"], p["_details"] = (_segments(row, c) for c in ("targets", "controls", "details"))
    p["_call"] = _call(row)
    p["_readings"] = _readings(row)
    if not p["_targets"]:
        raise CatalogError(f"{loc}: no target segment")
    if row["input"] == "bedcov":
        if not p["_controls"]:
            raise CatalogError(f"{loc}: a bedcov row needs control segments")
        if row["target_depth"] not in ("sum", "pooled") or row["control_depth"] not in ("flanks", "pooled"):
            raise CatalogError(f"{loc}: target_depth sum|pooled and control_depth flanks|pooled")
        p["_scale"], p["_offset"] = _num(row, "scale"), _num(row, "offset", 0.0)
        if not p["_scale"] or p["_scale"] <= 0:
            raise CatalogError(f"{loc}: scale must be a positive number")
        p["_min_control"] = _num(row, "min_control_depth", 0.0)
        p["_max_ratio"] = _num(row, "max_control_ratio")
        p["_mapq0_min"], p["_mapq0_max"] = _num(row, "mapq0_min"), _num(row, "mapq0_max")
        p["_passes"] = row["bedcov"].split(",")
        if (p["_mapq0_min"] is not None or p["_mapq0_max"] is not None) and len(p["_passes"]) != 2:
            raise CatalogError(f"{loc}: a MAPQ 0 guard needs two bedcov passes (all reads, MAPQ >= 1)")
        if p["_call"][0] == "class" and row["interval"] != "cv":
            raise CatalogError(f"{loc}: call class needs interval cv")
        if row["interval"] not in ("cv", "."):
            raise CatalogError(f"{loc}: interval {row['interval']!r} (cv or '.')")
    names = [n for n, *_ in p["_targets"] + p["_controls"] + p["_details"]]
    if len(names) != len(set(names)):
        raise CatalogError(f"{loc}: a segment name is used twice")
    return p


def load_catalog(path=None):
    """{locus: parsed row} of the catalog, in file order."""
    path = path or CATALOG
    with open(path) as f:
        lines = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or list(rows[0].keys()) != COLUMNS:
        raise CatalogError(f"{path}: the header must be {' '.join(COLUMNS)}")
    out, seen = {}, set()
    for r in rows:
        if None in r or None in r.values():
            raise CatalogError(f"{path}: row {r.get('locus')!r} has the wrong number of columns")
        if r["locus"] in out:
            raise CatalogError(f"{path}: locus {r['locus']} is listed twice")
        p = parse_row(r)
        for n, *_ in p["_targets"] + p["_controls"] + p["_details"]:
            if n in seen:
                raise CatalogError(f"{path}: segment name {n} is used by two loci")
            seen.add(n)
        out[r["locus"]] = p
    return out


def regions(row):
    """[(chrom, start, end, name)] of a row's segments, sorted by position."""
    segs = row["_controls"] + row["_targets"] + row["_details"]
    return sorted(((row["chrom"], s, e, n) for n, s, e, _ in segs), key=lambda r: (r[1], r[2]))


def read_bedcov(path, row):
    """{name: mean depth} from `samtools bedcov` output over the row's
    segments, which must all be there and nothing else; {} for an empty table."""
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
    if seen and seen != set(regions(row)):
        raise ValueError(f"{path}: its regions differ from the ones `copy_number.py bed --locus {row['locus']}` prints")
    return means


def _combine(mode, segs, m):
    if mode == "sum":
        return sum(m[n] for n, *_ in segs)
    if mode == "flanks":
        return sum(m[n] for n, *_ in segs) / len(segs)
    return sum(m[n] * (e - s) for n, s, e, _ in segs) / sum(e - s for _, s, e, _ in segs)


def _cv(values):
    """Coefficient of variation (sample standard deviation over the mean)."""
    mu = sum(values) / len(values)
    if mu <= 0 or len(values) < 2:
        return 0.0
    return math.sqrt(sum((v - mu) ** 2 for v in values) / (len(values) - 1)) / mu


def cut_copies(value, cuts):
    """The number of cut points a value reaches (cut:0.5,1.5 gives 0, 1 or 2)."""
    return sum(value >= c for c in cuts)


def reading(row, key):
    r = row["_readings"]
    if str(key) in r:
        return r[str(key)]
    if isinstance(key, int) and "more" in r and key > max(int(k) for k in r if k.isdigit()):
        return r["more"].format(n=key)
    return ""


def assess(row, m, q1=None):
    """The locus's block for mean depths m (every read) and q1 (reads of
    MAPQ >= 1, or None for a row without a MAPQ 0 guard).

    status is 'called' or 'not assessed'; reason names the guard that
    stopped it: no_data, low_control_depth, controls_disagree,
    low_mapq0_share, high_mapq0_share, below_floor, between_counts."""
    q1 = q1 if q1 is not None else {}
    tg, ct = row["_targets"], row["_controls"]
    b = {"locus": row["locus"], "status": "not assessed", "reason": None, "estimate": None, "total": None,
         "copies": None, "ci95": None, "class": None, "target_depth": None, "control_depth": None,
         "depths": dict(m), "mapq0_share": None, "details": {}, "reading": ""}
    if not m or any(n not in m for n, *_ in tg + ct):
        b["reason"] = "no_data"
        return b
    for n, *_ in ct:
        rest = [c for c in ct if c[0] != n]
        if m[n] == 0 and rest and _combine(row["control_depth"], rest, m) >= row["_min_control"] > 0:
            raise BrokenControl(f"BROKEN_CONTROL: control segment {n} of {row['locus']} has no read while the "
                                f"other controls have {_combine(row['control_depth'], rest, m):.1f}x; check its "
                                "coordinates in bin/cn_loci.tsv (or the sample lacks that stretch)")
    T, C = _combine(row["target_depth"], tg, m), _combine(row["control_depth"], ct, m)
    tq = _combine(row["target_depth"], tg, {n: q1.get(n, 0.0) for n, *_ in tg})
    share = 1 - tq / T if T > 0 else 0.0
    b.update(target_depth=T, control_depth=C, mapq0_share=share)
    low = min(m[n] for n, *_ in ct) if row["control_depth"] == "flanks" else C
    if low < row["_min_control"]:
        b["reason"] = "low_control_depth"
        return b
    if row["control_depth"] == "flanks" and row["_max_ratio"] is not None and len(ct) > 1:
        vals = [m[n] for n, *_ in ct]
        if max(vals) / min(vals) > row["_max_ratio"]:
            b["reason"] = "controls_disagree"
            return b
    total = row["_scale"] * T / C
    est = total - row["_offset"]
    b.update(total=total, estimate=est,
             details={n: row["_scale"] * m[n] / C for n, *_ in row["_details"] if n in m})
    # With no copy left there is no read to judge the aligner by (a row with
    # an offset has a floor of its own below, after this guard).
    no_copy = row["_offset"] == 0 and round(est) == 0
    if row["_mapq0_min"] is not None and not no_copy and share < row["_mapq0_min"]:
        b["reason"] = "low_mapq0_share"
        return b
    if row["_mapq0_max"] is not None and share > row["_mapq0_max"]:
        b["reason"] = "high_mapq0_share"
        return b
    if row["_offset"] > 0 and est <= 0:
        b["reason"] = "below_floor"
        return b
    mode, vals = row["_call"]
    if row["interval"] == "cv":
        rel = math.sqrt(_cv([m[n] for n, *_ in tg]) ** 2 / len(tg) + _cv([m[n] for n, *_ in ct]) ** 2 / len(ct))
        b["ci95"] = (max(0.0, total * math.exp(-Z95 * rel) - row["_offset"]),
                     total * math.exp(Z95 * rel) - row["_offset"])
    if mode == "round":
        n = round(est)
        if abs(est - n) > vals[0]:
            b["reason"] = "between_counts"
            return b
        b["copies"] = n
        b["reading"] = reading(row, n)
    elif mode == "cut":
        b["copies"] = cut_copies(est, vals)
        b["reading"] = reading(row, b["copies"])
    else:
        lo, hi = b["ci95"]
        b["class"] = "low" if hi < vals[0] else "high" if lo > vals[1] else "typical"
        b["reading"] = reading(row, b["class"])
    b["status"] = "called"
    return b


def _bed(args):
    cat = load_catalog(args.catalog)
    loci = args.locus or [k for k, r in cat.items() if r["input"] == "bedcov"]
    for loc in loci:
        if loc not in cat or cat[loc]["input"] != "bedcov":
            raise CatalogError(f"no bedcov locus {loc!r} in the catalog")
        for r in regions(cat[loc]):
            print("\t".join(map(str, r)))
    return 0


def _call_cmd(args):
    try:
        cat = load_catalog(args.catalog)
        if args.locus not in cat or cat[args.locus]["input"] != "bedcov":
            raise CatalogError(f"no bedcov locus {args.locus!r} in the catalog")
        row = cat[args.locus]
        m = read_bedcov(args.bedcov, row)
        guarded = row["_mapq0_min"] is not None or row["_mapq0_max"] is not None
        if guarded and not args.bedcov_mapq1:
            raise ValueError(f"{args.locus} has a MAPQ 0 guard: --bedcov-mapq1 is needed")
        q1 = read_bedcov(args.bedcov_mapq1, row) if args.bedcov_mapq1 else None
        if q1 is not None and bool(m) != bool(q1):
            raise ValueError(f"{args.bedcov_mapq1 if m else args.bedcov} is empty but the other table is not: "
                             "both must cover the regions, or both be empty (no reads on that contig)")
        block = assess(row, m, q1)
    except BrokenControl as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 3
    except (OSError, ValueError, KeyError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    rows, line = table(row, block)
    try:
        if args.json:
            with open(args.json, "w") as f:
                json.dump(block, f, indent=1, sort_keys=True)
                f.write("\n")
        if args.table:
            with open(args.table, "w") as f:
                f.write("metric\tvalue\n")
                for k, v in rows:
                    f.write(f"{k}\t{v}\n")
    except OSError as e:
        print(f"ERROR: cannot write: {e}", file=sys.stderr)
        return 2
    print(line)
    return 0


def table(row, block):
    """([(metric, value)], one summary line) of a locus's step table."""
    mod = TABLES.get(row["locus"])
    if mod:
        return importlib.import_module(mod).table(block)
    f = (lambda x: "." if x is None else f"{x:.2f}")
    rows = [("status", block["status"]), ("copies", "." if block["copies"] is None else block["copies"]),
            ("estimate", f(block["estimate"])), ("target_depth", f(block["target_depth"])),
            ("control_depth", f(block["control_depth"])), ("reason", block["reason"] or "."),
            ("reading", block["reading"] or ".")]
    return rows, f"{row['locus']}: {block['status']} ({block['reason'] or block['reading']})"


def main(argv=None):
    ap = argparse.ArgumentParser(prog="copy_number.py", description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bed", help="print the segments of the loci")
    b.add_argument("--catalog")
    b.add_argument("--locus", action="append")
    c = sub.add_parser("call", help="the copy number of one locus from samtools bedcov over its segments")
    c.add_argument("--catalog")
    c.add_argument("--locus", required=True)
    c.add_argument("--bedcov", required=True)
    c.add_argument("--bedcov-mapq1", help="samtools bedcov with the row's MAPQ >= 1 pass")
    c.add_argument("--json")
    c.add_argument("--table")
    args = ap.parse_args(argv)
    if args.cmd == "bed":
        try:
            return _bed(args)
        except (OSError, CatalogError) as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
    return _call_cmd(args)


if __name__ == "__main__":
    sys.exit(main())
