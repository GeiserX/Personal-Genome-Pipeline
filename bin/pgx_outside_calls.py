#!/usr/bin/env python3
"""pgx_outside_calls.py: the PharmCAT outside-call file of one sample.

  pgx_outside_calls.py --calls OUT_CALLS.tsv --consensus OUT_CONSENSUS.tsv
                       [--hla T1K_genotype.tsv] [--pypgx PYPGX_summary.tsv]
                       [--cyrius CYRIUS.tsv] [--depth-check CYP2D6_depth_check.tsv]

writes the outside-call file PharmCAT reads with -po (gene<TAB>diplotype, one
line per gene, empty when nothing is passed on) and a consensus table (Gene,
Result, Outside_call, Reason, Evidence) that says, for HLA-A, HLA-B and
CYP2D6, what each caller said and what reached PharmCAT.

  HLA-A, HLA-B  T1K's two alleles (step 08), cut to two fields (*57:01), when
                both have a quality above 0; one allele reads as homozygous,
                as T1K reports a homozygous gene.
  CYP2D6        only when pypgx (step 32) and Cyrius (step 21) give the same
                diplotype and the depth check (bin/cyp2d6_depth_check.py)
                passed. One caller alone, a disagreement, a failed or missing
                depth check: 'indeterminate', and nothing reaches PharmCAT.
                PharmCAT itself calls no CYP2D6 from a VCF (3.4.0 reports it
                with callSource NONE), so it is not a second caller here.

Used by scripts/36-pgx-consensus.sh and the PGX_CONSENSUS module (bin/ is on
the task PATH); tests/test_pgx_outside_calls.py checks it. Standard library
only.
"""
import argparse
import csv
import re
import sys

CONSENSUS_HEADER = ["Gene", "Result", "Outside_call", "Reason", "Evidence"]
HLA_GENES = ("HLA-A", "HLA-B")
NO_CALL = {"", "none", "n/a", "na", "failed", "indeterminate", "not called", "unknown", "."}
ONE_CALLER = "one caller only: a second caller must agree (Cyrius, --tools cyrius)"


# --- HLA ---------------------------------------------------------------------

def two_fields(allele, gene):
    """'HLA-B*57:01:01:02' -> '*57:01'."""
    name = allele[len(gene):] if allele.startswith(gene) else allele
    if not name.startswith("*"):
        name = "*" + name
    return ":".join(name.split(":")[:2])


def read_t1k(path):
    """{gene: [(allele, quality)]} from T1K's genotype TSV: gene, number of
    distinct alleles, then allele, abundance, quality twice ('.', 0, -1 when
    there is no second allele)."""
    out = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0] not in HLA_GENES:
                continue
            alleles = []
            for i in (2, 5):
                if i + 2 < len(fields) and fields[i] not in ("", "."):
                    try:
                        q = float(fields[i + 2])
                    except ValueError:
                        q = -1.0
                    alleles.append((fields[i], q))
            out[fields[0]] = alleles
    return out


def hla_row(gene, t1k):
    """(consensus row, outside-call diplotype or None) for one HLA gene."""
    if t1k is None:
        return [gene, "not typed", "no", "HLA typing (step 08) did not run", "-"], None
    alleles = t1k.get(gene) or []
    if not alleles:
        return [gene, "not typed", "no", "T1K reported no allele", "-"], None
    evidence = "; ".join(f"{a} (quality {q:g})" for a, q in alleles)
    if any(q <= 0 for _, q in alleles):
        return [gene, "indeterminate", "no",
                "T1K gives an allele quality 0 or below (T1K says to ignore it)", evidence], None
    names = [two_fields(a, gene) for a, _ in alleles]
    reason = "T1K (step 08)"
    if len(names) == 1:
        names *= 2
        reason += ": one allele, read as homozygous"
    dip = "/".join(names)
    return [gene, dip, "yes", reason, evidence], dip


# --- CYP2D6 ------------------------------------------------------------------

def is_call(diplotype):
    d = (diplotype or "").strip()
    return d.lower() not in NO_CALL and "/" in d and ";" not in d


def same_diplotype(a, b):
    """Equal up to the order of the two alleles."""
    def key(d):
        return sorted(x.strip() for x in d.split("/"))
    return key(a) == key(b)


def read_pypgx(path):
    """pypgx's CYP2D6 diplotype ('' when its summary has no CYP2D6 row)."""
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if (row.get("Gene") or "").strip() == "CYP2D6":
                return (row.get("Diplotype") or "").strip()
    return ""


def read_cyrius(path):
    """(genotype, filter) of the first sample row of Cyrius's TSV."""
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            return (row.get("Genotype") or "").strip(), (row.get("Filter") or "").strip()
    return "", ""


def read_depth_check(path):
    """(status, message) from bin/cyp2d6_depth_check.py's TSV."""
    kv = {}
    with open(path) as f:
        for line in f:
            k, _, v = line.rstrip("\n").partition("\t")
            kv[k] = v
    return kv.get("status", ""), kv.get("message", "")


def cyp2d6_row(pypgx_path, cyrius_path, depth_path):
    """(consensus row, outside-call diplotype or None) for CYP2D6."""
    pg = read_pypgx(pypgx_path) if pypgx_path else None
    cy = read_cyrius(cyrius_path) if cyrius_path else None
    dc = read_depth_check(depth_path) if depth_path else None

    pg_ok = pg is not None and is_call(pg)
    cy_ok = cy is not None and is_call(cy[0]) and cy[1] == "PASS"
    if pg is None:
        pg_txt = "not run"
    elif pg_ok:
        pg_txt = pg
    else:
        pg_txt = f"no call ({pg or 'no CYP2D6 row'})"
    if cy is None:
        cy_txt = "not run"
    elif cy_ok:
        cy_txt = cy[0]
    else:
        cy_txt = f"no call ({', '.join(x for x in cy if x) or 'no row'})"
    dc_txt = "not run" if dc is None else (dc[0] or "unreadable")
    evidence = f"pypgx: {pg_txt}; Cyrius: {cy_txt}; depth check: {dc_txt}"

    def indeterminate(reason):
        return ["CYP2D6", "indeterminate", "no", reason, evidence], None

    if dc is not None and dc[0] != "ok":
        return indeterminate(dc[1] or "the CYP2D6 depth check failed")
    if not pg_ok and not cy_ok:
        return indeterminate("no caller made a call")
    if not (pg_ok and cy_ok):
        return indeterminate(ONE_CALLER)
    if not same_diplotype(pg, cy[0]):
        return indeterminate("pypgx and Cyrius disagree")
    if dc is None:
        return indeterminate("the CYP2D6 depth check did not run")
    return ["CYP2D6", pg, "yes", "pypgx and Cyrius agree; the depth check passed", evidence], pg


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--calls", required=True, help="outside-call file to write (PharmCAT -po)")
    p.add_argument("--consensus", required=True, help="consensus table to write")
    p.add_argument("--hla", help="T1K genotype TSV (step 08 / HLA_TYPING)")
    p.add_argument("--pypgx", help="pypgx summary TSV (step 32 / PYPGX)")
    p.add_argument("--cyrius", help="Cyrius TSV (step 21 / CYRIUS)")
    p.add_argument("--depth-check", help="CYP2D6 depth check TSV (bin/cyp2d6_depth_check.py)")
    args = p.parse_args(argv)

    try:
        t1k = read_t1k(args.hla) if args.hla else None
        rows, calls = [], []
        for gene in HLA_GENES:
            row, dip = hla_row(gene, t1k)
            rows.append(row)
            if dip:
                calls.append((gene, dip))
        row, dip = cyp2d6_row(args.pypgx, args.cyrius, args.depth_check)
        rows.append(row)
        if dip:
            calls.append(("CYP2D6", dip))
    except (OSError, UnicodeDecodeError, csv.Error) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    for c in calls:
        if not re.fullmatch(r"[^\t\n]+", c[1]):
            print(f"ERROR: refusing to write a malformed call {c!r}", file=sys.stderr)
            return 1
    with open(args.calls, "w") as f:
        f.writelines(f"{g}\t{d}\n" for g, d in calls)
    with open(args.consensus, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(CONSENSUS_HEADER)
        w.writerows(rows)
    for r in rows:
        print(f"{r[0]}: {r[1]} ({'passed to PharmCAT' if r[2] == 'yes' else 'not passed'}; {r[3]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
