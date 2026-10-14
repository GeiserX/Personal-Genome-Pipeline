#!/usr/bin/env python3
"""traits_report.py: the traits table of step 42, read from the gVCF genotypes.

  traits_report.py --table assets/traits.tsv --genotypes G.tsv --sample S --out S_traits.tsv
  traits_report.py --check-table assets/traits.tsv

G.tsv is what bcftools writes in step 42 (scripts/42-traits.sh and the
Nextflow TRAITS_SITES process run the same commands): one line per gVCF record
at a table position, after the reference blocks are expanded to one record per
base, as CHROM POS REF ALT GENOTYPE, the genotype as bases (C/T, or C|T when
phased). A position with no line had no reads: it is "not_covered", never read
as the reference.

Each table row becomes one output row with the genotype (the reference allele
first), a status (called, not_covered, other_allele) and the reading for that
genotype. After the two APOE rows comes one APOE line that combines them into
the e2, e3 and e4 alleles. A gVCF whose REF differs from the table's reference
base at a position is not GRCh38 (or the table is wrong): the run stops.

--check-table checks the schema of the table and exits 1 on the first list of
problems; tests/test_traits_table.py runs the same check. Standard library
only (it runs in the plain python image).
"""
import argparse
import os
import re
import sys

COLUMNS = ["chrom", "pos", "rsid", "gene", "ref", "alt", "trait", "hom_ref", "het", "hom_alt",
           "ancestry", "note", "pmids"]
OUT_COLUMNS = ["rsid", "gene", "chrom", "pos", "trait", "genotype", "status", "reading", "ancestry",
               "note", "pmids"]
CHROMS = {f"chr{n}" for n in range(1, 23)}
BASES = set("ACGT")

# APOE: rs429358 C (Arg112) and rs7412 C (Arg158) make e4; T and T make e2.
APOE_E4_SITE, APOE_E2_SITE = "rs429358", "rs7412"
APOE_ALLELES = {APOE_E4_SITE: ("T", "C"), APOE_E2_SITE: ("C", "T")}
# Odds ratios against e3/e3 for Alzheimer's disease, people of European
# descent from clinic and autopsy studies (Farrer et al., JAMA 1997, PMID 9343467).
APOE_READINGS = {
    (0, 0): ("e3/e3", "e3/e3, the most common pair: the reference group for Alzheimer's risk."),
    (0, 1): ("e2/e3", "e2/e3: odds of Alzheimer's disease about 0.6 times those of e3/e3."),
    (0, 2): ("e2/e2", "e2/e2: odds of Alzheimer's disease about 0.6 times those of e3/e3 (not significant in that study)."),
    (1, 0): ("e3/e4", "One e4: odds of late-onset Alzheimer's disease about 3.2 times those of e3/e3. "
                      "A risk factor, not a diagnosis."),
    (1, 1): ("e2/e4", "e2/e4: odds of Alzheimer's disease about 2.6 times those of e3/e3. Without phase this "
                      "pair could also be the rare e1/e3. A risk factor, not a diagnosis."),
    (2, 0): ("e4/e4", "Two e4: odds of late-onset Alzheimer's disease about 14.9 times those of e3/e3, the highest. "
                      "Still a risk factor, not a diagnosis."),
    (2, 1): ("e1/e4", "A pair with the rare e1 allele: not read."),
    (1, 2): ("e1/e2", "A pair with the rare e1 allele: not read."),
    (2, 2): ("e1/e1", "A pair with the rare e1 allele: not read."),
}
APOE_NOTE = ("Odds ratios from a meta-analysis of people of European descent (Farrer 1997); the e4 effect was "
             "weaker in African American and Hispanic people and stronger in Japanese people. Step 25's "
             "Alzheimer's score (PGS000334) includes the same two sites; docs/42-traits.md says why step 06 "
             "does not list APOE.")
NOT_COVERED = "No reads at this position in the gVCF: not read."
OTHER_ALLELE = "An allele this table does not describe at this position: not read."


class TableError(Exception):
    pass


def read_table(path):
    """The rows of the traits table as dicts. Lines starting with # are comments;
    the first other line is the header, which must be COLUMNS."""
    rows, header = [], None
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cells = line.split("\t")
            if header is None:
                if cells != COLUMNS:
                    raise TableError(f"{path}:{n}: header is {cells}, expected {COLUMNS}")
                header = cells
                continue
            if len(cells) != len(COLUMNS):
                raise TableError(f"{path}:{n}: {len(cells)} columns, expected {len(COLUMNS)}")
            row = dict(zip(COLUMNS, cells))
            row["line"] = n
            rows.append(row)
    if header is None:
        raise TableError(f"{path}: no header line")
    return rows


def check_rows(rows):
    """Every schema problem of the table rows, as messages (empty when it is valid)."""
    errs, ids, sites = [], set(), set()
    for r in rows:
        where = f"line {r['line']} ({r['rsid']})"
        if r["chrom"] not in CHROMS:
            errs.append(f"{where}: chrom '{r['chrom']}' is not chr1 to chr22")
        if not re.fullmatch(r"[1-9][0-9]*", r["pos"]):
            errs.append(f"{where}: pos '{r['pos']}' is not a positive integer")
        if not re.fullmatch(r"rs[1-9][0-9]*", r["rsid"]):
            errs.append(f"{where}: rsid is not rs<number>")
        if r["rsid"] in ids:
            errs.append(f"{where}: rsid appears twice")
        ids.add(r["rsid"])
        if (r["chrom"], r["pos"]) in sites:
            errs.append(f"{where}: position {r['chrom']}:{r['pos']} appears twice")
        sites.add((r["chrom"], r["pos"]))
        for k in ("ref", "alt"):
            if r[k] not in BASES:
                errs.append(f"{where}: {k} '{r[k]}' is not one of A, C, G, T")
        if r["ref"] == r["alt"]:
            errs.append(f"{where}: ref and alt are the same base")
        for k in ("gene", "trait", "hom_ref", "het", "hom_alt", "note"):
            if not r[k].strip():
                errs.append(f"{where}: {k} is empty")
        if r["ancestry"] not in ("yes", "no"):
            errs.append(f"{where}: ancestry '{r['ancestry']}' is not yes or no")
        if not re.fullmatch(r"[1-9][0-9]*(,[1-9][0-9]*)*", r["pmids"]):
            errs.append(f"{where}: pmids '{r['pmids']}' is not a comma-separated list of PubMed ids")
        text = " ".join(r[k] for k in COLUMNS)
        if "snpedia" in text.lower():
            errs.append(f"{where}: mentions SNPedia; the table takes no text from it")
        if '"' in text or not text.isascii():
            errs.append(f"{where}: a double quote or a non-ASCII character (the reports read the table as plain TSV)")
    by_id = {r["rsid"]: r for r in rows}
    for rs, (ref, alt) in APOE_ALLELES.items():
        r = by_id.get(rs)
        if r is None:
            errs.append(f"APOE row {rs} is missing")
        elif (r["ref"], r["alt"]) != (ref, alt) or r["chrom"] != "chr19":
            errs.append(f"APOE row {rs} must be chr19 {ref}>{alt}, is {r['chrom']} {r['ref']}>{r['alt']}")
    return errs


def read_genotypes(path):
    """{(chrom, pos): [(REF, [allele, allele]), ...]} from bcftools' CHROM POS REF ALT GT lines."""
    recs = {}
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            c = line.split("\t")
            if len(c) < 5:
                raise TableError(f"{path}:{n}: {len(c)} columns, expected CHROM POS REF ALT GENOTYPE")
            recs.setdefault((c[0], c[1]), []).append((c[2], re.split(r"[/|]", c[4])))
    return recs


def call(row, recs):
    """(genotype, status, reading) of one table row from its gVCF records."""
    if not recs:
        return ".", "not_covered", NOT_COVERED
    for ref, _ in recs:
        if ref[:1] != row["ref"]:
            raise TableError(f"{row['rsid']} at {row['chrom']}:{row['pos']}: the gVCF has REF {ref}, the table "
                             f"{row['ref']} (the GRCh38 base). Is the gVCF on GRCh38?")
    snv = [r for r in recs if all(len(a) == 1 for a in r[1])]
    _, alleles = (snv or recs)[0]
    order = {row["ref"]: 0, row["alt"]: 1}
    gt = "/".join(sorted(alleles, key=lambda a: (order.get(a, 2), a)))
    if len(alleles) != 2 or any(a not in order for a in alleles):
        return gt, "other_allele", OTHER_ALLELE
    n = alleles.count(row["alt"])
    return gt, "called", row[("hom_ref", "het", "hom_alt")[n]]


def apoe_line(rows, calls):
    """The APOE line: the e2/e3/e4 pair from rs429358 and rs7412."""
    by_id = {r["rsid"]: (r, c) for r, c in zip(rows, calls)}
    (r4, c4), (r2, c2) = by_id[APOE_E4_SITE], by_id[APOE_E2_SITE]
    line = {"rsid": f"{APOE_E4_SITE}+{APOE_E2_SITE}", "gene": "APOE", "chrom": "chr19",
            "pos": f"{r4['pos']},{r2['pos']}", "trait": "APOE alleles (e2, e3, e4)", "ancestry": "yes",
            "note": APOE_NOTE, "pmids": "8346443,9343467"}
    if c4[1] != "called" or c2[1] != "called":
        status = "not_covered" if "not_covered" in (c4[1], c2[1]) else "other_allele"
        return dict(line, genotype=".", status=status,
                    reading=f"Not read: needs a genotype at both {APOE_E4_SITE} and {APOE_E2_SITE}.")
    a = c4[0].split("/").count(r4["alt"])   # Arg112 copies (e4 change)
    b = c2[0].split("/").count(r2["alt"])   # Cys158 copies (e2 change)
    gt, reading = APOE_READINGS[(a, b)]
    return dict(line, genotype=gt, status="called", reading=reading)


def report(rows, recs):
    """The output rows: one per table row, then the APOE line after the second APOE row."""
    calls = [call(r, recs.get((r["chrom"], r["pos"]), [])) for r in rows]
    out = []
    apoe_done = False
    seen_apoe = set()
    for r, (gt, status, reading) in zip(rows, calls):
        out.append({"rsid": r["rsid"], "gene": r["gene"], "chrom": r["chrom"], "pos": r["pos"],
                    "trait": r["trait"], "genotype": gt, "status": status, "reading": reading,
                    "ancestry": r["ancestry"], "note": r["note"], "pmids": r["pmids"]})
        if r["rsid"] in APOE_ALLELES:
            seen_apoe.add(r["rsid"])
            if seen_apoe == set(APOE_ALLELES) and not apoe_done:
                out.append(apoe_line(rows, calls))
                apoe_done = True
    return out


def write_tsv(path, out):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write("\t".join(OUT_COLUMNS) + "\n")
        for o in out:
            f.write("\t".join(str(o[k]) for k in OUT_COLUMNS) + "\n")
    os.replace(tmp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table", help="assets/traits.tsv")
    ap.add_argument("--genotypes", help="bcftools' CHROM POS REF ALT GENOTYPE lines at the table positions")
    ap.add_argument("--sample")
    ap.add_argument("--out", help="the sample's traits table (TSV)")
    ap.add_argument("--check-table", metavar="TABLE", help="check the table's schema and exit")
    a = ap.parse_args(argv)
    try:
        if a.check_table:
            errs = check_rows(read_table(a.check_table))
            for e in errs:
                print(f"ERROR: {e}", file=sys.stderr)
            if not errs:
                print(f"{a.check_table}: schema OK")
            return 1 if errs else 0
        if not (a.table and a.genotypes and a.sample and a.out):
            ap.error("give --table, --genotypes, --sample and --out (or --check-table)")
        rows = read_table(a.table)
        errs = check_rows(rows)
        if errs:
            for e in errs:
                print(f"ERROR: {a.table}: {e}", file=sys.stderr)
            return 1
        out = report(rows, read_genotypes(a.genotypes))
    except (TableError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    write_tsv(a.out, out)
    counts = {}
    for o in out:
        counts[o["status"]] = counts.get(o["status"], 0) + 1
    print(f"{a.sample}: {len(rows)} table rows; " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
          + f" (APOE line included). Written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
