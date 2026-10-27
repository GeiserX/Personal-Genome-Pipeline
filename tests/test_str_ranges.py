#!/usr/bin/env python3
"""The repeat ranges: the generated str_ranges.tsv, the report's cited subset and the docs.

1. str_ranges.tsv, which bin/strchive_convert.py writes beside the STRchive
   catalog (built here from the pinned files by tests/test_str_catalog.py):
   one row per ExpansionHunter locus of the catalog (82), each traced to its
   STRchive id; normal_max and pathogenic_min are the catalog's NormalMax and
   PathologicMin, empty where the catalog has no thresholds and the reason in
   no_thresholds; intermediate_min and intermediate_max are STRchive's own.
   Controls: a table with one row dropped, or one row twice, is caught.
2. bin/str_intermediate_ranges.tsv, the GeneReviews ranges the report names
   (HTT, FMR1, DMPK): the columns bin/collect_summary.py reads, one row per
   locus, whole numbers with min <= max ending right below the next class, a
   label and a GeneReviews citation; docs/interpreting-results.md prints the
   same range in the locus's row; every locus has a row in str_ranges.tsv,
   and the range lies above its normal_max and below its pathogenic_min, so a
   count the report calls intermediate is one Stranger calls pre_mutation
   (control: a planted HTT range of 20-35 is caught); collect_summary reads
   every row.
3. The "Key Disease Thresholds" table of docs/09-str-expansions.md: each
   row's normal maximum, pre-mutation band and full-mutation minimum are the
   catalog's (control: the table with HTT's minimum changed to 37 is caught).

Run: python3 tests/test_str_ranges.py
"""
import json
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests"))
import collect_summary  # noqa: E402
import test_str_catalog  # noqa: E402

TSV = os.path.join(REPO, "bin", "str_intermediate_ranges.tsv")
DOC = os.path.join(REPO, "docs", "interpreting-results.md")
STEP_DOC = os.path.join(REPO, "docs", "09-str-expansions.md")
COLUMNS = ["locus", "intermediate_min", "intermediate_max", "next_min", "next_class", "label", "citation"]
RANGE_COLUMNS = ["locus", "source_id", "normal_max", "intermediate_min", "intermediate_max", "pathogenic_min",
                 "no_thresholds"]
CITATION = re.compile(r"^GeneReviews: .+, https://www\.ncbi\.nlm\.nih\.gov/books/NBK\d+/$")

FAILS = []


def check(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


def table(path):
    head, out = None, []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        if head is None:
            head = c
            continue
        out.append(c)
    return head, out


def doc_rows():
    """Locus -> its row of the repeat table in docs/interpreting-results.md."""
    out = {}
    for line in open(DOC):
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(c) >= 6:
            out.setdefault(c[0].upper(), line)
    return out


def num(x):
    return int(x) if x != "" else None


def range_problems(rows, cat):
    """What is wrong with the generated ranges rows against the catalog (empty: nothing)."""
    bad = []
    by = {e["LocusId"]: e for e in cat}
    loci = [r[0] for r in rows]
    if sorted(loci) != sorted(by):
        bad.append(f"loci differ from the catalog's: missing {sorted(set(by) - set(loci))}, "
                   f"extra {sorted(set(loci) - set(by))}")
    if len(loci) != len(set(loci)):
        bad.append(f"twice: {sorted({x for x in loci if loci.count(x) > 1})}")
    for r in rows:
        if len(r) != len(RANGE_COLUMNS) or r[0] not in by:
            bad.append(f"bad row {r}")
            continue
        d, e = dict(zip(RANGE_COLUMNS, r)), by[r[0]]
        if d["source_id"] != e["SourceId"]:
            bad.append(f"{r[0]}: source_id {d['source_id']} is not {e['SourceId']}")
        if (num(d["normal_max"]), num(d["pathogenic_min"])) != (e.get("NormalMax"), e.get("PathologicMin")):
            bad.append(f"{r[0]}: {d['normal_max']}/{d['pathogenic_min']} are not the catalog's thresholds")
        if d["no_thresholds"] != e.get("NoThresholds", ""):
            bad.append(f"{r[0]}: no_thresholds differs from the catalog's NoThresholds")
    return bad


def part_generated(cat, cat_path):
    path = os.path.join(os.path.dirname(cat_path), "str_ranges.tsv")
    if not os.path.exists(path):
        check("the converter writes str_ranges.tsv beside the catalog", False, path)
        return {}
    head, rows = table(path)
    check("str_ranges.tsv: the header", head == RANGE_COLUMNS, head)
    check("str_ranges.tsv: 82 rows, one per ExpansionHunter locus of STRchive", len(rows) == 82, len(rows))
    bad = range_problems(rows, cat)
    check("str_ranges.tsv: one row per catalog locus, with the catalog's thresholds", not bad, bad)
    check("control: a table with one row dropped is caught", bool(range_problems(rows[1:], cat)))
    check("control: a table with one row twice is caught", bool(range_problems(rows + rows[:1], cat)))
    by = {r[0]: dict(zip(RANGE_COLUMNS, r)) for r in rows}
    loci = test_str_catalog.files()[2]
    rec = {x["id"]: x for x in json.load(open(loci))}
    same = all((num(d["intermediate_min"]), num(d["intermediate_max"]))
               == (rec[d["source_id"]].get("intermediate_min"), rec[d["source_id"]].get("intermediate_max"))
               for d in by.values())
    check("str_ranges.tsv: the intermediate ranges are STRchive's", same)
    check("ATXN2 (SCA2): normal to 28, intermediate 29-34, pathogenic from 35 (STRchive)",
          [by.get("ATXN2", {}).get(k) for k in RANGE_COLUMNS[2:6]] == ["28", "29", "34", "35"], by.get("ATXN2"))
    check("MUC1: no thresholds, and the reason", by.get("MUC1", {}).get("normal_max") == ""
          and "short reads cannot size it" in by.get("MUC1", {}).get("no_thresholds", ""), by.get("MUC1"))
    return by


def part_subset(gen):
    head, body = table(TSV)
    check("the header is the columns collect_summary.py reads", head == COLUMNS, head)
    check("at least HTT, FMR1 and DMPK are listed", {"HTT", "FMR1", "DMPK"} <= {r[0] for r in body},
          [r[0] for r in body])
    loci = [r[0] for r in body]
    check("one row per locus", len(loci) == len(set(loci)), loci)
    docs = doc_rows()

    def inside(locus, lo, hi):
        g = gen.get(locus, {})
        return bool(g.get("normal_max")) and int(g["normal_max"]) < lo and hi < int(g["pathogenic_min"])
    for r in body:
        if len(r) != len(COLUMNS):
            check(f"{r[0]}: {len(COLUMNS)} columns", False, r)
            continue
        d = dict(zip(COLUMNS, r))
        nums = [d["intermediate_min"], d["intermediate_max"], d["next_min"]]
        if not all(x.isdigit() for x in nums):
            check(f"{d['locus']}: whole-number range", False, nums)
            continue
        lo, hi, nxt = map(int, nums)
        check(f"{d['locus']}: min <= max ({lo}-{hi})", lo <= hi, (lo, hi))
        check(f"{d['locus']}: the range ends right below {d['next_class']} ({nxt})", hi + 1 == nxt, (hi, nxt))
        check(f"{d['locus']}: a label and a next class", bool(d["label"].strip() and d["next_class"].strip()), d)
        check(f"{d['locus']}: cited to a GeneReviews chapter", bool(CITATION.match(d["citation"])), d["citation"])
        row = docs.get(d["locus"].upper(), "")
        check(f"{d['locus']}: docs/interpreting-results.md prints {lo}-{hi} in its row",
              bool(re.search(rf"(?<![\d-]){lo}-{hi}(?![\d-])", row)), row.strip())
        g = gen.get(d["locus"], {})
        check(f"{d['locus']}: {lo}-{hi} has a row in str_ranges.tsv and lies above its normal maximum "
              f"({g.get('normal_max')}) and below its pathogenic minimum ({g.get('pathogenic_min')})",
              inside(d["locus"], lo, hi), g)
    check("control: a planted HTT range of 20-35 is outside STRchive's pre-mutation band", not inside("HTT", 20, 35))
    got = collect_summary.str_ranges(TSV)
    check("collect_summary.str_ranges reads every row", sorted(got) == sorted(loci), sorted(got))


def threshold_rows(text):
    """Gene -> (normal max, pre-mutation band, full-mutation minimum) of the Key Disease Thresholds table."""
    out, on = {}, False
    for line in text.splitlines():
        if line.startswith("## "):
            on = line.strip() == "## Key Disease Thresholds"
            continue
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        if on and line.startswith("|") and len(c) == 6 and c[1] not in ("Gene", "---"):
            out[c[1]] = (c[3], c[4], c[5])
    return out


def threshold_problems(rows, cat):
    by = {e["LocusId"]: e for e in cat}
    bad = [] if rows else ["no rows"]
    for gene, got in rows.items():
        e = by.get(gene, {})
        if "NormalMax" not in e:
            bad.append(f"{gene}: no thresholds in the catalog")
            continue
        n, p = e["NormalMax"], e["PathologicMin"]
        want = (f"≤{n}", f"{n + 1}-{p - 1}", f"≥{p}")
        if got != want:
            bad.append(f"{gene}: the docs say {got}, the catalog {want}")
    return bad


def part_docs(cat):
    text = open(STEP_DOC).read()
    rows = threshold_rows(text)
    check("docs/09-str-expansions.md: the Key Disease Thresholds table lists HTT, FMR1, FXN, C9ORF72, DMPK, ATXN1, "
          "ATXN2", {"HTT", "FMR1", "FXN", "C9ORF72", "DMPK", "ATXN1", "ATXN2"} <= set(rows), sorted(rows))
    bad = threshold_problems(rows, cat)
    check("docs/09-str-expansions.md: every row gives the catalog's thresholds", not bad, bad)
    planted = threshold_rows(text.replace("| HTT | CAG | ≤26 | 27-35 | ≥36 |", "| HTT | CAG | ≤26 | 27-35 | ≥37 |"))
    check("control: the table with HTT's pathogenic minimum changed to 37 is caught",
          planted != rows and bool(threshold_problems(planted, cat)), planted.get("HTT"))


def main():
    cat, cat_path = test_str_catalog.build_catalog()
    gen = part_generated(cat, cat_path)
    part_subset(gen)
    part_docs(cat)
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
