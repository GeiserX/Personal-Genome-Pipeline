#!/usr/bin/env python3
"""bin/str_intermediate_ranges.tsv: the intermediate repeat ranges the report names.

Checks every row:
  - the columns are the ones bin/collect_summary.py reads, one row per locus;
  - the range is whole numbers with min <= max, and it ends right below the
    next class (next_min = intermediate_max + 1), so it neither overlaps
    the disease range nor leaves a gap before it;
  - a label and a GeneReviews citation with its NCBI Bookshelf URL;
  - docs/interpreting-results.md prints the same range in the locus's row
    of its table, so the docs and the report cannot disagree;
and that bin/collect_summary.py reads the file (str_ranges) as written.

Run: python3 tests/test_str_ranges.py
"""
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402

TSV = os.path.join(REPO, "bin", "str_intermediate_ranges.tsv")
DOC = os.path.join(REPO, "docs", "interpreting-results.md")
COLUMNS = ["locus", "intermediate_min", "intermediate_max", "next_min", "next_class", "label", "citation"]
CITATION = re.compile(r"^GeneReviews: .+, https://www\.ncbi\.nlm\.nih\.gov/books/NBK\d+/$")

FAILS = []


def check(desc, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'} {desc}" + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(desc)


def rows():
    head, out = None, []
    for line in open(TSV):
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


def main():
    head, body = rows()
    check("the header is the columns collect_summary.py reads", head == COLUMNS, head)
    check("at least HTT, FMR1 and DMPK are listed", {"HTT", "FMR1", "DMPK"} <= {r[0] for r in body},
          [r[0] for r in body])
    loci = [r[0] for r in body]
    check("one row per locus", len(loci) == len(set(loci)), loci)
    docs = doc_rows()
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
    got = collect_summary.str_ranges(TSV)
    check("collect_summary.str_ranges reads every row", sorted(got) == sorted(loci), sorted(got))
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
