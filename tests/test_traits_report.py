#!/usr/bin/env python3
"""bin/traits_report.py (step 42) on synthetic bcftools output and on assets/traits.tsv.

Checks:
  1. a reference block expanded at a position (ALT '.', genotype G/G) reads as
     ref/ref with the table's hom_ref text; ref/alt and alt/alt read het and
     hom_alt, phased (C|T) as unphased, the reference allele printed first;
  2. a planted alternate allele changes the reading: the same position read
     as ref/ref and then as ref/alt gives two different readings, each the
     table's own;
  3. a position with no line is not_covered (never read as the reference); an
     allele the table does not describe and a deletion starting at the
     position are other_allele, not read; a half call that reaches the script
     (the bcftools genotype filter of step 42 drops it first, so in a run its
     position is not_covered) is other_allele too;
  4. a gVCF REF that is not the table's reference base stops the run (exit 1,
     the row named), so a gVCF on another build is not read;
  5. the APOE line: every pair of rs429358 and rs7412 genotypes gives its
     e-allele pair (e3/e3, e2/e3, e2/e2, e3/e4, e4/e4, e2/e4 with the e1/e3
     caveat, and the e1 pairs as not read); one APOE site not covered leaves
     the line not read; the line comes right after the second APOE row;
  6. the command line writes one row per table row plus the APOE line, the
     columns in order, and exits 1 on a table that fails the schema check.

Run: python3 tests/test_traits_report.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import traits_report as tr  # noqa: E402

TABLE = os.path.join(REPO, "assets", "traits.tsv")
SCRIPT = os.path.join(REPO, "bin", "traits_report.py")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def recs(*lines):
    """bcftools' CHROM POS REF ALT GENOTYPE lines -> read_genotypes' records."""
    d = tempfile.mkdtemp()
    try:
        p = os.path.join(d, "g.tsv")
        with open(p, "w") as f:
            f.write("".join("\t".join(map(str, l)) + "\n" for l in lines))
        return tr.read_genotypes(p)
    finally:
        shutil.rmtree(d)


def main():
    rows = tr.read_table(TABLE)
    by_id = {r["rsid"]: r for r in rows}
    lct = by_id["rs4988235"]   # chr2:135851076 G>A
    c, p, ref, alt = lct["chrom"], lct["pos"], lct["ref"], lct["alt"]
    at = lambda *ls: recs(*ls).get((c, p), [])  # noqa: E731

    # 1. genotypes and their readings
    g = tr.call(lct, at((c, p, ref, ".", f"{ref}/{ref}")))
    check("expanded reference block: ref/ref, hom_ref reading", g == (f"{ref}/{ref}", "called", lct["hom_ref"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{alt}/{ref}")))
    check("ref/alt: het reading, reference allele first", g == (f"{ref}/{alt}", "called", lct["het"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{ref}|{alt}")))
    check("phased ref|alt reads as ref/alt", g == (f"{ref}/{alt}", "called", lct["het"]), g)
    g = tr.call(lct, at((c, p, ref, alt, f"{alt}/{alt}")))
    check("alt/alt: hom_alt reading", g == (f"{alt}/{alt}", "called", lct["hom_alt"]), g)

    # 2. a planted alternate allele changes the reading
    before = tr.call(lct, at((c, p, ref, ".", f"{ref}/{ref}")))[2]
    after = tr.call(lct, at((c, p, ref, alt, f"{ref}/{alt}")))[2]
    check("planted alt: the reading changes from hom_ref to het", before == lct["hom_ref"] and after == lct["het"]
          and before != after, (before, after))

    # 3. not covered and other alleles
    check("no line: not_covered", tr.call(lct, [])[1] == "not_covered")
    other = next(b for b in "ACGT" if b not in (ref, alt))
    g = tr.call(lct, at((c, p, ref, other, f"{ref}/{other}")))
    check("an allele the table does not describe: other_allele, not read", g[1] == "other_allele"
          and g[2] == tr.OTHER_ALLELE, g)
    g = tr.call(lct, at((c, p, ref, alt, f"./{alt}")))
    check("a half call: other_allele", g[1] == "other_allele", g)
    g = tr.call(lct, at((c, p, ref + "T", ref, f"{ref}T/{ref}")))
    check("a deletion starting at the position: other_allele", g[1] == "other_allele", g)
    g = tr.call(lct, at((c, p, ref + "T", ref, f"{ref}T/{ref}"), (c, p, ref, alt, f"{ref}/{alt}")))
    check("a SNV record beside an indel at the position is the one read", g == (f"{ref}/{alt}", "called", lct["het"]), g)

    # 4. a REF that is not the GRCh38 base stops
    wrong = other
    try:
        tr.call(lct, at((c, p, wrong, alt, f"{wrong}/{alt}")))
        check("a gVCF REF that differs from the table stops the run", False, "no error")
    except tr.TableError as e:
        check("a gVCF REF that differs from the table stops the run, naming the row", "rs4988235" in str(e), e)

    # 5. APOE
    r4, r2 = by_id["rs429358"], by_id["rs7412"]
    want = {("T/T", "C/C"): "e3/e3", ("T/T", "C/T"): "e2/e3", ("T/T", "T/T"): "e2/e2",
            ("T/C", "C/C"): "e3/e4", ("C/C", "C/C"): "e4/e4", ("T/C", "C/T"): "e2/e4",
            ("C/C", "C/T"): "e1/e4", ("T/C", "T/T"): "e1/e2", ("C/C", "T/T"): "e1/e1"}
    for (g4, g2), pair in want.items():
        rr = recs((r4["chrom"], r4["pos"], "T", "C", g4), (r2["chrom"], r2["pos"], "C", "T", g2))
        out = tr.report(rows, rr)
        line = next(o for o in out if o["rsid"] == "rs429358+rs7412")
        check(f"APOE: rs429358 {g4} and rs7412 {g2} -> {pair}", line["genotype"] == pair and line["status"] == "called",
              line)
    rr = recs((r4["chrom"], r4["pos"], "T", "C", "T/C"), (r2["chrom"], r2["pos"], "C", "T", "C/T"))
    line = next(o for o in tr.report(rows, rr) if o["rsid"] == "rs429358+rs7412")
    check("APOE e2/e4: says it could also be the rare e1/e3", "e1/e3" in line["reading"], line)
    rr = recs((r4["chrom"], r4["pos"], "T", "C", "T/C"))
    out = tr.report(rows, rr)
    line = next(o for o in out if o["rsid"] == "rs429358+rs7412")
    check("APOE: one site not covered leaves the line not read", line["genotype"] == "."
          and line["status"] == "not_covered" and "Not read" in line["reading"], line)
    ids = [o["rsid"] for o in out]
    check("APOE: the line comes right after the second APOE row",
          ids.index("rs429358+rs7412") == max(ids.index("rs429358"), ids.index("rs7412")) + 1, ids)
    check("APOE: one line only", ids.count("rs429358+rs7412") == 1, ids)
    e4 = next(o for o in tr.report(rows, recs((r4["chrom"], r4["pos"], "T", "C", "C/C"),
                                                (r2["chrom"], r2["pos"], "C", "T", "C/C"))) if o["gene"] == "APOE"
              and "+" in o["rsid"])
    check("APOE e4/e4: a risk factor, not a diagnosis", "not a diagnosis" in e4["reading"], e4)

    # 6. the command line
    work = tempfile.mkdtemp()
    try:
        g = os.path.join(work, "g.tsv")
        with open(g, "w") as f:
            f.write(f"{c}\t{p}\t{ref}\t{alt}\t{ref}/{alt}\n")
        out = os.path.join(work, "S_traits.tsv")
        res = subprocess.run([sys.executable, SCRIPT, "--table", TABLE, "--genotypes", g, "--sample", "S", "--out", out],
                             capture_output=True, text=True)
        check("CLI: exits 0", res.returncode == 0, res.stderr)
        lines = open(out).read().splitlines()
        check("CLI: the header is the output columns", lines[0].split("\t") == tr.OUT_COLUMNS, lines[0])
        check("CLI: one row per table row plus the APOE line", len(lines) - 1 == len(rows) + 1, len(lines))
        row = next(l.split("\t") for l in lines[1:] if l.startswith("rs4988235\t"))
        check("CLI: the planted het reads het", row[5] == f"{ref}/{alt}" and row[6] == "called" and row[7] == lct["het"], row)
        check("CLI: a row with no line is not_covered",
              next(l for l in lines if l.startswith("rs671\t")).split("\t")[6] == "not_covered")
        bad = os.path.join(work, "bad.tsv")
        with open(TABLE) as f, open(bad, "w") as b:
            b.write(f.read().replace("\trs671\tALDH2\tG\tA\t", "\trs671\tALDH2\tG\tG\t"))
        res = subprocess.run([sys.executable, SCRIPT, "--table", bad, "--genotypes", g, "--sample", "S",
                              "--out", os.path.join(work, "x.tsv")], capture_output=True, text=True)
        check("CLI: a table that fails the schema check exits 1, naming the row",
              res.returncode == 1 and "rs671" in res.stderr and not os.path.exists(os.path.join(work, "x.tsv")),
              (res.returncode, res.stderr))
    finally:
        shutil.rmtree(work)

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
