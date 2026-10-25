#!/usr/bin/env python3
"""bin/lpa_kiv2.py: step 46's LPA KIV-2 repeat count from short-read depth,
and its Lp(a) section in both reports.

On `samtools bedcov` tables made up for each case (summed depth = mean depth
x exon length), 30x over the single-copy LPA exons, with 90% of the KIV-2
depth from reads of MAPQ 0 unless a case says otherwise:
  1. the arithmetic: the KIV-2 block holds N + 2 copies of the exon pair over
     seven reference pairs, so a block depth of 30 x (N + 2) / 14 gives N, for
     N = 12 (GRCh38's own six units on each chromosome), 20 and 40;
  2. the interval: equal depths on every exon give a zero-width interval;
     planted spread gives the interval sqrt(cv_kiv2^2 / 14 + cv_flank^2 / 15)
     on N + 2, computed here by hand, and the interval contains the estimate;
     twice the spread gives a wider interval;
  3. the class: a count whose whole interval is under 30 is low (29 with a
     narrow interval too), over 45 is high, one whose interval reaches into
     30 to 45 is typical (29 and 47 with a wide interval);
  4. under 10x over the single-copy exons is 'too little depth', with no
     count; no read of MAPQ 0 over KIV-2 (novoalign's default) is 'not
     assessed', naming the cause; two empty tables (a BAM with no chr6) is
     'not assessed';
  5. the CLI: `bed` prints the 29 regions, 14 KIV-2 exons of 2,394 bp and 15
     single-copy exons; `estimate` writes the metric table and exits 0; a
     table with other regions or a short row exits 2;
  6. assets/traits.tsv has the two LPA risk variants with the alleles
     bin/collect_summary.py counts, and the blood-test advice;
  7. both reports: a low count with one risk allele is flagged with both
     reasons and the advice; a typical count with no risk allele is not
     flagged; without step 42's table the variants are 'not read'; 'too
     little depth' is shown as such; no card and no text section without
     step 46's table (only the line in the list of steps that did not run).

Run: python3 tests/test_lpa_kiv2.py
"""
import math
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "lpa_kiv2.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import lpa_kiv2 as lpa  # noqa: E402

FAILS = 0
# The two LPA rows; named once here so no line of this file reads as an rsID with a genotype.
R1, R2 = "rs10455872", "rs3798220"
KIV = [n for *_, n in lpa.REGIONS if n.startswith("kiv2_")]
FLANK = [n for *_, n in lpa.REGIONS if n.startswith("flank_")]


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS += 1


def depths(n_copies, flank=30.0, kiv_spread=0.0, flank_spread=0.0):
    """Mean depth per region: the flanks at `flank`, the KIV-2 block at the
    depth n_copies gives; a spread alternates +x and -x around each mean."""
    k = flank * (n_copies + 2) / 14
    m = {}
    for i, n in enumerate(KIV):
        m[n] = k * (1 + (kiv_spread if i % 2 == 0 else -kiv_spread))
    for i, n in enumerate(FLANK):
        m[n] = flank * (1 + (flank_spread if i % 2 == 0 else -flank_spread))
    return m


def q1(m, share=0.9):
    return {n: v * (1 - share) if n.startswith("kiv2_") else v for n, v in m.items()}


def bedcov(path, m):
    with open(path, "w") as f:
        for c, s, e, n in lpa.REGIONS:
            f.write(f"{c}\t{s}\t{e}\t{n}\t{round(m[n] * (e - s))}\n")


def cv(vals):
    mu = sum(vals) / len(vals)
    return math.sqrt(sum((v - mu) ** 2 for v in vals) / (len(vals) - 1)) / mu


def main():
    # 1. arithmetic
    for n in (12, 20, 40):
        r = lpa.assess(depths(n), q1(depths(n)))
        check(f"arithmetic: a block at 30 x {n + 2} / 14 is {n} KIV-2 copies",
              r["status"] == "estimated" and r["kiv2_copies"] == f"{n:.1f}", r)
    # 2. interval
    r = lpa.assess(depths(20), q1(depths(20)))
    check("interval: equal depths give a zero-width interval", r["ci95_low"] == r["ci95_high"] == "20.0", r)
    m = depths(20, kiv_spread=0.10, flank_spread=0.15)
    r = lpa.assess(m, q1(m))
    rel = math.sqrt(cv([m[n] for n in KIV]) ** 2 / 14 + cv([m[n] for n in FLANK]) ** 2 / 15)
    klen = sum(e - s for _, s, e, n in lpa.REGIONS if n in KIV)
    flen = sum(e - s for _, s, e, n in lpa.REGIONS if n in FLANK)
    k = sum(m[n] * (e - s) for _, s, e, n in lpa.REGIONS if n in KIV) / klen
    f = sum(m[n] * (e - s) for _, s, e, n in lpa.REGIONS if n in FLANK) / flen
    total = 14 * k / f
    want = (f"{total - 2:.1f}", f"{total * math.exp(-1.96 * rel) - 2:.1f}", f"{total * math.exp(1.96 * rel) - 2:.1f}")
    got = (r["kiv2_copies"], r["ci95_low"], r["ci95_high"])
    check(f"interval: by hand {want}", got == want, got)
    check("interval: it contains the estimate",
          float(r["ci95_low"]) < float(r["kiv2_copies"]) < float(r["ci95_high"]), r)
    m2 = depths(20, kiv_spread=0.20, flank_spread=0.30)
    r2 = lpa.assess(m2, q1(m2))
    check("interval: twice the spread is wider",
          float(r2["ci95_high"]) - float(r2["ci95_low"]) > float(r["ci95_high"]) - float(r["ci95_low"]), (r, r2))
    # 3. class
    for n, sp, cls in ((20, 0.03, "low"), (29, 0.03, "low"), (29, 0.25, "typical"), (38, 0.03, "typical"),
                       (47, 0.25, "typical"), (60, 0.03, "high")):
        m = depths(n, kiv_spread=sp, flank_spread=sp)
        r = lpa.assess(m, q1(m))
        check(f"class: {n} copies ({r['ci95_low']} to {r['ci95_high']}) is {cls}", r["class"] == cls, r)
    r = lpa.assess(depths(20), q1(depths(20)))
    check("class: the low message says what a low count means",
          "a low count goes with small apo(a) isoforms and a higher Lp(a) level" in r["message"], r)
    # 4. not assessed
    r = lpa.assess(depths(20, flank=8.0), q1(depths(20, flank=8.0)))
    check("8x over the single-copy exons: too little depth, no count",
          r["status"] == "too little depth" and r["kiv2_copies"] == "." and "8.0x" in r["message"], r)
    r = lpa.assess(depths(20, flank=10.0), q1(depths(20, flank=10.0)))
    check("10x: estimated", r["status"] == "estimated", r)
    r = lpa.assess(depths(20), q1(depths(20), share=0.0))
    check("no read of MAPQ 0 over KIV-2: not assessed, naming the cause",
          r["status"] == "not assessed" and r["kiv2_copies"] == "." and "novoalign" in r["message"], r)
    r = lpa.assess(depths(20), q1(depths(20), share=0.25))
    check("a quarter of the depth from MAPQ 0: estimated", r["status"] == "estimated", r)
    r = lpa.assess({}, {})
    check("no chr6: not assessed", r["status"] == "not assessed" and "no chr6" in r["message"], r)

    # 5. CLI
    out = subprocess.run([sys.executable, SCRIPT, "bed"], capture_output=True, text=True)
    rows = [x.split("\t") for x in out.stdout.splitlines()]
    check("bed: 29 regions", out.returncode == 0 and len(rows) == 29, out)
    check("bed: 14 KIV-2 exons of 2,394 bp in total",
          sum(int(e) - int(s) for c, s, e, n in rows if n.startswith("kiv2_")) == 2394
          and sum(n.startswith("kiv2_") for *_, n in rows) == 14, rows)
    check("bed: 15 single-copy exons", sum(n.startswith("flank_") for *_, n in rows) == 15, rows)
    with tempfile.TemporaryDirectory() as t:
        a, b, o = (os.path.join(t, x) for x in ("a.tsv", "b.tsv", "o.tsv"))
        bedcov(a, depths(20))
        bedcov(b, q1(depths(20)))
        cli = [sys.executable, SCRIPT, "estimate", "--bedcov", a, "--bedcov-mapq1", b, "--out", o]
        out = subprocess.run(cli, capture_output=True, text=True)
        tab = dict(x.split("\t", 1) for x in open(o).read().splitlines()[1:]) if os.path.exists(o) else {}
        check("estimate: exit 0 and the metric table", out.returncode == 0 and tab.get("kiv2_copies") == "20.0"
              and tab.get("status") == "estimated" and "blood" in tab.get("advice", ""), (out, tab))
        open(a, "w").close()
        open(b, "w").close()
        out = subprocess.run(cli, capture_output=True, text=True)
        tab = dict(x.split("\t", 1) for x in open(o).read().splitlines()[1:])
        check("estimate: two empty tables (no chr6) exit 0, not assessed",
              out.returncode == 0 and tab.get("status") == "not assessed", (out, tab))
        bedcov(a, depths(20))
        with open(a, "a") as f:
            f.write("chr6\t1\t100\tother\t500\n")
        out = subprocess.run(cli, capture_output=True, text=True)
        check("estimate: other regions exit 2", out.returncode == 2 and "differ" in out.stderr, out)
        with open(a, "w") as f:
            f.write("chr6\t160646213\t160646395\n")
        out = subprocess.run(cli, capture_output=True, text=True)
        check("estimate: a short row exits 2", out.returncode == 2 and "expected chrom" in out.stderr, out)

    # 6. the traits table
    import collect_summary
    import render_report
    import traits_report
    table = {r["rsid"]: r for r in traits_report.read_table(os.path.join(REPO, "assets", "traits.tsv"))}
    for rs, alt in collect_summary.LPA_RISK_ALLELES.items():
        r = table.get(rs, {})
        check(f"traits.tsv: {rs} is an LPA row whose alt ({alt}) is the risk allele the report counts",
              r.get("gene") == "LPA" and r.get("alt") == alt, r)
        check(f"traits.tsv: {rs} says the actionable test is a blood Lp(a) measurement",
              "measure" in (r.get("note") or "").lower() and "blood" in (r.get("note") or ""), r)

    # 7. reports
    def sample(t, name, n_copies, flank=30.0, genotypes=None):
        d = os.path.join(t, name, "S")
        os.makedirs(os.path.join(d, "lpa"))
        a, b = os.path.join(t, name, "a.tsv"), os.path.join(t, name, "b.tsv")
        m = depths(n_copies, flank=flank, kiv_spread=0.03, flank_spread=0.03)
        bedcov(a, m)
        bedcov(b, q1(m))
        lpa.main(["estimate", "--bedcov", a, "--bedcov-mapq1", b, "--out", os.path.join(d, "lpa", "S_lpa_kiv2.tsv")])
        if genotypes:
            os.makedirs(os.path.join(d, "traits"))
            with open(os.path.join(d, "traits", "S_traits.tsv"), "w") as f:
                f.write("\t".join(traits_report.OUT_COLUMNS) + "\n")
                for rs, gt in genotypes.items():
                    row = {k: "." for k in traits_report.OUT_COLUMNS}
                    row.update(rsid=rs, gene="LPA", genotype=gt, status="called" if gt != "." else "not_covered")
                    f.write("\t".join(row[k] for k in traits_report.OUT_COLUMNS) + "\n")
        summ = collect_summary.collect("S", d)
        return summ, render_report.text_report(summ), render_report.html_report(summ)

    with tempfile.TemporaryDirectory() as t:
        summ, txt, html = sample(t, "low", 20, genotypes={R1: "A/G", R2: "T/T"})
        sec = summ["sections"].get("lpa", {})
        check("report: the section is read", sec.get("state") not in (None, "missing", "unreadable"), sec)
        check("report: a low count with one risk allele is flagged with both reasons",
              sec.get("data", {}).get("flag") == "high"
              and "1 risk allele at rs10455872 and rs3798220 and a low KIV-2 repeat count"
              in sec.get("data", {}).get("flag_text", ""), sec)
        check("report: the text section has the count, both variants and the advice",
              "## Lp(a): LPA KIV-2 Repeat Count and Risk Variants" in txt and "KIV-2 repeats (both chromosomes): 20.0" in txt
              and f"{R1}: A/G (1 risk allele)" in txt and f"{R2}: T/T (0 risk alleles)" in txt
              and lpa.ADVICE in txt, txt)
        check("report: the HTML card has the flag and the advice",
              "<h2>Lp(a): LPA KIV-2 Repeat Count</h2>" in html and "Genetic sign of high Lp(a)" in html
              and render_report.E(lpa.ADVICE) in html, html[-3000:])
        summ, txt, html = sample(t, "typ", 38, genotypes={R1: "A/A", R2: "T|T"})
        check("report: a typical count with no risk allele is not flagged",
              summ["sections"]["lpa"]["data"]["flag"] == "." and "No genetic sign of a high Lp(a)" in txt
              and "Genetic sign of high Lp(a)" not in html, txt)
        summ, txt, html = sample(t, "nov", 38)
        check("report: without step 42 the variants are not read",
              "rs10455872 and rs3798220: not read (step 42, the traits panel, did not run)" in txt, txt)
        summ, txt, html = sample(t, "thin", 38, flank=3.0, genotypes={R1: ".", R2: "C/C"})
        check("report: too little depth is said, and a not-covered variant is not read",
              "KIV-2 repeats: too little depth" in txt and f"{R1}: . (not read)" in txt
              and f"{R2}: C/C (2 risk alleles)" in txt, txt)
        os.makedirs(os.path.join(t, "none", "S"))
        nothing = collect_summary.collect("S", os.path.join(t, "none", "S"))
        check("report: no card and no text section without step 46",
              "<h2>Lp(a)" not in render_report.html_report(nothing) and "## Lp(a)" not in render_report.text_report(nothing))

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
