#!/usr/bin/env python3
"""bin/genotype_rules_indels.py (APOL1 G2 beside step 42's APOL1 rule), and its path
through bin/collect_summary.py and bin/render_report.py.

Checks:
  1. the G2 marker matches the pinned dbSNP record and the GRCh38 context; controls:
     a marker with ref and alt swapped, one at another position and one whose REF is
     not the GRCh38 bases each fail, naming the mistake;
  2. read_g2 on synthetic VCF records: the left-aligned dbSNP record and a
     right-shifted record of the same deletion each count (het 1, hom 2); a gVCF
     reference block over the span gives 0 copies, called, with its GQ; a
     variant-only VCF with no record gives no_call_record; a gVCF with no record
     gives not_covered; a no-call, a filtered G2 call, an SNV inside the span, a
     5-base deletion, a REF that is not GRCh38 and two het records of the one
     deletion give not_read; a record outside
     the span is skipped; control: a 6-base deletion just outside the span is not G2;
  3. combine: every G1 x G2 count gives its label and risk-copy reading; three risk
     copies are not read; a G1 rule not called is not read, unless G2 has two copies
     (G2/G2); G2 not read keeps the G1 result alone; a '.' result renders 'not read';
     no_call_record says the count rests on the absence of a call;
  4. collect_summary reads G2 from the sample's gVCF (bgzipped with a .tbi, and
     plain gzip without one) or its VCF, and puts it on the APOL1 rule line only;
     render_report shows "With G2" in the HTML and the text report; control: the
     same folder with a reference gVCF gives another result.

Run: python3 tests/test_genotype_rules_indels.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import genotype_rules_indels as m  # noqa: E402
import render_report as rr  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def ctx(a, b):
    """GRCh38 bases a..b (1-based, inclusive) from the pinned context."""
    return m.GRCh38_CONTEXT[a - m.CONTEXT_START:b - m.CONTEXT_START + 1]


def rec(pos, ref, alt, gt, filt="PASS", info=".", gq="30"):
    return ["chr22", str(pos), ".", ref, alt, "30", filt, info, "GT:GQ", f"{gt}:{gq}"]


LEFT = (36265995, "AATAATT", "A")                            # dbSNP's record
RIGHT = (36265999, ctx(36265999, 36266005), ctx(36265999, 36265999))   # the same deletion, shifted right
BLOCK = rec(36265900, "A", "<*>", "0/0", ".", "END=36266100", "50")


def marker_checks():
    check("1. the G2 marker matches dbSNP and GRCh38", m.marker_errors() == [], m.marker_errors())
    check("the G2 haplotype is GRCh38 with 6 bases deleted", len(m.GRCh38_CONTEXT) - len(m.G2_HAPLOTYPE) == 6)
    swapped = dict(m.G2, ref="A", alt="AATAATT")
    e = m.marker_errors(swapped)
    check("control: ref and alt swapped fails, named as swapped", any("ref and alt are swapped" in x for x in e), e)
    moved = dict(m.G2, pos=36265996)
    e = m.marker_errors(moved)
    check("control: another position fails against dbSNP and GRCh38",
          any("dbSNP chr22:36265995" in x for x in e) and any("not the GRCh38 bases" in x for x in e), e)
    e = m.marker_errors(dict(m.G2, ref="AATAATA"), dbsnp=("rs71785313", "chr22", 36265995, "AATAATA", "A"))
    check("control: a REF that is not the GRCh38 bases fails even when dbSNP agrees",
          any("not the GRCh38 bases" in x for x in e), e)


def read_checks():
    for name, (p, r, a) in (("left-aligned", LEFT), ("right-shifted", RIGHT)):
        for gt, n in (("0/1", 1), ("1|1", 2)):
            g = m.read_g2([rec(p, r, a, gt)], "vcf")
            check(f"2. {name} G2 {gt}: {n} copies", (g["status"], g["copies"]) == ("called", n), g)
    g = m.read_g2([rec(*LEFT, "0/1"), rec(*RIGHT, "0/1")], "vcf")
    check("two het records of the one deletion (left-aligned and shifted): not read, never two copies",
          (g["status"], g["copies"]) == ("not_read", None) and "more than one record" in g["text"], g)
    g = m.read_g2([BLOCK], "gvcf")
    check("a gVCF reference block over the span: 0 copies, called, its GQ",
          (g["status"], g["copies"], g["gq"]) == ("called", 0, "50"), g)
    g = m.read_g2([], "vcf")
    check("a VCF with no record: no_call_record, 0 copies", (g["status"], g["copies"]) == ("no_call_record", 0), g)
    g = m.read_g2([], "gvcf")
    check("a gVCF with no record: not_covered", (g["status"], g["copies"]) == ("not_covered", None), g)
    half = rec(36265900, "A", "<*>", "0/0", ".", "END=36266000", "50")
    g = m.read_g2([half], "gvcf")
    check("a gVCF block that stops inside the span: not_covered", g["status"] == "not_covered", g)
    snv = ctx(36266001, 36266001)
    other = "G" if snv != "G" else "C"
    bad = {"a no-call": rec(*LEFT, "./."), "a filtered G2 call": rec(*LEFT, "0/1", "RefCall"),
           "an SNV inside the span": rec(36266001, snv, other, "0/1"),
           "a 5-base deletion": rec(36265995, ctx(36265995, 36266000), "A", "0/1"),
           "a REF that is not GRCh38": rec(36265995, "CATAATT", "C", "0/1")}
    for desc, r in bad.items():
        g = m.read_g2([BLOCK, r], "gvcf")
        check(f"{desc}: not_read", (g["status"], g["copies"]) == ("not_read", None), g)
    g = m.read_g2([BLOCK, rec(36266100, "A", "G", "0/1"), ["chr21", "36265995", ".", "AATAATT", "A", ".",
                                                          "PASS", ".", "GT", "1/1"]], "gvcf")
    check("records outside the span or on another contig are skipped", (g["status"], g["copies"]) == ("called", 0), g)
    # control: a 6-base deletion right after the span (anchor 36266006) is not G2
    out = rec(36266006, ctx(36266006, 36266012), ctx(36266006, 36266006), "1/1")
    g = m.read_g2([BLOCK, out], "gvcf")
    check("control: a 6-base deletion outside the span does not count as G2", g["copies"] == 0, g)


def combine_checks():
    want = {("G0/G0", 0): "G0/G0", ("G0/G0", 1): "G0/G2", ("G0/G0", 2): "G2/G2", ("G0/G1", 0): "G0/G1",
            ("G0/G1", 1): "G1/G2", ("G1/G1", 0): "G1/G1"}
    for (g1, c), label in want.items():
        res = m.combine({"status": "called", "result": g1}, {"status": "called", "copies": c, "text": "t"})
        risk = label.count("G1") + label.count("G2")
        word = {0: "No G1 and no G2", 1: "One APOL1 risk copy", 2: "Two APOL1 risk copies"}[risk]
        check(f"3. {g1} + {c} G2: {label}", res["result"] == label and word in res["reading"], res)
    res = m.combine({"status": "called", "result": "G1/G1"}, {"status": "called", "copies": 1, "text": "t"})
    check("G1/G1 with a G2 copy (three risk copies): not read", res["result"] == "." and "never share" in res["reading"],
          res)
    res = m.combine({"status": "not_covered", "result": "."}, {"status": "called", "copies": 1, "text": "t"})
    check("the G1 rule not called: not read", res["result"] == ".", res)
    res = m.combine({"status": "not_covered", "result": "."}, {"status": "called", "copies": 2, "text": "t"})
    check("the G1 rule not called with two G2 copies: G2/G2, two risk copies",
          res["result"] == "G2/G2" and "Two APOL1 risk copies" in res["reading"], res)
    res = m.combine({"status": "called", "result": "G0/G1"}, {"status": "not_read", "copies": None, "text": "G2 x."})
    check("G2 not read: the G1 result stands alone, without repeating G2's text",
          res["result"] == "." and "stands on its own" in res["reading"] and "G2 x." not in res["reading"], res)
    for res in (m.combine({"status": "called", "result": "G1/G1"}, {"status": "called", "copies": 1, "text": "t"}),
                m.combine({"status": "not_covered", "result": "."}, {"status": "called", "copies": 1, "text": "t"})):
        shown = rr.g2_result({"g2": dict(res, status="called", text="t")})
        check(f"render_report: a '.' result shows 'not read' ({res['reading'][:30]})", shown == "not read", shown)
    res = m.combine({"status": "called", "result": "G0/G1"}, m.read_g2([], "vcf"))
    check("no_call_record: the label says no G2 call and the reading why",
          res["result"] == "G0/G1 (no G2 call)" and "no reference blocks" in res["reading"], res)


def bgzip_with_tbi(path, lines):
    """A bgzipped VCF and its .tbi from bgzip and tabix when they exist, else None
    (tests/e2e/genotype-rules-2-g2.sh covers the .tbi path with the pinned image)."""
    if not (shutil.which("bgzip") and shutil.which("tabix")):
        return None
    with open(path[:-3], "w") as f:
        f.writelines(lines)
    subprocess.run(["bgzip", "-f", path[:-3]], check=True)
    subprocess.run(["tabix", "-p", "vcf", path], check=True)
    return path


def report_checks():
    import collect_summary as cs
    import render_report as rr
    import traits_report as tr
    rules = tr.read_rules(os.path.join(REPO, "assets", "genotype_rules.tsv"))
    head = "##fileformat=VCFv4.2\n##contig=<ID=chr22>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    g1 = [("chr22", "36265860", "A", "G", "A/G", "40"), ("chr22", "36265988", "T", "G", "T/G", "40")]
    work = tempfile.mkdtemp()
    try:
        def folder(name, gvcf_lines, kind):
            d = os.path.join(work, name)
            os.makedirs(os.path.join(d, "traits"))
            os.makedirs(os.path.join(d, "vcf"))
            p = os.path.join(d, "traits", "S_genotype_rules.tsv")
            g = os.path.join(work, name + "_g.tsv")
            with open(g, "w") as f:
                f.writelines("\t".join(x) + "\n" for x in g1)
            tr.write_tsv(p, tr.rule_lines(rules, tr.read_genotypes(g)), tr.RULE_OUT_COLUMNS)
            v = os.path.join(d, "vcf", "S.g.vcf.gz" if kind == "gvcf" else "S.vcf.gz")
            with gzip.open(v, "wt") as f:
                f.write(head + "".join("\t".join(r) + "\n" for r in gvcf_lines))
            return d

        het = rec(*LEFT, "0/1")
        d = folder("het", [het], "gvcf")
        rows = cs.genotype_rules(d, "S")
        apol1 = next(r for r in rows if r["rule"] == "APOL1")
        others = [r["rule"] for r in rows if "g2" in r and r["rule"] != "APOL1"]
        check("4. collect_summary: G0/G1 plus a G2 copy in the gVCF reads G1/G2",
              apol1["g2"]["result"] == "G1/G2" and apol1["g2"]["source"] == "vcf/S.g.vcf.gz", apol1.get("g2"))
        check("only the APOL1 line carries g2", others == [], others)
        page = "\n".join(rr.rules_html({"rules": rows}, "S"))
        check("render_report HTML: With G2: G1/G2 and the high-risk reading",
              "With G2: G1/G2" in page and "Two APOL1 risk copies" in page, page[:400])
        text = "\n".join(rr.traits_text({"called": 0, "not_covered": 0, "other_allele": 0, "rows": [],
                                         "rules": rows}, "S"))
        check("render_report text: With G2 (rs71785313): G1/G2", "With G2 (rs71785313): G1/G2" in text, text)
        d = folder("ref", [BLOCK], "gvcf")
        apol1 = next(r for r in cs.genotype_rules(d, "S") if r["rule"] == "APOL1")
        check("control: a reference gVCF gives G0/G1", apol1["g2"]["result"] == "G0/G1", apol1.get("g2"))
        d = folder("vcf", [], "vcf")
        apol1 = next(r for r in cs.genotype_rules(d, "S") if r["rule"] == "APOL1")
        check("a VCF with no G2 record: G0/G1 (no G2 call), from vcf/S.vcf.gz",
              apol1["g2"]["result"] == "G0/G1 (no G2 call)" and apol1["g2"]["source"] == "vcf/S.vcf.gz",
              apol1.get("g2"))
        d = os.path.join(work, "none")
        os.makedirs(os.path.join(d, "traits"))
        shutil.copy(os.path.join(work, "het", "traits", "S_genotype_rules.tsv"), os.path.join(d, "traits"))
        apol1 = next(r for r in cs.genotype_rules(d, "S") if r["rule"] == "APOL1")
        check("no VCF at all: G2 not read", apol1["g2"]["status"] == "not_read" and apol1["g2"]["result"] == ".",
              apol1.get("g2"))
        tb = os.path.join(work, "tbi", "vcf", "S.g.vcf.gz")
        os.makedirs(os.path.dirname(tb))
        pad = [rec(36000000 + i * 1000, "A", "<*>", "0/0", ".", f"END={36000000 + i * 1000 + 999}") for i in range(300)]
        lines = [head] + ["\t".join(r) + "\n" for r in sorted(pad + [rec(*LEFT, "1/1")], key=lambda r: int(r[1]))]
        if bgzip_with_tbi(tb, lines):
            got = m.read_g2(cs.vcf_window(tb, "chr22", m.CONTEXT_START, m.SPAN[1]), "gvcf")
            check("vcf_window through the .tbi finds the G2 record (2 copies)", got["copies"] == 2, got)
        else:
            print("[SKIP] no bgzip and tabix here: the .tbi path runs where they are installed")
    finally:
        shutil.rmtree(work)


def main():
    marker_checks()
    read_checks()
    combine_checks()
    report_checks()
    print(f"\n{'FAILED' if FAILS else 'ALL PASSED'}: {FAILS} failure(s)")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
