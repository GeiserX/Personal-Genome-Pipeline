#!/usr/bin/env python3
"""The HLA section of the reports (step 08): bin/hla_alleles.tsv, what
bin/collect_summary.py derives from T1K's genotype TSV, and what
bin/render_report.py prints. Synthetic T1K rows only.

Checks:
  1. bin/hla_alleles.tsv: its columns, an allele group per row (B*27,
     B*57:01), kind drug or disease, PubMed ids as numbers, a citation, no
     duplicate allele, and every allele the report must cover (the drug
     alleles of docs/08-hla-typing.md, B*27, C*06:02, DQB1*06:02,
     DRB1*15:01 and the shared-epitope DRB1 alleles); each planted problem
     (a wrong kind, a PubMed id that is not a number, a duplicate, a missing
     allele, a three-field allele) is caught;
  2. the coeliac DQ line from DQA1 and DQB1: DQ2.5 likely cis with DRB1*03,
     DQ2.5 possible trans without it, DQ2.2, DQ8, none; the red control
     swaps DQA1*05:01 for DQA1*01:01 so DQ2.5 must not be derived; a
     low-confidence DQA1 gives "not determined";
  3. the DQB1*02 copies: 0, 1, 2 from two alleles, and 2 from T1K's one
     homozygous allele;
  4. the allele table: B*57:01 present and absent, B*27 matched by group,
     a candidate list of which only some are B*57:01 is "uncertain", a
     low-confidence locus is "not determined";
  5. a quality-0 allele with 12 candidates is stored with three and printed
     as "no confident call (12 candidates", in both reports;
  6. the HFE row carries its label in both reports;
  7. both reports print the derived lines and the not-a-diagnosis note, and
     every summary validates against tests/schema/summary.schema.json.

Run: python3 tests/test_hla_report.py
"""
import json
import os
import re
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.join(REPO, "tests", "schema"))
import collect_summary as cs  # noqa: E402
import render_report  # noqa: E402
import validate  # noqa: E402

TABLE = os.path.join(REPO, "bin", "hla_alleles.tsv")
COLUMNS = ["allele", "kind", "association", "pmids", "citation"]
REQUIRED = ["A*31:01", "B*15:02", "B*57:01", "B*58:01", "B*27", "C*06:02", "DQB1*06:02", "DRB1*15:01",
            "DRB1*01:01", "DRB1*04:01", "DRB1*04:04", "DRB1*04:05", "DRB1*04:08", "DRB1*14:02"]
ALLELE = re.compile(r"^[A-Z]+[0-9]?\*[0-9]{2}(:[0-9]{2,3})?$")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:800]}")
    if not ok:
        FAILS += 1


def table_errors(path):
    """Every schema problem of an allele table, as a list of strings."""
    errs, head, seen = [], None, set()
    with open(path) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if head is None:
                head = c
                if c != COLUMNS:
                    errs.append(f"header {c} is not {COLUMNS}")
                continue
            if len(c) != len(COLUMNS):
                errs.append(f"line {n}: {len(c)} columns")
                continue
            r = dict(zip(COLUMNS, c))
            if not ALLELE.match(r["allele"]):
                errs.append(f"line {n}: allele {r['allele']!r} is not a one- or two-field group")
            if r["kind"] not in ("drug", "disease"):
                errs.append(f"line {n}: kind {r['kind']!r}")
            if not r["association"].strip():
                errs.append(f"line {n}: no association")
            if not r["pmids"] or not all(p.isdigit() for p in r["pmids"].split(",")):
                errs.append(f"line {n}: pmids {r['pmids']!r}")
            if not r["citation"].strip():
                errs.append(f"line {n}: no citation")
            if r["allele"] in seen:
                errs.append(f"line {n}: duplicate {r['allele']}")
            seen.add(r["allele"])
    errs += [f"missing {a}" for a in REQUIRED if a not in seen]
    return errs


def t1k(rows):
    """T1K's genotype TSV from (gene, allele1, q1, allele2, q2); allele2 None
    writes T1K's one-allele row ('.', 0, -1)."""
    out = []
    for gene, a1, q1, a2, q2 in rows:
        if a2 is None:
            out.append(f"{gene}\t1\t{a1}\t30\t{q1}\t.\t0\t-1\n")
        else:
            out.append(f"{gene}\t2\t{a1}\t30\t{q1}\t{a2}\t20\t{q2}\n")
    return "".join(out)


def summary(work, name, rows):
    d = os.path.join(work, name, "S")
    os.makedirs(os.path.join(d, "hla_t1k"), exist_ok=True)
    with open(os.path.join(d, "hla_t1k", "S_hla_genotype.tsv"), "w") as f:
        f.write(t1k(rows))
    return cs.collect("S", d)


def main():
    schema = json.load(open(os.path.join(REPO, "tests", "schema", "summary.schema.json")))
    work = tempfile.mkdtemp()
    try:
        # 1. the table
        errs = table_errors(TABLE)
        check("bin/hla_alleles.tsv passes its schema and lists every required allele", not errs, errs)
        check("collect_summary reads the same rows", [r["allele"] for r in cs.hla_read_table()]
              == [l.split("\t")[0] for l in open(TABLE) if not l.startswith("#") and l.strip()][1:])
        good = open(TABLE).read()
        for desc, bad in (("a wrong kind", good.replace("\tdrug\t", "\tdrugs\t", 1)),
                          ("a PubMed id that is not a number", good.replace("\t21428769\t", "\tPMID21428769\t")),
                          ("a duplicate allele", good + "B*27\tdisease\tx\t1\ty\n"),
                          ("a missing required allele", good.replace("C*06:02\t", "C*07:02\t")),
                          ("a three-field allele", good.replace("B*58:01\t", "B*58:01:01\t"))):
            p = os.path.join(work, "bad.tsv")
            with open(p, "w") as f:
                f.write(bad)
            check(f"control: the table check catches {desc}", bool(table_errors(p)))

        # 2. coeliac DQ
        cis = summary(work, "cis", [("HLA-DQA1", "DQA1*05:01:01", 40, "DQA1*01:01:01", 40),
                                    ("HLA-DQB1", "DQB1*02:01:01", 40, "DQB1*05:01:01", 40),
                                    ("HLA-DRB1", "DRB1*03:01:01", 40, "DRB1*01:01:01", 40),
                                    ("HLA-B", "B*57:01:01", 40, "B*08:01:01", 40)])
        c = cis["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("DQA1*05:01 with DQB1*02:01 and DRB1*03 -> DQ2.5 present, likely cis; one DQB1*02",
              c["status"] == "derived" and c["dq2_5"] == "present" and c["dq2_5_phase"].startswith("likely cis")
              and c["dq2_2"] == "absent" and c["dq8"] == "absent" and c["dqb1_02_dose"] == 1, c)
        check("cis summary validates", not validate.validate(schema, json.loads(json.dumps(cis))))
        trans = summary(work, "trans", [("HLA-DQA1", "DQA1*05:05:01", 40, "DQA1*02:01:01", 40),
                                        ("HLA-DQB1", "DQB1*03:01:01", 40, "DQB1*02:02:01", 40),
                                        ("HLA-DRB1", "DRB1*11:01:01", 40, "DRB1*07:01:01", 40)])
        c = trans["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("DQA1*05:05 and DQB1*02:02 without DRB1*03 -> DQ2.5 possible trans (and DQ2.2 from the DR7 pair)",
              c["dq2_5"] == "present" and c["dq2_5_phase"].startswith("possible trans") and c["dq2_2"] == "present", c)
        dq22 = summary(work, "dq22", [("HLA-DQA1", "DQA1*02:01:01", 40, "DQA1*01:02:01", 40),
                                      ("HLA-DQB1", "DQB1*02:02:01", 40, "DQB1*06:02:01", 40)])
        c = dq22["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("DQA1*02:01 with DQB1*02:02 and no DQA1*05 -> DQ2.2 only",
              c["dq2_5"] == "absent" and c["dq2_2"] == "present" and c["dq8"] == "absent" and c["dq2_5_phase"] is None, c)
        dq8 = summary(work, "dq8", [("HLA-DQA1", "DQA1*03:01:01", 40, "DQA1*01:01:01", 40),
                                    ("HLA-DQB1", "DQB1*03:02:01", 40, "DQB1*05:01:01", 40)])
        c = dq8["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("DQA1*03 with DQB1*03:02 -> DQ8, no DQB1*02",
              c["dq8"] == "present" and c["dq2_5"] == "absent" and c["dqb1_02_dose"] == 0, c)
        red = summary(work, "red", [("HLA-DQA1", "DQA1*01:01:01", 40, "DQA1*01:02:01", 40),
                                    ("HLA-DQB1", "DQB1*02:01:01", 40, "DQB1*05:01:01", 40),
                                    ("HLA-DRB1", "DRB1*03:01:01", 40, "DRB1*01:01:01", 40)])
        c = red["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("red control: DQA1*05:01 swapped for DQA1*01:01 -> DQ2.5 must not be derived",
              c["dq2_5"] == "absent" and c["dq2_5_phase"] is None, c)
        low = summary(work, "low", [("HLA-DQA1", "DQA1*05:01:01", 0, "DQA1*01:01:01", 40),
                                    ("HLA-DQB1", "DQB1*02:01:01", 40, "DQB1*05:01:01", 40)])
        c = low["sections"]["hla"]["data"]["derived"]["coeliac"]
        check("a low-confidence DQA1 -> coeliac DQ not determined, the DQB1*02 count still given",
              c["status"] == "not_determined" and c["dq2_5"] is None and "HLA-DQA1 low confidence" in c["reason"]
              and c["dqb1_02_dose"] == 1, c)

        # 3. DQB1*02 copies
        two = summary(work, "two", [("HLA-DQB1", "DQB1*02:01:01", 40, "DQB1*02:02:01", 40)])
        hom = summary(work, "hom", [("HLA-DQB1", "DQB1*02:01:01", 40, None, None)])
        check("DQB1*02:01 / DQB1*02:02 -> two DQB1*02",
              two["sections"]["hla"]["data"]["derived"]["coeliac"]["dqb1_02_dose"] == 2)
        check("T1K's one homozygous DQB1*02:01 -> two DQB1*02",
              hom["sections"]["hla"]["data"]["derived"]["coeliac"]["dqb1_02_dose"] == 2)
        check("without DQA1 the coeliac line says DQA1 was not typed",
              hom["sections"]["hla"]["data"]["derived"]["coeliac"]["reason"] == "HLA-DQA1 not typed")

        # 4. the allele table
        def state(summ, allele):
            return {a["allele"]: a["state"] for a in summ["sections"]["hla"]["data"]["derived"]["alleles"]}[allele]
        check("B*57:01:01 -> B*57:01 present", state(cis, "B*57:01") == "present")
        b7 = summary(work, "b7", [("HLA-B", "B*07:02:01", 40, "B*08:01:01", 40)])
        check("B*57:01 absent with B*07:02 / B*08:01, and B*58:01 absent in the same sample",
              state(b7, "B*57:01") == "absent" and state(b7, "B*58:01") == "absent")
        grp = summary(work, "b27", [("HLA-B", "B*27:05:02", 40, "B*57:03:01,B*57:01:01", 40),
                                    ("HLA-C", "C*06:02:01", 0, "C*07:01:01", 40)])
        check("B*27:05:02 -> the B*27 group present", state(grp, "B*27") == "present")
        check("a candidate list of B*57:03 and B*57:01 -> B*57:01 uncertain", state(grp, "B*57:01") == "uncertain")
        check("a low-confidence HLA-C -> C*06:02 not determined", state(grp, "C*06:02") == "not_determined")
        check("an untyped DQB1 -> DQB1*06:02 not determined", state(grp, "DQB1*06:02") == "not_determined")

        # 5. quality-0 candidate list, 6. HFE label
        many = ",".join(f"DRB5*01:{i:02d}" for i in range(1, 13))
        cap = summary(work, "cap", [("HLA-DRB5", many, 0, None, None),
                                    ("HLA-HFE", "HFE*001:01", 40, "HFE*002:01", 40),
                                    ("HLA-DQA1", "DQA1*05:01:01", 40, "DQA1*01:01:01", 40),
                                    ("HLA-DQB1", "DQB1*02:01:01", 40, "DQB1*05:01:01", 40),
                                    ("HLA-B", "B*57:01:01", 40, "B*08:01:01", 40)])
        loci = {l["gene"]: l for l in cap["sections"]["hla"]["data"]["loci"]}
        check("a quality-0 allele with 12 candidates is stored with three and the count",
              loci["HLA-DRB5"]["alleles"] == ["DRB5*01:01,DRB5*01:02,DRB5*01:03"]
              and loci["HLA-DRB5"]["candidates"] == [12], loci["HLA-DRB5"])
        check("the HFE row has its label", loci["HLA-HFE"].get("label") == cs.HLA_HFE_LABEL, loci["HLA-HFE"])
        check("cap summary validates", not validate.validate(schema, json.loads(json.dumps(cap))),
              validate.validate(schema, json.loads(json.dumps(cap))))
        txt, html = render_report.text_report(cap), render_report.html_report(cap)
        check("text: the quality-0 locus is one short line",
              "HLA-DRB5 no confident call (12 candidates, quality 0); low confidence" in txt
              and "DRB5*01:04" not in txt, [l for l in txt.splitlines() if "DRB5" in l])
        check("HTML: the same", "no confident call (12 candidates, quality 0)" in html and "DRB5*01:04" not in html)
        check("both reports label the HFE row", txt.count("not the HFE C282Y/H63D genotype") == 1
              and html.count("not the HFE C282Y/H63D genotype") == 1)

        # 7. derived lines and the note
        want = ["Coeliac DQ: DQ2.5 present, phase unknown (HLA-DRB1 not typed); DQ2.2 absent; DQ8 absent",
                "DQB1*02 copies: 1",
                "B*57:01 present: abacavir hypersensitivity; flucloxacillin liver injury",
                "C*06:02 not determined (HLA-C not typed)"]
        check("text: the derived lines", all(w in txt for w in want), [l for l in txt.splitlines() if "DQ" in l or "*" in l])
        check("HTML: the derived lines", all(render_report.E(w) in html for w in want))
        check("both reports say it is information about alleles, not a diagnosis",
              txt.count("not a diagnosis") >= 1 and render_report.E(render_report.HLA_DERIVED_NOTE) in html)
        words = re.compile(r"\byou (should|must)\b|\bconsult\b|\btreat", re.I)
        check("no clinical advice in the derived lines or the note",
              not words.search(" ".join(render_report.hla_derived_lines(cap["sections"]["hla"]["data"]["derived"]))
                               + render_report.HLA_DERIVED_NOTE))
        check("control: the advice check finds advice", bool(words.search("you should consult")))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
