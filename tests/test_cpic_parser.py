#!/usr/bin/env python3
"""bin/pgx_parse.py, the one PharmCAT reader (step 27, the CPIC_LOOKUP module
and this test import or run the same file; there is no copy to keep in sync).

Checks:
  1. the real PharmCAT 3.2.0 report.json of the HG002 fixture
     (tests/fixtures/pharmcat/report-3.2.0.json) parses to genes; a gene PharmCAT
     lists several possible diplotypes for, with different phenotypes (CYP2C19,
     528 of them), is 'ambiguous' and never read as its first diplotype; a
     report that parses to zero genes makes `pgx_parse.py cpic-report` exit 1
     with a PARSING FAILED report, never an all-clear;
  2. on PharmCAT's own example report (pharmcat-docs-example.json, 10 genes
     with a non-normal phenotype) every such gene gets a drug list, with the
     recommendation for its own diplotype; synthetic genes check the fallbacks:
     the report's `drugs` section first, then `relatedDrugs`, then the static
     table, and otherwise a line saying it is not in the drug table;
  3. the flat (3.x), nested (2.x) and list layouts all parse;
  4. a gene PharmCAT could not call (or could not resolve) while pypgx did
     prints a warning naming every drug PharmCAT links to the gene, or a note
     when pypgx's call is normal, and the comparison marks it 'pypgx only'.

Run: python3 tests/test_cpic_parser.py
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import pgx_parse  # noqa: E402

FIXTURE = os.path.join(REPO, "tests", "fixtures", "pharmcat", "report-3.2.0.json")
EXAMPLE = os.path.join(REPO, "tests", "fixtures", "pharmcat", "pharmcat-docs-example.json")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:400]}")
    if not ok:
        FAILS += 1


def dip(a1, a2, phen):
    return {"sourceDiplotypes": [{"allele1": {"name": a1}, "allele2": {"name": a2},
                                  "label": f"{a1}/{a2}", "phenotypes": [phen]}]}


def run_report(work, name, data, pypgx=None):
    """Run the cpic-report command; returns (exit code, recommendations, phenotypes rows, comparison rows)."""
    d = os.path.join(work, name)
    os.makedirs(d)
    rep = os.path.join(d, "report.json")
    with open(rep, "w") as f:
        json.dump(data, f)
    args = ["cpic-report", "--sample", name, "--report", rep, "--outdir", d]
    comp = os.path.join(d, "comparison.tsv")
    if pypgx is not None:
        p = os.path.join(d, "pypgx.tsv")
        with open(p, "w") as f:
            f.write("Gene\tDiplotype\tPhenotype\tCNV_call\tSource\n")
            for g, (dl, ph) in pypgx.items():
                f.write(f"{g}\t{dl}\t{ph}\t.\tbam\n")
        args += ["--pypgx", p, "--comparison", comp]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        rc = pgx_parse.main(args)
    rec = open(os.path.join(d, f"{name}_cpic_recommendations.txt")).read()
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, f"{name}_phenotypes.tsv"))][1:]
    crow = [l.rstrip("\n").split("\t") for l in open(comp)][1:] if os.path.exists(comp) else []
    return rc, rec, rows, crow


def block(rec, gene):
    """The medications block of one gene in the recommendations text."""
    lines, on = [], False
    for line in rec.splitlines():
        if line.startswith(f"  {gene} -- "):
            on = True
        if on:
            lines.append(line)
            if "Action:" in line:
                break
    return "\n".join(lines)


def main():
    work = tempfile.mkdtemp()
    try:
        # 1. the real report
        with open(FIXTURE) as f:
            real = json.load(f)
        status, calls = pgx_parse.parse_genes(real)
        print(f"fixture: PharmCAT {real.get('pharmcatVersion')}, {len(calls)} genes: "
              + ", ".join(f"{c.gene} {c.diplotype} ({c.status})" for c in calls))
        check("the HG002 report parses (status OK)", status == "OK", status)
        check("the HG002 report yields genes", len(calls) >= 10, len(calls))
        check("the HG002 report calls at least one gene", any(c.called for c in calls))
        rc, rec, rows, _ = run_report(work, "HG002", real)
        check("cpic-report on the HG002 report exits 0", rc == 0, rc)
        check("one phenotypes row per gene, with a Status", len(rows) == len(calls) and all(len(r) == 4 for r in rows), rows[:3])
        by = {c.gene: c for c in calls}
        c19 = by.get("CYP2C19")
        check("HG002 CYP2C19: 528 possible diplotypes, status ambiguous",
              c19 is not None and len(c19.labels) == 528 and c19.status == "ambiguous", c19 and (len(c19.labels), c19.status))
        check("HG002 CYP2C19: not reported as its first diplotype's Normal Metabolizer",
              c19 is not None and c19.phenotype.startswith("ambiguous: ") and "(1 of 528 possible)" in c19.diplotype,
              c19 and (c19.diplotype, c19.phenotype))
        check("HG002 CYP2B6 is ambiguous too", by.get("CYP2B6") is not None and by["CYP2B6"].status == "ambiguous")
        check("the recommendations list the ambiguous genes in their own section",
              "Genes With More Than One Possible Result:" in rec and "  CYP2C19 -- 528 possible diplotypes" in rec, rec[:2000])
        check("the phenotypes table says ambiguous for CYP2C19",
              any(r[0] == "CYP2C19" and r[3] == "ambiguous" for r in rows), rows)
        guidance = pgx_parse.drug_guidance(real)
        check("the report's drugs section is read", sum(len(v) for v in guidance.values()) > 0)

        with open(EXAMPLE) as f:
            example = json.load(f)
        _, ecalls = pgx_parse.parse_genes(example)
        nonnormal = [c for c in ecalls if c.status == "non-normal"]
        check("example report: 10 genes with a non-normal phenotype", len(nonnormal) == 10, [c.gene for c in nonnormal])
        rc, rec, _, _ = run_report(work, "EXAMPLE", example)
        check("cpic-report on the example report exits 0", rc == 0, rc)
        for c in nonnormal:
            b = block(rec, c.gene)
            check(f"example {c.gene} ({c.phenotype}): a drug list", "    Drugs" in b and "not in the drug table" not in b, b)
        check("example CYP2D6 *1/*3 (Intermediate): codeine's recommendation for that diplotype",
              "- codeine [CPIC, Moderate]: Use codeine label recommended age- or weight-specific dosing. "
              "If no response and opioid use is warranted, consider a non-tramadol opioid." in block(rec, "CYP2D6"),
              block(rec, "CYP2D6"))

        stripped = dict(real, genes={g: {"geneSymbol": g} for g in real.get("genes", {})})
        rc, rec, rows, _ = run_report(work, "EMPTY", stripped)
        check("a report whose genes carry no diplotype exits 1", rc == 1, rc)
        check("and writes PARSING FAILED, never an all-clear",
              "PARSING FAILED" in rec and "all genes were successfully called" not in rec and rows == [], rec[:200])
        rc, rec, _, _ = run_report(work, "NOKEYS", {"metadata": {}})
        check("a report with no gene section at all exits 1", rc == 1 and "PARSING FAILED" in rec, rc)

        # 2. drug-list fallbacks
        data = {"genes": {
            "CYP2C19": dip("*2", "*2", "Poor Metabolizer"),
            "ABCG2": dict(dip("rs2231142 variant (T)", "rs2231142 variant (T)", "Poor Function"),
                          relatedDrugs=[{"name": "rosuvastatin"}]),
            "NEWGENE": dip("*1", "*7", "Decreased Function"),
            "CYP4F2": dip("*1", "*3", "Decreased Function"),
            "SLCO1B1": dip("*1", "*1", "Normal Function")},
            "drugs": {"CPIC Guideline Annotation": {"warfarin": {"name": "warfarin", "guidelines": [{"annotations": [
                {"phenotypes": {"CYP4F2": "Decreased Function"}, "classification": "Optional",
                 "drugRecommendation": "Increase the dose by 5-10%."}]}]}}}}
        rc, rec, rows, _ = run_report(work, "FALLBACK", data)
        check("fallbacks: exits 0", rc == 0, rc)
        check("CYP4F2: the report's own guidance (warfarin, CPIC Optional)",
              "- warfarin [CPIC, Optional]: Increase the dose by 5-10%." in block(rec, "CYP4F2"), block(rec, "CYP4F2"))
        check("ABCG2: relatedDrugs when the report has no guidance for it",
              "Drugs PharmCAT links to ABCG2" in block(rec, "ABCG2") and "rosuvastatin" in block(rec, "ABCG2"))
        check("CYP2C19: the static table when the report names no drug",
              "pipeline fallback table" in block(rec, "CYP2C19") and "clopidogrel" in block(rec, "CYP2C19"))
        check("NEWGENE: a line saying it is not in the drug table, not silence",
              "NEWGENE has a non-normal phenotype but is not in the drug table" in block(rec, "NEWGENE"), rec)
        check("SLCO1B1 normal: no medications block", block(rec, "SLCO1B1") == "")

        # PharmCAT lists an annotation for every diplotype the sample may have:
        # the recommendation shown must be the one for the called diplotype.
        ann = lambda label, phen, text: {"phenotypes": {"CYP2C19": phen}, "classification": "Strong",
                                         "drugRecommendation": text,
                                         "genotypes": [{"diplotypes": [{"gene": "CYP2C19", "label": label}]}]}
        data = {"genes": {"CYP2C19": dip("*1", "*2", "Intermediate Metabolizer")},
                "drugs": {"CPIC Guideline Annotation": {"clopidogrel": {"name": "clopidogrel", "guidelines": [{"annotations": [
                    ann("*1/*1", "Normal Metabolizer", "Standard dose."),
                    ann("*1/*2", "Intermediate Metabolizer", "Use an alternative antiplatelet."),
                    ann("*2/*2", "Poor Metabolizer", "Avoid clopidogrel.")]}]}}}}
        rc, rec, _, _ = run_report(work, "MATCH", data)
        b = block(rec, "CYP2C19")
        check("the recommendation is the called diplotype's (*1/*2), not another possible one's",
              "- clopidogrel [CPIC, Strong]: Use an alternative antiplatelet." in b
              and "Standard dose." not in b and "Avoid clopidogrel." not in b, b)

        # 3. layouts
        g2 = {"CYP2C19": dip("*1", "*2", "Intermediate Metabolizer"), "CYP2D6": dip("*1", "*1", "Normal Metabolizer")}
        for label, d, want in (
                ("flat 3.x", {"genes": g2}, 2),
                ("nested 2.x", {"genes": {"CPIC": g2, "DPWG": g2}}, 2),
                ("list", {"genes": [dict(g2["CYP2C19"], gene="CYP2C19")]}, 1),
                ("recommendationDiplotypes", {"genes": {"DPYD": {"recommendationDiplotypes": [
                    {"label": "c.1905+1G>A/Reference", "phenotypes": ["Intermediate Metabolizer"]}]}}}, 1)):
            st, cl = pgx_parse.parse_genes(d)
            check(f"layout {label}: {want} gene(s)", st == "OK" and len(cl) == want, (st, [c.gene for c in cl]))

        # 4. PharmCAT no result, pypgx call
        data = {"genes": {"CYP2D6": dict(dip("Unknown", "Unknown", "No Result"),
                                         relatedDrugs=[{"name": "codeine"}, {"name": "tramadol"}]),
                          "CYP2C9": dip("Unknown", "Unknown", "No Result"),
                          "CYP2C19": dip("*1", "*1", "Normal Metabolizer")},
                # a two-gene annotation that names CYP2D6 for one drug only
                "drugs": {"CPIC Guideline Annotation": {"amitriptyline": {"name": "amitriptyline", "guidelines": [
                    {"annotations": [{"phenotypes": {"CYP2C19": "Normal Metabolizer", "CYP2D6": "No Result"},
                                      "drugRecommendation": "x", "classification": "Optional"}]}]}}}}
        rc, rec, rows, comp = run_report(work, "PYPGX", data,
                                         pypgx={"CYP2D6": ("*1/*4", "Intermediate Metabolizer"),
                                                "CYP2C9": ("*1/*1", "Normal Metabolizer"),
                                                "CYP2C19": ("*1/*1", "Normal Metabolizer")})
        check("pypgx: exits 0", rc == 0, rc)
        check("pypgx: the warning names CYP2D6 and pypgx's call",
              "WARNING: PharmCAT has no result for CYP2D6, but pypgx (step 32) called *1/*4" in rec, rec)
        check("pypgx: the warning lists every drug PharmCAT links to CYP2D6, not only the one another gene matched",
              "Drugs affected by CYP2D6: codeine, tramadol" in rec, rec)
        check("pypgx: a normal pypgx call is a NOTE, not a WARNING",
              "NOTE: PharmCAT has no result for CYP2C9; pypgx (step 32) called *1/*1 (Normal Metabolizer)." in rec
              and "WARNING: PharmCAT has no result for CYP2C9" not in rec, rec)
        check("pypgx: CYP2D6 is in the uncallable list", "CYP2D6 -- No Result (not callable" in rec)
        check("pypgx: comparison marks CYP2D6 'pypgx only' and CYP2C19 concordant",
              ["CYP2D6", "Unknown/Unknown", "*1/*4", "pypgx only", "pypgx only"] in comp
              and ["CYP2C19", "*1/*1", "*1/*1", "Yes", "both"] in comp, comp)

        # an ambiguous PharmCAT gene is not compared as if it were its first diplotype
        data = {"genes": {"CYP2C19": {"sourceDiplotypes": [
            {"allele1": {"name": "*1"}, "allele2": {"name": "*1"}, "label": "*1/*1", "phenotypes": ["Normal Metabolizer"]},
            {"allele1": {"name": "*1"}, "allele2": {"name": "*2"}, "label": "*1/*2", "phenotypes": ["Intermediate Metabolizer"]}]}}}
        rc, rec, rows, comp = run_report(work, "AMBIG", data, pypgx={"CYP2C19": ("*1/*2", "Intermediate Metabolizer")})
        check("ambiguous: the phenotypes row says ambiguous", rows and rows[0][3] == "ambiguous", rows)
        check("ambiguous: compared as ambiguous, pypgx's call warned about",
              comp and comp[0][1] == "ambiguous (2 possible diplotypes)" and comp[0][3] == "pypgx only"
              and "WARNING: PharmCAT has no single result (2 possible diplotypes) for CYP2C19" in rec, (comp, rec))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
