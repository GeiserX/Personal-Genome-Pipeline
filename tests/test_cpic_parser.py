#!/usr/bin/env python3
"""bin/pgx_parse.py, the one PharmCAT reader (step 27, the CPIC_LOOKUP module
and this test import or run the same file; there is no copy to keep in sync).

Checks:
  1. the real PharmCAT report.json of the HG002 fixture, from each PharmCAT
     release in REAL_REPORTS (tests/fixtures/pharmcat/report-<version>.json),
     parses to genes; a gene PharmCAT lists several possible diplotypes for,
     with different phenotypes (CYP2C19, 528 of them in 3.2.0), is 'ambiguous'
     and never read as its first diplotype; a report that parses to zero genes
     makes `pgx_parse.py cpic-report` exit 1 with a PARSING FAILED report,
     never an all-clear;
  2. on PharmCAT's own example report (pharmcat-docs-example.json) exactly 3
     genes count as non-normal (CYP2D6, HLA-B, NAT2); the 7 whose phenotype is
     n/a, 'no phenotype assigned', a genotype or 'Uncertain Susceptibility' are
     'unclassified', listed on their own and given no drug list; every
     non-normal gene gets a drug list, with the recommendation for its own
     diplotype; synthetic genes check the fallbacks:
     the report's `drugs` section first, then `relatedDrugs`, then the static
     table, and otherwise a line saying it is not in the drug table;
  3. the flat (3.x), nested (2.x) and list layouts all parse;
  4. a gene PharmCAT could not call (or could not resolve) while pypgx did
     prints a warning naming every drug PharmCAT links to the gene, or a note
     when pypgx's call is normal, and the comparison marks it 'pypgx only';
     an unreadable pypgx summary is reported, not fatal, and makes no comparison;
  5. a gene PharmCAT leaves 'Indeterminate' at gene level while its drug
     annotations give it one phenotype shows that phenotype, with a note saying
     where it came from, and takes its status; with no annotation, two that
     disagree, one that is Indeterminate itself, or possible diplotypes whose
     guidance differs, it stays Indeterminate and unclassified;
  6. a single-position gene PharmCAT labels 'rsX reference (C)/rsX variant
     (T)' matches pypgx's 'Reference/rsX', and a real difference still does
     not; the comparison carries both tools' phenotypes, so a gene pypgx has
     no phenotype rule for shows PharmCAT's; a CPIC recommendation over 200
     characters prints whole; CYP4F2 (and VKORC1) with no PharmCAT result
     takes the genotype PharmCAT read at the *3 (-1639) position as a
     labelled single-position fallback, and a position with no genotype is
     'not callable'.

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

FIXTURES = os.path.join(REPO, "tests", "fixtures", "pharmcat")
# One real HG002 report per PharmCAT release the parser is tested on: the
# version, and the genes it lists more than one possible diplotype for, with
# how many.
REAL_REPORTS = (
    ("3.2.0", {"CYP2C19": 528, "CYP2B6": 4}),
    ("3.4.0", {"CYP2C19": 528, "CYP2B6": 4}),
)
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
        # 1. the real reports, one per PharmCAT release
        for version, ambiguous in REAL_REPORTS:
            with open(os.path.join(FIXTURES, f"report-{version}.json")) as f:
                real = json.load(f)
            tag = f"{version}:"
            status, calls = pgx_parse.parse_genes(real)
            print(f"fixture: PharmCAT {real.get('pharmcatVersion')}, {len(calls)} genes: "
                  + ", ".join(f"{c.gene} {c.diplotype} ({c.status})" for c in calls))
            check(f"{tag} the report is PharmCAT {version}", real.get("pharmcatVersion") == version, real.get("pharmcatVersion"))
            check(f"{tag} the HG002 report parses (status OK)", status == "OK", status)
            check(f"{tag} the HG002 report yields genes", len(calls) >= 10, len(calls))
            check(f"{tag} the HG002 report calls at least one gene", any(c.called for c in calls))
            rc, rec, rows, _ = run_report(work, f"HG002_{version}", real)
            check(f"{tag} cpic-report on the HG002 report exits 0", rc == 0, rc)
            check(f"{tag} one phenotypes row per gene, with a Status",
                  len(rows) == len(calls) and all(len(r) == 4 for r in rows), rows[:3])
            by = {c.gene: c for c in calls}
            check(f"{tag} the genes with more than one possible diplotype are {sorted(ambiguous)}",
                  sorted(c.gene for c in calls if len(c.labels) > 1) == sorted(ambiguous),
                  [(c.gene, len(c.labels)) for c in calls if len(c.labels) > 1])
            for gene, n in sorted(ambiguous.items()):
                g = by.get(gene)
                check(f"{tag} HG002 {gene}: {n} possible diplotypes, status ambiguous",
                      g is not None and len(g.labels) == n and g.status == "ambiguous", g and (len(g.labels), g.status))
                check(f"{tag} HG002 {gene}: not reported as its first diplotype's phenotype",
                      g is not None and g.phenotype.startswith("ambiguous: ") and f"(1 of {n} possible)" in g.diplotype,
                      g and (g.diplotype, g.phenotype))
                check(f"{tag} the recommendations list {gene} with the ambiguous genes",
                      "Genes With More Than One Possible Result:" in rec and f"  {gene} -- {n} possible diplotypes" in rec,
                      rec[:2000])
                check(f"{tag} the phenotypes table says ambiguous for {gene}",
                      any(r[0] == gene and r[3] == "ambiguous" for r in rows), rows)
            guidance = pgx_parse.drug_guidance(real)
            check(f"{tag} the report's drugs section is read", sum(len(v) for v in guidance.values()) > 0)

        with open(EXAMPLE) as f:
            example = json.load(f)
        _, ecalls = pgx_parse.parse_genes(example)
        nonnormal = [c for c in ecalls if c.status == "non-normal"]
        check("example report: 3 genes with a non-normal phenotype (CYP2D6, HLA-B, NAT2)",
              sorted(c.gene for c in nonnormal) == ["CYP2D6", "HLA-B", "NAT2"], [c.gene for c in nonnormal])
        unclassified = ["CACNA1S", "CFTR", "CYP4F2", "IFNL3", "MT-RNR1", "RYR1", "VKORC1"]
        check("example report: n/a, 'no phenotype assigned', a genotype and 'Uncertain Susceptibility' are unclassified",
              sorted(c.gene for c in ecalls if c.status == "unclassified") == unclassified,
              [(c.gene, c.phenotype, c.status) for c in ecalls])
        rc, rec, erows, _ = run_report(work, "EXAMPLE", example)
        check("cpic-report on the example report exits 0", rc == 0, rc)
        meds = rec.split("Affected Medications:")[1].split("\n\n", 1)[1].split("Called Genes Without")[0]
        check("example: the medications section lists only the 3 non-normal genes",
              sorted(l.split(" -- ")[0].strip() for l in meds.splitlines() if l.startswith("  ") and " -- " in l
                     and not l.startswith("   ")) == ["CYP2D6", "HLA-B", "NAT2"], meds[:1500])
        own = rec.split("Called Genes Without a Function Phenotype:")[-1].split("Uncallable Genes:")[0]
        check("example: the unclassified genes are listed in their own section, VKORC1 with its genotype",
              all(f"  {g} -- " in own for g in unclassified) and "  VKORC1 -- -1639 GG (" in own, own)
        check("example: the phenotypes table says unclassified for them",
              sorted(r[0] for r in erows if r[3] == "unclassified") == unclassified, erows)
        for phen, want in (("Likely Poor Metabolizer", True), ("Possible Intermediate Metabolizer", True),
                           ("Ultrarapid Metabolizer", True), ("Possible Decreased Function", True),
                           ("Increased Function", True), ("Deficient", True), ("Deficient with CNSHA", True),
                           ("Malignant Hyperthermia Susceptibility", True), ("*57:01 positive", True),
                           ("increased risk of aminoglycoside-induced hearing loss", True),
                           ("ivacaftor responsive in CF patients", True),
                           ("ivacaftor non-responsive in CF patients", False), ("Indeterminate", False),
                           ("n/a", False), ("no phenotype assigned", False), ("-1639 GA", True),
                           ("-1639 AA", True), ("-1639 GG", False),
                           ("Uncertain Susceptibility", False), ("*58:01 negative", False)):
            check(f"is_non_normal({phen!r}) is {want}", pgx_parse.is_non_normal(phen) is want)
        ga = json.loads(json.dumps(example))
        ga["genes"]["VKORC1"]["sourceDiplotypes"][0].update(
            label="rs9923231 reference (C)/rs9923231 variant (T)", phenotypes=["-1639 GA"])
        _, gcalls = pgx_parse.parse_genes(ga)
        check("example with VKORC1 -1639 GA: VKORC1 counts as non-normal",
              [c.status for c in gcalls if c.gene == "VKORC1"] == ["non-normal"],
              [(c.gene, c.phenotype, c.status) for c in gcalls if c.gene == "VKORC1"])
        rc, grec, _, _ = run_report(work, "EXAMPLE_GA", ga)
        check("example with VKORC1 -1639 GA: its drug list names warfarin",
              rc == 0 and "    Drugs" in block(grec, "VKORC1") and "warfarin" in block(grec, "VKORC1"),
              block(grec, "VKORC1"))
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

        # The phenotype fallback (annotations naming no diplotype) matches the
        # called phenotype exactly: "Intermediate Metabolizer" guidance is not
        # the guidance for "Likely Intermediate Metabolizer".
        pann = lambda phen, text: {"phenotypes": {"CYP2C19": phen}, "classification": "Strong",
                                   "drugRecommendation": text}
        data = {"genes": {"CYP2C19": dip("*1", "*9", "Likely Intermediate Metabolizer")},
                "drugs": {"CPIC Guideline Annotation": {"clopidogrel": {"name": "clopidogrel", "guidelines": [{"annotations": [
                    pann("Intermediate Metabolizer", "Guidance for IM."),
                    pann("Likely Intermediate Metabolizer", "Guidance for likely IM.")]}]}}}}
        rc, rec, _, _ = run_report(work, "PHENMATCH", data)
        b = block(rec, "CYP2C19")
        check("phenotype fallback: the exact phenotype's guidance only, not a substring match",
              "Guidance for likely IM." in b and "Guidance for IM." not in b, b)

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
              ["CYP2D6", "Unknown/Unknown", "*1/*4", "pypgx only", "pypgx only",
               "No Result", "Intermediate Metabolizer"] in comp
              and ["CYP2C19", "*1/*1", "*1/*1", "Yes", "both", "Normal Metabolizer", "Normal Metabolizer"] in comp,
              comp)

        # an ambiguous PharmCAT gene is not compared as if it were its first diplotype
        data = {"genes": {"CYP2C19": {"sourceDiplotypes": [
            {"allele1": {"name": "*1"}, "allele2": {"name": "*1"}, "label": "*1/*1", "phenotypes": ["Normal Metabolizer"]},
            {"allele1": {"name": "*1"}, "allele2": {"name": "*2"}, "label": "*1/*2", "phenotypes": ["Intermediate Metabolizer"]}]}}}
        rc, rec, rows, comp = run_report(work, "AMBIG", data, pypgx={"CYP2C19": ("*1/*2", "Intermediate Metabolizer")})
        check("ambiguous: the phenotypes row says ambiguous", rows and rows[0][3] == "ambiguous", rows)
        check("ambiguous: compared as ambiguous, pypgx's call warned about",
              comp and comp[0][1] == "ambiguous (2 possible diplotypes)" and comp[0][3] == "pypgx only"
              and "WARNING: PharmCAT has no single result (2 possible diplotypes) for CYP2C19" in rec, (comp, rec))

        # an unreadable pypgx summary: the CPIC report is still written, it says
        # the comparison could not be made, and no comparison file is written
        d = os.path.join(work, "BADPYPGX")
        os.makedirs(d)
        rep, bad, comp = (os.path.join(d, n) for n in ("report.json", "pypgx.tsv", "comparison.tsv"))
        with open(rep, "w") as f:
            json.dump({"genes": {"CYP2C19": dip("*2", "*2", "Poor Metabolizer")}}, f)
        with open(bad, "wb") as f:
            f.write(b"Gene\tDiplotype\n\xff\xfe\x00broken\n")
        with open(comp, "w") as f:  # a comparison left by an earlier run
            f.write("Gene\tPharmCAT_diplotype\tpypgx_diplotype\tMatch\tCalled_by\nCYP2C19\t*1/*1\t*1/*1\tYes\tboth\n")
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                rc = pgx_parse.main(["cpic-report", "--sample", "BADPYPGX", "--report", rep, "--outdir", d,
                                     "--pypgx", bad, "--comparison", comp])
            except Exception as e:  # the pre-fix behaviour: the step dies
                rc = f"raised {type(e).__name__}: {e}"
        rp = os.path.join(d, "BADPYPGX_cpic_recommendations.txt")
        rec = open(rp).read() if os.path.exists(rp) else ""
        check("unreadable pypgx summary: exits 0 with the CPIC report written", rc == 0 and "CYP2C19 -- Poor" in rec, rc)
        check("unreadable pypgx summary: the report says no comparison was made, and an earlier one is removed",
              "Could not read the pypgx summary" in rec and not os.path.exists(comp), rec[-600:])

        # 5. Indeterminate at gene level, a phenotype in the drug annotations
        iann = lambda phen, label="*1/*80", rec=None: {"phenotypes": {"UGT1A1": phen}, "classification": "Strong",
                                                       "drugRecommendation": rec or f"Guidance for {phen}.",
                                             "genotypes": [{"diplotypes": [{"gene": "UGT1A1", "label": label}]}]}
        drugs = lambda *anns: {"CPIC Guideline Annotation": {"atazanavir": {"name": "atazanavir", "guidelines": [
            {"annotations": list(anns)}]}}}
        gene = {"UGT1A1": dict(dip("*1", "*80", "Indeterminate"), relatedDrugs=[{"name": "atazanavir"}])}
        note = "[gene level Indeterminate; phenotype from PharmCAT's drug annotations]"
        rc, rec, rows, _ = run_report(work, "INDET", {"genes": gene, "drugs": drugs(iann("Intermediate Metabolizer"))})
        check("Indeterminate gene, one drug-level phenotype: exits 0", rc == 0, rc)
        check("Indeterminate gene: the table shows the drug-level phenotype with its note, and its status",
              rows == [["UGT1A1", "*1/*80", f"Intermediate Metabolizer {note}", "non-normal"]], rows)
        check("Indeterminate gene: the gene results show the note",
              f"Intermediate Metabolizer {note}" in rec.split("Affected Medications:")[0], rec[:1500])
        check("Indeterminate gene: listed with its drugs, the guidance for the called diplotype",
              f"  UGT1A1 -- Intermediate Metabolizer {note}:" in rec
              and "- atazanavir [CPIC, Strong]: Guidance for Intermediate Metabolizer." in block(rec, "UGT1A1"),
              block(rec, "UGT1A1"))
        check("Indeterminate gene: not in the section of genes without a function phenotype",
              "Called Genes Without a Function Phenotype:" not in rec, rec)
        for label, d in (("no drug annotation", {}),
                         ("two drug annotations that disagree",
                          drugs(iann("Intermediate Metabolizer"), iann("Poor Metabolizer"))),
                         ("a drug annotation that is Indeterminate too", drugs(iann("Indeterminate"))),
                         ("one drug annotation with a phenotype and one Indeterminate",
                          drugs(iann("Intermediate Metabolizer"), iann("Indeterminate")))):
            rc, rec, rows, _ = run_report(work, "INDET" + str(len(label)), dict({"genes": gene}, **({"drugs": d} if d else {})))
            check(f"Indeterminate gene, {label}: stays Indeterminate and unclassified, no note",
                  rc == 0 and rows == [["UGT1A1", "*1/*80", "Indeterminate", "unclassified"]] and note not in rec, rows)
        rc, rec, rows, _ = run_report(work, "INDETNORMAL", {"genes": gene, "drugs": drugs(iann("Normal Metabolizer"))})
        check("Indeterminate gene, drug level Normal Metabolizer: normal, with the note",
              rows == [["UGT1A1", "*1/*80", f"Normal Metabolizer {note}", "normal"]], rows)
        # two possible diplotypes, both Indeterminate, whose drug annotations
        # disagree: the gene is undecided and is not read as its first one
        two = dip("*1", "*80", "Indeterminate")
        two["sourceDiplotypes"] += dip("*80", "*80", "Indeterminate")["sourceDiplotypes"]
        rc, rec, rows, _ = run_report(work, "INDETTWO", {"genes": {"UGT1A1": two}, "drugs": drugs(
            iann("Intermediate Metabolizer", "*1/*80"), iann("Normal Metabolizer", "*80/*80"))})
        check("Indeterminate gene, two possible diplotypes whose drug annotations disagree: stays Indeterminate, no note",
              rc == 0 and rows == [["UGT1A1", "*1/*80 (1 of 2 possible)", "Indeterminate", "unclassified"]]
              and note not in rec, rows)
        # the same phenotype for both, but different guidance: the drug list
        # would show the first diplotype's, so the gene stays undecided
        rc, rec, rows, _ = run_report(work, "INDETTWOREC", {"genes": {"UGT1A1": two}, "drugs": drugs(
            iann("Intermediate Metabolizer", "*1/*80"),
            iann("Intermediate Metabolizer", "*80/*80", "Other guidance."))})
        check("Indeterminate gene, two possible diplotypes with one phenotype but different guidance: stays Indeterminate",
              rc == 0 and rows == [["UGT1A1", "*1/*80 (1 of 2 possible)", "Indeterminate", "unclassified"]]
              and note not in rec, rows)
        # the same phenotype and the same guidance for both: it takes the phenotype
        rc, rec, rows, _ = run_report(work, "INDETTWOSAME", {"genes": {"UGT1A1": two}, "drugs": drugs(
            iann("Intermediate Metabolizer", "*1/*80"), iann("Intermediate Metabolizer", "*80/*80"))})
        check("Indeterminate gene, two possible diplotypes with the same phenotype and guidance: takes it, with the note",
              rows == [["UGT1A1", "*1/*80 (1 of 2 possible)", f"Intermediate Metabolizer {note}", "non-normal"]], rows)

        # 6. label-only differences, phenotypes, long text, single-position fallback
        vk = "rs9923231 reference (C)/rs9923231 variant (T)"
        data = {"genes": {"VKORC1": {"sourceDiplotypes": [{"label": vk, "phenotypes": ["-1639 GA"]}]},
                          "IFNL3": {"sourceDiplotypes": [{"label": "rs12979860 variant (T)/rs12979860 variant (T)",
                                                          "phenotypes": ["unfavorable response genotype"]}]}}}
        rc, rec, _, comp = run_report(work, "LABELS", data, pypgx={
            "VKORC1": ("Reference/rs9923231", "no pypgx phenotype rule"),
            "IFNL3": ("Reference/rs12979860", "Unfavorable Response")})
        check("VKORC1 'rsX reference (C)/rsX variant (T)' against pypgx 'Reference/rsX' is a Match",
              ["VKORC1", vk, "Reference/rs9923231", "Yes", "both", "-1639 GA", "no pypgx phenotype rule"] in comp, comp)
        check("IFNL3 'rsX variant (T)/rsX variant (T)' against pypgx 'Reference/rsX' is still a conflict",
              [r[3] for r in comp if r[0] == "IFNL3"] == ["No"], comp)
        check("the comparison shows PharmCAT's phenotype beside pypgx's 'no pypgx phenotype rule', no Indeterminate",
              all("Indeterminate" not in "\t".join(r) for r in comp), comp)
        with open(os.path.join(work, "LABELS", "comparison.tsv")) as f:
            head = f.readline().rstrip("\n").split("\t")
        check("the comparison header ends with the two phenotype columns",
              head == ["Gene", "PharmCAT_diplotype", "pypgx_diplotype", "Match", "Called_by",
                       "PharmCAT_phenotype", "pypgx_phenotype"], head)

        long = ("Initiate therapy with 50% of the starting dose followed by titration of dose based on "
                "toxicity or therapeutic drug monitoring. " * 4 + "If plasma concentrations are not available, "
                "consider a lower dose of a tricyclic in this patient.")
        data = {"genes": {"CYP2C19": dip("*2", "*2", "Poor Metabolizer")},
                "drugs": {"CPIC Guideline Annotation": {"amitriptyline": {"name": "amitriptyline", "guidelines": [
                    {"annotations": [{"phenotypes": {"CYP2C19": "Poor Metabolizer"}, "classification": "Optional",
                                      "drugRecommendation": long}]}]}}}}
        rc, rec, _, _ = run_report(work, "LONG", data)
        check(f"a recommendation of {len(long)} characters prints whole",
              len(long) > 400 and f"- amitriptyline [CPIC, Optional]: {long}" in rec and "..." not in block(rec, "CYP2C19"),
              block(rec, "CYP2C19"))

        def variants(gene, call, rsid):  # call before rsid: an rsID next to a genotype reads as a leak
            return {"sourceDiplotypes": [{"allele1": {"name": "Unknown"}, "allele2": {"name": "Unknown"},
                                          "label": "Unknown/Unknown", "phenotypes": []}],
                    "variants": [{"position": 1, "dbSnpId": "rs0000001", "call": "G/A"},
                                 {"position": 2, "dbSnpId": rsid, "call": call, "referenceAllele": "C"}]}
        fb = "Single-Position Fallback (not a PharmCAT call):"
        rc, rec, rows, _ = run_report(work, "CYP4F2HET", {"genes": {"CYP4F2": variants("CYP4F2", "C/T", "rs2108622")}})
        tail = rec.split(fb)[-1]
        check("CYP4F2 No Result with C/T at rs2108622: the fallback names the genotype and the *3 rule, labelled",
              rc == 0 and fb in rec
              and "WARNING: PharmCAT has no result for CYP4F2; single-position fallback, not a PharmCAT call: "
                  "rs2108622 reads C/T, one copy of the *3 variant (T, p.V433M)." in tail
              and "Warfarin (CPIC): the CPIC warfarin dosing algorithm" in tail and "5-10%" in tail, rec[-1500:])
        check("CYP4F2 fallback: PharmCAT's own result stays not called in the phenotypes table",
              rows == [["CYP4F2", "Unknown/Unknown", "No Result", "not called"]], rows)
        rc, rec, _, _ = run_report(work, "CYP4F2HOM", {"genes": {"CYP4F2": variants("CYP4F2", "T|T", "rs2108622")}})
        check("CYP4F2 T|T (phased): two copies", "rs2108622 reads T|T, two copies of the *3 variant" in rec, rec[-1200:])
        rc, rec, _, _ = run_report(work, "CYP4F2REF", {"genes": {"CYP4F2": variants("CYP4F2", "C/C", "rs2108622")}})
        check("CYP4F2 C/C: no copy, the rule does not apply, and no WARNING",
              "rs2108622 reads C/C, no copy of the *3 variant (T, p.V433M). The CPIC warfarin rule for this position does not apply."
              in rec and "WARNING: PharmCAT has no result for CYP4F2" not in rec, rec[-1200:])
        rc, rec, _, _ = run_report(work, "CYP4F2MISS", {"genes": {"CYP4F2": variants("CYP4F2", None, "rs2108622")}})
        check("CYP4F2 with no genotype at rs2108622: not callable",
              "CYP4F2 -- PharmCAT has no result; rs2108622 has no genotype in the VCF PharmCAT read: not callable."
              in rec and "WARNING" not in rec.split(fb)[-1], rec[-1200:])
        rc, rec, _, _ = run_report(work, "VKORC1HET", {"genes": {"VKORC1": variants("VKORC1", "C/T", "rs9923231")}})
        check("VKORC1 No Result with C/T at rs9923231: the fallback reads -1639 GA",
              "single-position fallback, not a PharmCAT call: rs9923231 reads C/T, -1639 GA." in rec, rec[-1200:])
        called = {"genes": {"CYP4F2": dict(dip("*1", "*3", "Decreased Function"),
                                           variants=[{"dbSnpId": "rs2108622", "call": "C/T"}])}}
        rc, rec, _, _ = run_report(work, "CYP4F2CALLED", called)
        check("CYP4F2 called by PharmCAT: no fallback section", fb not in rec, rec[-1200:])
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
