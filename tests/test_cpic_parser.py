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
  7. an MT-RNR1 outside call (step 07 makes it from step 20's chrM calls) is
     not checked against step 36's consensus table, which never lists it:
     with a table that passes CYP2D6 only, MT-RNR1 stays called, gets its
     drug list and no stale warning, while an unconfirmed CYP2D6 outside call
     is still held back; `step36-outside` lists the step 36 genes the report
     calls from an outside call and never MT-RNR1.
  8. phenoconversion (opt-in --medications with assets/pgx_inhibitors.tsv):
     every row of the inhibitor table has the shape the parser reads (gene,
     strength, CPIC's factor for that strength, an FDA source, a CPIC source
     for CYP2D6); a strong CYP2D6 inhibitor on the list turns a normal
     metabolizer into 'Poor Metabolizer (phenoconverted by <drug>)' beside the
     unchanged genotype phenotype, a moderate one multiplies the activity
     score by 0.5, a moderate one without an activity score converts nothing,
     a CYP2C19 inhibitor gets a note and no converted phenotype, a line
     naming no drug of the table is listed as not checked, and a list that
     names no drug leaves the phenotypes table byte for byte as without one.
     Control: the same table with the inhibitor's factor set to 1 must not
     give a poor metabolizer.

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


INHIBITORS = os.path.join(REPO, "assets", "pgx_inhibitors.tsv")
FDA_TABLE = ("https://www.fda.gov/drugs/drug-interactions-labeling/"
             "healthcare-professionals-fdas-examples-drugs-interact-cyp-enzymes-and-transporter-systems")


def run_report(work, name, data, pypgx=None, consensus=None, medications=None, inhibitors=INHIBITORS):
    """Run the cpic-report command; returns (exit code, recommendations, phenotypes rows, comparison rows).
    consensus: rows (Gene, Result, Outside_call) of a step 36 table to pass."""
    d = os.path.join(work, name)
    os.makedirs(d)
    rep = os.path.join(d, "report.json")
    with open(rep, "w") as f:
        json.dump(data, f)
    args = ["cpic-report", "--sample", name, "--report", rep, "--outdir", d]
    if consensus is not None:
        cp = os.path.join(d, "consensus.tsv")
        with open(cp, "w") as f:
            f.write("Gene\tResult\tOutside_call\tReason\tEvidence\n")
            for g, res, out in consensus:
                f.write(f"{g}\t{res}\t{out}\ttest\ttest\n")
        args += ["--consensus", cp]
    if medications is not None:
        mp = os.path.join(d, "medications.txt")
        with open(mp, "w") as f:
            f.write(medications)
        args += ["--medications", mp, "--inhibitors", inhibitors]
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


def pheno_report(nm_score="2.0"):
    """A synthetic report: CYP2D6 *1/*2 normal metabolizer (activity score
    nm_score, None for none) and CYP2C19 *1/*1 normal metabolizer."""
    d6 = dip("*1", "*2", "Normal Metabolizer")
    d6["sourceDiplotypes"][0]["activityScore"] = nm_score
    return {"pharmcatVersion": "3.4.0", "genes": {"CYP2D6": d6, "CYP2C19": dip("*1", "*1", "Normal Metabolizer")}}


def section(rec):
    """The phenoconversion section of the recommendations text."""
    part = rec.split("Medications and Phenoconversion (CYP2D6, CYP2C19):", 1)
    return part[1].split("Affected Medications:", 1)[0] if len(part) == 2 else ""


def phenoconversion_checks(work):
    # 8a. the inhibitor table
    with open(INHIBITORS) as f:
        lines = [line for line in f if not line.startswith("#")]
    header = lines[0].rstrip("\n").split("\t")
    check("inhibitors: the header", header == ["Drug", "Gene", "Strength", "Activity_score_factor",
                                                "Class_source", "Rule_source"], header)
    rows = [line.rstrip("\n").split("\t") for line in lines[1:] if line.strip()]
    factor = {("CYP2D6", "strong"): "0", ("CYP2D6", "moderate"): "0.5",
              ("CYP2C19", "strong"): "-", ("CYP2C19", "moderate"): "-"}
    bad = [r for r in rows if len(r) != 6 or (r[1], r[2]) not in factor or r[3] != factor[(r[1], r[2])]
           or r[0] != r[0].strip().lower() or not r[0] or r[4] != FDA_TABLE
           or r[5] != ("https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
                       if r[1] == "CYP2D6" else "-")]
    check("inhibitors: every row has a gene, a strength, CPIC's factor for it and its sources", not bad, bad)
    keys = [(r[0], r[1]) for r in rows]
    check("inhibitors: no drug twice for one gene", len(keys) == len(set(keys)), keys)
    for gene in pgx_parse.PHENOCONVERSION_GENES:
        for st in ("strong", "moderate"):
            check(f"inhibitors: {gene} has {st} inhibitors", any(r[1] == gene and r[2] == st for r in rows))
    check("inhibitors: the parser reads the table", len(pgx_parse.read_inhibitors(INHIBITORS)) == len(rows))

    # 8b. a strong CYP2D6 inhibitor, a CYP2C19 inhibitor, a name not in the table
    rc, rec, rows, _ = run_report(work, "PHENO", pheno_report(),
                                  medications="Paroxetine 20 mg daily\nfluconazole\nnotadrug\n# a comment\n\n")
    got = {r[0]: r for r in rows}
    check("phenoconversion: exits 0", rc == 0, rc)
    check("phenoconversion: strong CYP2D6 inhibitor turns the NM into 'Poor Metabolizer (phenoconverted by paroxetine)'",
          got.get("CYP2D6", [""] * 5)[4:] == ["Poor Metabolizer (phenoconverted by paroxetine)"], got.get("CYP2D6"))
    check("phenoconversion: the genotype phenotype is kept beside it",
          got.get("CYP2D6", [""] * 5)[2:4] == ["Normal Metabolizer", "normal"], got.get("CYP2D6"))
    sec = section(rec)
    check("phenoconversion: a WARNING line with both phenotypes",
          "WARNING: CYP2D6 genotype phenotype Normal Metabolizer; phenoconverted phenotype "
          "Poor Metabolizer (phenoconverted by paroxetine)" in sec, sec)
    check("phenoconversion: CPIC's wording for a strong inhibitor",
          "the CYP2D6 activity score is adjusted to 0 and the predicted phenotype is a poor metabolizer" in sec, sec)
    check("phenoconversion: CYP2C19 inhibitor gets a note and no converted phenotype",
          got.get("CYP2C19", [""] * 5)[4].startswith("not computed: fluconazole")
          and "CYP2C19 genotype phenotype Normal Metabolizer; on fluconazole (strong CYP2C19 inhibitor)" in sec,
          (got.get("CYP2C19"), sec))
    check("phenoconversion: a name not in the table is listed, not dropped",
          "Not in the inhibitor table: notadrug." in sec, sec)
    check("phenoconversion: the docs line says it is an aid for the clinician",
          "an aid for your clinician, not advice" in sec, sec)

    # 8b'. a line that starts with a separator still names its drug
    _, rec, rows, _ = run_report(work, "PHENOSEP", pheno_report(), medications=", paroxetine\n;\n")
    got = {r[0]: r for r in rows}
    check("phenoconversion: ', paroxetine' is read as paroxetine",
          got.get("CYP2D6", [""] * 5)[4:] == ["Poor Metabolizer (phenoconverted by paroxetine)"]
          and "Not in the inhibitor table" not in section(rec), (got.get("CYP2D6"), section(rec)))

    # 8c. moderate: activity score x 0.5; without a score nothing is converted
    _, rec, rows, _ = run_report(work, "PHENOMOD", pheno_report("2.0"), medications="duloxetine\n")
    got = {r[0]: r for r in rows}
    check("phenoconversion: moderate inhibitor, AS 2 x 0.5 = 1 is an intermediate metabolizer",
          got["CYP2D6"][4] == "Intermediate Metabolizer (phenoconverted by duloxetine)"
          and "activity score 2 x 0.5 = 1" in section(rec), (got["CYP2D6"], section(rec)))
    check("phenoconversion: a gene without an inhibitor on the list gets '-'", got["CYP2C19"][4] == "-", got["CYP2C19"])
    _, rec, rows, _ = run_report(work, "PHENONOAS", pheno_report(None), medications="duloxetine\n")
    got = {r[0]: r for r in rows}
    check("phenoconversion: moderate inhibitor without an activity score converts nothing",
          got["CYP2D6"][4] == "not converted: the report gives no activity score", got["CYP2D6"])

    # 8d. a list that names no drug changes nothing
    _, rec0, _, _ = run_report(work, "PHENONONE", pheno_report())
    _, rec1, _, _ = run_report(work, "PHENOEMPTY", pheno_report(), medications="# nothing taken\n\n")
    t0 = open(os.path.join(work, "PHENONONE", "PHENONONE_phenotypes.tsv")).read()
    t1 = open(os.path.join(work, "PHENOEMPTY", "PHENOEMPTY_phenotypes.tsv")).read()
    check("phenoconversion: an empty list leaves the phenotypes table as without one", t0 == t1, (t0, t1))
    check("phenoconversion: an empty list says so and converts nothing",
          "The medication list names no drug: no phenotype is converted." in section(rec1)
          and "phenoconverted" not in rec1.replace("PHENOEMPTY", ""), section(rec1))
    check("phenoconversion: no list, no section", "Medications and Phenoconversion" not in rec0, rec0[:600])

    # 8e. control: the inhibitor's factor set to 1 must not make a poor metabolizer
    broken = os.path.join(work, "inhibitors-factor-1.tsv")
    with open(INHIBITORS) as f, open(broken, "w") as out:
        for line in f:
            out.write(line.replace("paroxetine\tCYP2D6\tstrong\t0\t", "paroxetine\tCYP2D6\tstrong\t1\t"))
    _, _, rows, _ = run_report(work, "PHENOBROKEN", pheno_report(), medications="paroxetine\n", inhibitors=broken)
    got = {r[0]: r for r in rows}
    check("control: with the factor set to 1 the NM is not converted to a poor metabolizer",
          not got["CYP2D6"][4].startswith("Poor Metabolizer"), got["CYP2D6"])
    # 8f. a NaN activity score converts nothing; a table row with a missing cell is reported, not a crash
    _, rec, rows, _ = run_report(work, "PHENONAN", pheno_report("NaN"), medications="duloxetine\n")
    got = {r[0]: r for r in rows}
    check("phenoconversion: a NaN activity score converts nothing",
          got["CYP2D6"][4] == "not converted: the report gives no activity score", got["CYP2D6"])
    short = os.path.join(work, "inhibitors-short-row.tsv")
    with open(short, "w") as out:
        out.write("Drug\tGene\tStrength\tActivity_score_factor\tClass_source\tRule_source\nparoxetine\n")
    rc, rec, rows, _ = run_report(work, "PHENOSHORT", pheno_report(), medications="paroxetine\n", inhibitors=short)
    check("phenoconversion: a table row with no Gene cell is reported as unreadable, the report still written",
          rc == 0 and "Could not read the medication list or the inhibitor table" in section(rec), (rc, section(rec)))
    nofac = os.path.join(work, "inhibitors-no-factor.tsv")
    with open(INHIBITORS) as f, open(nofac, "w") as out:
        for line in f:
            out.write(line.replace("paroxetine\tCYP2D6\tstrong\t0\t", "paroxetine\tCYP2D6\tstrong\t\t"))
    rc, rec, rows, _ = run_report(work, "PHENONOFAC", pheno_report(), medications="paroxetine\n", inhibitors=nofac)
    check("phenoconversion: a CYP2D6 row with an empty factor is reported as unreadable, the report still written",
          rc == 0 and "Could not read the medication list or the inhibitor table" in section(rec), (rc, section(rec)))
    with contextlib.redirect_stderr(io.StringIO()):
        try:
            pgx_parse.main(["cpic-report", "--sample", "X", "--report", "r.json", "--outdir", work,
                            "--medications", "m.txt"])
            code = 0
        except SystemExit as e:
            code = e.code
    check("--medications without --inhibitors is refused", code == 2, code)


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
        check("VKORC1 No Result with C/T at rs9923231: a WARNING reads -1639 GA, with the VKORC1 warfarin rule",
              "WARNING: PharmCAT has no result for VKORC1; single-position fallback, not a PharmCAT call: "
              "rs9923231 reads C/T, -1639 GA." in rec
              and "Warfarin (CPIC): the CPIC warfarin guideline doses by a pharmacogenetic algorithm "
                  "that takes VKORC1 -1639G>A with CYP2C9; each A allele lowers the expected dose." in rec,
              rec[-1200:])
        called = {"genes": {"CYP4F2": dict(dip("*1", "*3", "Decreased Function"),
                                           variants=[{"dbSnpId": "rs2108622", "call": "C/T"}])}}
        rc, rec, _, _ = run_report(work, "CYP4F2CALLED", called)
        check("CYP4F2 called by PharmCAT: no fallback section", fb not in rec, rec[-1200:])

        # 7. MT-RNR1 from step 20's chrM calls next to step 36's consensus
        def outside(gene, a1, a2, phen, drugs):
            g = dip(a1, a2, phen) if a2 else {"sourceDiplotypes": [
                {"allele1": {"name": a1}, "allele2": None, "label": a1, "phenotypes": [phen]}]}
            return dict(g, callSource="OUTSIDE", relatedDrugs=[{"name": n} for n in drugs])
        mito = {"genes": {
            "MT-RNR1": outside("MT-RNR1", "m.1555A>G", None,
                               "increased risk of aminoglycoside-induced hearing loss", ["amikacin", "gentamicin"]),
            "CYP2D6": outside("CYP2D6", "*1", "*4", "Intermediate Metabolizer", ["codeine"])}}
        rc, rec, rows, _ = run_report(work, "MTRNR1CONS", mito, consensus=[("CYP2D6", "*1/*1", "yes")])
        check("MT-RNR1 outside call, consensus without it: exits 0", rc == 0, rc)
        check("MT-RNR1 outside call, consensus without it: stays called and non-normal",
              ["MT-RNR1", "m.1555A>G", "increased risk of aminoglycoside-induced hearing loss", "non-normal"] in rows, rows)
        check("MT-RNR1 outside call, consensus without it: listed with its drugs",
              "  MT-RNR1 -- " in rec and "amikacin" in block(rec, "MT-RNR1") and "gentamicin" in block(rec, "MT-RNR1"),
              block(rec, "MT-RNR1") or rec[:1500])
        check("MT-RNR1 outside call, consensus without it: no stale warning for MT-RNR1",
              "outside call MT-RNR1" not in rec, rec[-1500:])
        check("control: the CYP2D6 outside call the consensus passes as another diplotype is still held back",
              "WARNING: PharmCAT's report has the outside call CYP2D6 *1/*4" in rec
              and "  CYP2D6 -- " not in rec.split("Affected Medications:", 1)[-1].split("Calls From Other Tools", 1)[0],
              rec[-1500:])
        rep = os.path.join(work, "MTRNR1CONS", "report.json")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = pgx_parse.main(["step36-outside", "--report", rep])
        check("step36-outside lists CYP2D6 and not MT-RNR1", rc == 0 and out.getvalue().split() == ["CYP2D6"],
              (rc, out.getvalue()))
        only = os.path.join(work, "MTRNR1ONLY.json")
        with open(only, "w") as f:
            json.dump({"genes": {"MT-RNR1": mito["genes"]["MT-RNR1"]}}, f)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = pgx_parse.main(["step36-outside", "--report", only])
        check("step36-outside lists nothing for a report whose only outside call is MT-RNR1",
              rc == 0 and out.getvalue() == "", (rc, out.getvalue()))
        phenoconversion_checks(work)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
