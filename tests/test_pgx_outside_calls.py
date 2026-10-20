#!/usr/bin/env python3
"""bin/pgx_outside_calls.py: what reaches PharmCAT as an outside call.

CYP2D6 reaches PharmCAT only when pypgx and a second caller (Cyrius) give the
same diplotype and the depth check passed; HLA-A and HLA-B come from T1K.
Each case runs the script on synthetic inputs (tests/fixtures/pgx/) and
compares both files it writes, byte for byte:

  1. disagree      pypgx and Cyrius differ: no outside call, 'indeterminate'
  2. one missing   pypgx only (Cyrius not run): no outside call
  3. agree         both give the same diplotype, depth ok: the call is passed
  4. depth         both agree, but the depth check found multi-mapped reads
  5. no call       Cyrius ran but made no call: one caller only
  6. HLA           T1K's HLA-A and HLA-B, truncated to two fields
  7. HLA edges     one allele (read as homozygous); an allele of quality 0

and what step 27 (bin/pgx_parse.py cpic-report --consensus) makes of it:

  8. CYP2D6 held back: the CPIC report says so, names the drugs it affects,
     and drops the pypgx-only warning it prints without the table
  9. HLA passed on: the CPIC report lists it in the outside-call section and
     marks the gene PharmCAT reports with callSource OUTSIDE
 10. a stale report: PharmCAT's report has outside calls the consensus does
     not pass on (another HLA-B, a CYP2D6 now held back, or no readable
     table at all); they give no drug guidance and the report says to rerun
     step 07

and step 32's summary (pgx_outside_calls.py pypgx-summary) and the probable
CYP2D6 call the reports show:

 11. the star alleles of a diplotype: a copy count (*2x2) is not an allele,
     a tandem (*36+*10) is two
 12. the defining variants of a call looked up in the VCF: a heterozygous
     defining variant at GQ 12 flags the call 'low-confidence: check the
     reads', the same at GQ 50 does not; the lowest GQ, the weakest variant
     with its AD, a variant missing from the VCF or not PASS; *1/*1 has none
 13. pypgx alone with the depth check ok and PASS defining variants: the
     reports' probable call, while the outside-call file stays empty and
     step 36's row stays indeterminate; no probable call with a failed depth
     check, a usable Cyrius call, a low-GQ or non-PASS defining variant, or
     a summary without the columns
 14. the whole subcommand on a fake pypgx and pysam: GSTT1 not callable on
     this reference (not FAILED), a gene without results FAILED, pypgx's CNV
     call for a BAM gene only, the evidence columns per gene

Run: python3 tests/test_pgx_outside_calls.py
"""
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "pgx_outside_calls.py")
FX = os.path.join(REPO, "tests", "fixtures", "pgx")
FAILS = []

HEADER = "Gene\tResult\tOutside_call\tReason\tEvidence\n"
NO_HLA = ("HLA-A\tnot typed\tno\tHLA typing (step 08) did not run\t-\n"
          "HLA-B\tnot typed\tno\tHLA typing (step 08) did not run\t-\n")
MULTIMAPPED = ("CYP2D6 copy number unreliable: reads are multi-mapped "
               "(was this BAM aligned to a reference with ALT contigs?)")


def run(**inputs):
    """(calls file, consensus file) the script writes for these inputs."""
    with tempfile.TemporaryDirectory() as tmp:
        calls, cons = os.path.join(tmp, "calls.tsv"), os.path.join(tmp, "consensus.tsv")
        cmd = [sys.executable, SCRIPT, "--calls", calls, "--consensus", cons]
        for opt, name in inputs.items():
            cmd += ["--" + opt.replace("_", "-"), os.path.join(FX, name)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return f"exit {r.returncode}: {r.stderr.strip()}", ""
        return open(calls).read(), open(cons).read()


def case(desc, want_calls, want_consensus, **inputs):
    calls, cons = run(**inputs)
    for what, got, want in (("outside calls", calls, want_calls), ("consensus", cons, want_consensus)):
        if got == want:
            print(f"PASS {desc}: {what}")
        else:
            print(f"FAIL {desc}: {what}\n--- got\n{got}--- want\n{want}---")
            FAILS.append(f"{desc}: {what}")


# 1. Written first: a naive "take pypgx" passes pypgx's call on and fails here.
case("1 disagree", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tpypgx and Cyrius disagree\tpypgx: *1/*4; Cyrius: *1/*2; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_4.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_ok.tsv")

case("2 one missing", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
     "\tpypgx: *1/*2; Cyrius: not run; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", depth_check="depth_check_ok.tsv")

case("3 agree", "CYP2D6\t*1/*2\n",
     HEADER + NO_HLA +
     "CYP2D6\t*1/*2\tyes\tpypgx and Cyrius agree; the depth check passed"
     "\tpypgx: *1/*2; Cyrius: *1/*2; depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_ok.tsv")

case("4 depth unreliable", "",
     HEADER + NO_HLA +
     f"CYP2D6\tindeterminate\tno\t{MULTIMAPPED}\tpypgx: *1/*2; Cyrius: *1/*2; depth check: unreliable\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv", depth_check="depth_check_multimapped.tsv")

case("4b no depth check", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tthe CYP2D6 depth check did not run"
     "\tpypgx: *1/*2; Cyrius: *1/*2; depth check: not run\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_1_2.tsv")

case("5 Cyrius made no call", "",
     HEADER + NO_HLA +
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius: setup.sh --cyrius)"
     "\tpypgx: *1/*2; Cyrius: no call (None, Not_assigned_to_haplotypes); depth check: ok\n",
     pypgx="pypgx_cyp2d6_1_2.tsv", cyrius="cyrius_none.tsv", depth_check="depth_check_ok.tsv")

case("6 HLA from T1K", "HLA-A\t*02:01/*24:02\nHLA-B\t*07:02/*44:02\n",
     HEADER +
     "HLA-A\t*02:01/*24:02\tyes\tT1K (step 08)\tHLA-A*02:01:01 (quality 60); HLA-A*24:02:01 (quality 60)\n"
     "HLA-B\t*07:02/*44:02\tyes\tT1K (step 08)\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 48)\n"
     "CYP2D6\tindeterminate\tno\tno caller made a call\tpypgx: not run; Cyrius: not run; depth check: not run\n",
     hla="t1k_genotype.tsv")

case("7 HLA one allele, and an allele of quality 0", "HLA-A\t*02:01/*02:01\n",
     HEADER +
     "HLA-A\t*02:01/*02:01\tyes\tT1K (step 08): one allele, read as homozygous\tHLA-A*02:01:01 (quality 60)\n"
     "HLA-B\tindeterminate\tno\tT1K gives an allele quality 0 or below (T1K says to ignore it)"
     "\tHLA-B*07:02:01 (quality 60); HLA-B*44:02:01 (quality 0)\n"
     "CYP2D6\tindeterminate\tno\tno caller made a call\tpypgx: not run; Cyrius: not run; depth check: not run\n",
     hla="t1k_genotype_homozygous_lowqual.tsv")

# --- step 27 with the consensus table ------------------------------------------
PARSE = os.path.join(REPO, "bin", "pgx_parse.py")
REPORTS = os.path.join(REPO, "tests", "fixtures", "pharmcat")


def cpic(report, **inputs):
    """The recommendations text of pgx_parse.py cpic-report, with a consensus
    table from pgx_outside_calls.py when consensus_from is given."""
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [sys.executable, PARSE, "cpic-report", "--sample", "T", "--outdir", tmp,
               "--report", os.path.join(REPORTS, report)]
        if "pypgx" in inputs:
            cmd += ["--pypgx", os.path.join(FX, inputs["pypgx"])]
        if "consensus_from" in inputs:
            cons = os.path.join(tmp, "consensus.tsv")
            gen = [sys.executable, SCRIPT, "--calls", os.path.join(tmp, "calls.tsv"), "--consensus", cons]
            for opt, name in inputs["consensus_from"].items():
                gen += ["--" + opt.replace("_", "-"), os.path.join(FX, name)]
            subprocess.run(gen, check=True, capture_output=True)
            cmd += ["--consensus", cons]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            return f"exit {r.returncode}: {r.stderr}"
        return open(os.path.join(tmp, "T_cpic_recommendations.txt")).read()


def report_has(desc, text, needle, present=True):
    if (needle in text) == present:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: {'missing' if present else 'unexpected'} {needle!r}\n{text}")
        FAILS.append(desc)


# 8. The 3.4.0 report of the fixture: PharmCAT has no CYP2D6 result.
plain = cpic("report-3.4.0.json", pypgx="pypgx_cyp2d6_1_4.tsv")
report_has("8 control: without the table, pypgx's call alone draws a warning", plain,
           "WARNING: PharmCAT has no result for CYP2D6, but pypgx")
held = cpic("report-3.4.0.json", pypgx="pypgx_cyp2d6_1_4.tsv",
            consensus_from={"pypgx": "pypgx_cyp2d6_1_4.tsv", "depth_check": "depth_check_ok.tsv"})
report_has("8 the report has the outside-call section", held, "Calls From Other Tools (outside calls, step 36):")
report_has("8 CYP2D6 is indeterminate and why", held,
           "CYP2D6   indeterminate            not passed to PharmCAT: " + "one caller only")
report_has("8 it names the drugs CYP2D6 affects", held, "No drug guidance is given for CYP2D6 here. Drugs affected by CYP2D6: ")
report_has("8 no pypgx-only warning for the held-back gene", held, "PharmCAT has no result for CYP2D6", False)

# 9. PharmCAT's example report has HLA-B from an outside call (and CYP2D6 and
# MT-RNR1, which this consensus does not pass on: case 10).
passed = cpic("pharmcat-docs-example.json", consensus_from={"hla": "t1k_genotype_docs_example.tsv"})
report_has("9 HLA-B is listed as passed from T1K", passed,
           "HLA-B    *15:02/*57:01            passed to PharmCAT from T1K (step 08); PharmCAT reports it as an outside call.")
report_has("9 the gene line marks PharmCAT's outside call", passed, "[outside call]")
report_has("9 the confirmed HLA-B keeps its drug guidance", passed, "  HLA-B -- ")
report_has("9 no stale warning for the confirmed HLA-B", passed, "outside call HLA-B", False)


def medications(text):
    """The Affected Medications section of a recommendations text."""
    return text.split("Affected Medications:", 1)[-1].split("Calls From Other Tools", 1)[0]


# 10. Outside calls the consensus does not pass on.
plain = cpic("pharmcat-docs-example.json")
report_has("10 control: without the table, the example's CYP2D6 has drug guidance", medications(plain), "  CYP2D6 -- ")
report_has("10 control: without the table, the example's HLA-B has drug guidance", medications(plain), "  HLA-B -- ")
other = cpic("pharmcat-docs-example.json", consensus_from={"hla": "t1k_genotype.tsv"})
report_has("10 HLA-B of another call gets no drug guidance", medications(other), "  HLA-B -- ", False)
report_has("10 the report says the HLA-B outside call is stale", other,
           "WARNING: PharmCAT's report has the outside call HLA-B *15:02/*57:01, which the consensus does not pass on")
report_has("10 the consensus HLA-B is not shown as PharmCAT's", other,
           "HLA-B    *07:02/*44:02            passed to PharmCAT from T1K (step 08); PharmCAT's report does not list it")
held = cpic("pharmcat-docs-example.json", consensus_from={"pypgx": "pypgx_cyp2d6_1_4.tsv", "depth_check": "depth_check_ok.tsv"})
report_has("10 a held-back CYP2D6 gets no drug guidance from an older report", medications(held), "  CYP2D6 -- ", False)
report_has("10 the report says the CYP2D6 outside call is stale", held,
           "WARNING: PharmCAT's report has the outside call CYP2D6 *1/*3, which the consensus does not pass on")
report_has("10 and says how to fix it", held, "Rerun step 07, then this step. No drug guidance is given for CYP2D6 here.")
with tempfile.TemporaryDirectory() as tmp:
    r = subprocess.run([sys.executable, PARSE, "cpic-report", "--sample", "T", "--outdir", tmp,
                        "--report", os.path.join(REPORTS, "pharmcat-docs-example.json"),
                        "--consensus", os.path.join(tmp, "missing.tsv")], capture_output=True, text=True)
    gone = open(os.path.join(tmp, "T_cpic_recommendations.txt")).read() if r.returncode == 0 else f"exit {r.returncode}"
report_has("10 with the table missing, no outside call gets drug guidance", medications(gone), "  HLA-B -- ", False)
report_has("10 with the table missing, the report says to run step 36 again", gone, "Run step 36 again, then steps 07 and 27.")

# --- step 32's summary and the probable CYP2D6 call -------------------------------
sys.path.insert(0, os.path.join(REPO, "bin"))
import pgx_outside_calls as poc  # noqa: E402


def same(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


# 11.
same("11 a copy count is not an allele", poc.star_alleles("*2x2/*4"), ["*2", "*4"])
same("11 a tandem is two alleles, each once", poc.star_alleles("*36+*10/*10"), ["*36", "*10"])
same("11 a no-call has none", poc.star_alleles("Indeterminate"), [])

# 12. CYP2C19 *2 is defined by two variants; *17 by one.
DEFS = {"*1": [], "*2": ["10-94775367-A-G", "10-94781859-G-A"], "*17": ["10-94761900-C-T"]}


def vcf(gq_781859=50, filt="PASS", drop=()):
    recs = {"10-94775367-A-G": {"filter": "PASS", "gq": 48, "ad": "20,19"},
            "10-94781859-G-A": {"filter": filt, "gq": gq_781859, "ad": "26,3"},
            "10-94761900-C-T": {"filter": "PASS", "gq": 60, "ad": "0,31"}}
    return lambda v: None if v in drop else recs.get(v)


ev = poc.defining_evidence("*2/*17", DEFS.get, vcf(gq_781859=12))
same("12 a defining heterozygote at GQ 12 flags the call", ev[4], poc.LOW_GQ)
same("12 the lowest GQ and the weakest variant with its AD", ev[1:3], ["12", "chr10:94781859 G>A GQ 12 AD 26,3"])
same("12 three of three in the VCF, all PASS", (ev[0], ev[3]), ("3/3 in the VCF", "PASS"))
same("12 the same variant at GQ 50 does not flag it", poc.defining_evidence("*2/*17", DEFS.get, vcf())[4], "ok")
same("12 the floor is GQ 20", poc.DEFAULT_GQ_FLOOR, 20)
same("12 GQ 19 is below the floor, GQ 20 is not",
     [poc.defining_evidence("*2/*17", DEFS.get, vcf(gq_781859=g))[4] for g in (19, 20)], [poc.LOW_GQ, "ok"])
same("12 a variant missing from the VCF is named, and the call is not PASS-backed",
     poc.defining_evidence("*2/*17", DEFS.get, vcf(drop=("10-94775367-A-G",)))[3],
     "not in the VCF: chr10:94775367 A>G")
same("12 a RefCall defining variant is not PASS",
     poc.defining_evidence("*2/*17", DEFS.get, vcf(filt="RefCall"))[3], "not PASS: chr10:94781859 G>A RefCall")
same("12 an unfiltered defining variant (FILTER .) is not PASS",
     poc.defining_evidence("*2/*17", DEFS.get, vcf(filt="."))[3], "not PASS: chr10:94781859 G>A no filter")
same("12 *1/*1 rests on no defining variant", poc.defining_evidence("*1/*1", DEFS.get, vcf()),
     ["none", "N/A", "N/A", "none", "N/A"])
same("12 an allele pypgx's table lacks is named, and the call is not PASS-backed",
     poc.defining_evidence("*2/*99", DEFS.get, vcf())[3], "not in pypgx's table: *99")
same("12 a call whose alleles are all unknown to the table rests on no lookup",
     poc.defining_evidence("*98/*99", DEFS.get, vcf()), ["unknown", "N/A", "N/A", "not in pypgx's table: *98, *99", "N/A"])
same("12 a defining variant with no GQ in the VCF says so, it is not 'ok'",
     poc.defining_evidence("*2/*17", DEFS.get, vcf(gq_781859=None))[4], poc.NO_GQ)
same("12 a Float GQ is read like an Integer one",
     poc.defining_evidence("*2/*17", DEFS.get, vcf(gq_781859=12.0))[4], poc.LOW_GQ)
same("12 a no-call has no evidence", poc.defining_evidence("Indeterminate", DEFS.get, vcf()), ["N/A"] * 5)

# 13.
ROW = {"Gene": "CYP2D6", "Diplotype": "*1/*2", "Phenotype": "Normal Metabolizer", "Defining_filter": "PASS",
       "Quality": "ok"}
same("13 pypgx alone, depth ok, PASS defining variants: probable",
     poc.probable_cyp2d6(ROW, None, "ok"),
     {"diplotype": "*1/*2", "phenotype": "Normal Metabolizer", "label": "probable (one caller: pypgx)"})
same("13 Cyrius ran and made no call: still one caller, probable",
     (poc.probable_cyp2d6(ROW, ("None", "Not_assigned_to_haplotypes"), "ok") or {}).get("label"), poc.PROBABLE)
same("13 red control: a failed depth check gives no probable call", poc.probable_cyp2d6(ROW, None, "unreliable"), None)
same("13 a usable Cyrius call: the two callers decide, no probable call",
     poc.probable_cyp2d6(ROW, ("*1/*2", "PASS"), "ok"), None)
same("13 a low-GQ defining variant: no probable call",
     poc.probable_cyp2d6(dict(ROW, Quality=poc.LOW_GQ), None, "ok"), None)
same("13 a defining variant not PASS: no probable call",
     poc.probable_cyp2d6(dict(ROW, Defining_filter="not PASS: chr22:42128945 C>T RefCall"), None, "ok"), None)
same("13 an allele pypgx's table lacks: no probable call",
     poc.probable_cyp2d6(dict(ROW, Defining_filter="not in pypgx's table: *99"), None, "ok"), None)
same("13 a defining variant with no GQ: no probable call",
     poc.probable_cyp2d6(dict(ROW, Quality=poc.NO_GQ), None, "ok"), None)
same("13 a summary without the columns: no probable call",
     poc.probable_cyp2d6({"Diplotype": "*1/*2", "Phenotype": "NM"}, None, "ok"), None)
with tempfile.TemporaryDirectory() as tmp:
    summ = os.path.join(tmp, "pypgx.tsv")
    with open(summ, "w") as f:
        f.write("\t".join(poc.SUMMARY_HEADER) + "\n"
                "CYP2D6\t*1/*2\tNormal Metabolizer\tNormal\tBAM\tcalled\t1/1 in the VCF\t50"
                "\tchr22:42126611 C>G GQ 50 AD 20,22\tPASS\tok\n")
    calls, cons = os.path.join(tmp, "c.tsv"), os.path.join(tmp, "k.tsv")
    subprocess.run([sys.executable, SCRIPT, "--calls", calls, "--consensus", cons, "--pypgx", summ,
                    "--depth-check", os.path.join(FX, "depth_check_ok.tsv")], check=True, capture_output=True)
    same("13 the outside-call file stays empty for a probable call", open(calls).read(), "")
    same("13 step 36's row stays indeterminate, one caller only", open(cons).read().splitlines()[-1].split("\t")[:4],
         ["CYP2D6", "indeterminate", "no", poc.ONE_CALLER])

# 14. The subcommand on a fake pypgx (module and print-data command) and pysam.
FAKE_PYPGX = """
def has_phenotype(g):
    return g != "GSTM1"
DEFS = {("CYP2C19", "*2"): ["10-94775367-A-G", "10-94781859-G-A"], ("CYP2C19", "*1"): [],
        ("GSTM1", "*A"): []}
def list_variants(gene, alleles=None, mode="all", assembly="GRCh37"):
    assert mode == "core" and assembly == "GRCh38"
    return DEFS[(gene, alleles[0])]
"""
FAKE_PRINT_DATA = """#!/usr/bin/env python3
import sys
gene = sys.argv[2].split("/")[-2]
rows = {"CYP2C19": ("*1/*2", "Intermediate Metabolizer", "N/A"), "GSTM1": ("*0/*0", "Indeterminate", "WholeDel2")}
g, p, c = rows[gene]
print("\\tGenotype\\tPhenotype\\tCNV")
print(f"S\\t{g}\\t{p}\\t{c}")
"""
FAKE_PYSAM = """
RECS = {("chr10", 94775367): ("A", ("G",), ["PASS"], 48, (20, 19)),
        ("chr10", 94781859): ("G", ("A",), ["PASS"], 12, (26, 3))}
class _Rec:
    def __init__(self, pos, r):
        self.pos, self.ref, self.alts = pos, r[0], r[1]
        self.filter = {k: None for k in r[2]}
        self.samples = [{"GQ": r[3], "AD": r[4]}]
class VariantFile:
    def __init__(self, path):
        self.header = type("H", (), {"contigs": {"chr10": None, "chr22": None}})()
    def fetch(self, c, s, e):
        return [_Rec(e, RECS[(c, e)])] if (c, e) in RECS else []
"""
with tempfile.TemporaryDirectory() as tmp:
    lib, binp, outbase = os.path.join(tmp, "lib"), os.path.join(tmp, "bin"), os.path.join(tmp, "out")
    for d in (lib, binp, os.path.join(outbase, "CYP2C19"), os.path.join(outbase, "GSTM1")):
        os.makedirs(d)
    open(os.path.join(lib, "pypgx.py"), "w").write(FAKE_PYPGX)
    open(os.path.join(lib, "pysam.py"), "w").write(FAKE_PYSAM)
    open(os.path.join(binp, "pypgx"), "w").write(FAKE_PRINT_DATA)
    os.chmod(os.path.join(binp, "pypgx"), 0o755)
    for g in ("CYP2C19", "GSTM1"):
        open(os.path.join(outbase, g, "results.zip"), "w").close()
    out = os.path.join(tmp, "summary.tsv")
    env = dict(os.environ, PYTHONPATH=lib, PATH=binp + os.pathsep + os.environ["PATH"])
    r = subprocess.run([sys.executable, SCRIPT, "pypgx-summary", "--outbase", outbase, "--bam-genes", "GSTM1 GSTT1 CYP2D6",
                        "--vcf-genes", "CYP2C19", "--vcf", os.path.join(tmp, "x.vcf.gz"), "--out", out,
                        "--not-callable", "GSTT1"], capture_output=True, text=True, env=env)
    rows = {line.split("\t")[0]: line.rstrip("\n").split("\t") for line in open(out)} if r.returncode == 0 else {}
    same("14 the subcommand exits 0", (r.returncode, r.stderr.strip()), (0, ""))
    same("14 the header", rows.get("Gene"), poc.SUMMARY_HEADER)
    same("14 GSTT1 is not callable on this reference, not FAILED", rows.get("GSTT1"),
         ["GSTT1", "N/A", "N/A", "N/A", "BAM", "not callable on this reference"] + ["N/A"] * 5)
    same("14 CYP2D6 without results is FAILED", (rows.get("CYP2D6") or [""] * 6)[1:6:4], ["FAILED", "failed"])
    same("14 GSTM1: pypgx's CNV call, its no-rule phenotype", (rows.get("GSTM1") or [""] * 6)[1:6],
         ["*0/*0", "no pypgx phenotype rule", "WholeDel2", "BAM", "called"])
    same("14 CYP2C19: no CNV call for a VCF gene, the GQ 12 defining variant flags it", rows.get("CYP2C19"),
         ["CYP2C19", "*1/*2", "Intermediate Metabolizer", "N/A", "VCF", "called", "2/2 in the VCF", "12",
          "chr10:94781859 G>A GQ 12 AD 26,3", "PASS", poc.LOW_GQ])
    same("14 the log names the low-confidence call", "check the reads: CYP2C19 *1/*2" in r.stdout, True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
