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
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius, --tools cyrius)"
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
     "CYP2D6\tindeterminate\tno\tone caller only: a second caller must agree (Cyrius, --tools cyrius)"
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

# 9. PharmCAT's example report has HLA-B from an outside call.
passed = cpic("pharmcat-docs-example.json", consensus_from={"hla": "t1k_genotype.tsv"})
report_has("9 HLA-B is listed as passed from T1K", passed,
           "HLA-B    *07:02/*44:02            passed to PharmCAT from T1K (step 08); PharmCAT reports it as an outside call.")
report_has("9 the gene line marks PharmCAT's outside call", passed, "[outside call]")

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
