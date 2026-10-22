#!/usr/bin/env python3
"""bin/pgx_outside_calls.py with OptiType, the opt-in second HLA caller.

With HLA_QC=true (bash) or --hla_qc (Nextflow) step 08 also runs OptiType on
the HLA-A, -B and -C reads, and step 36 passes HLA-A or HLA-B to PharmCAT
only when T1K and OptiType give the same two alleles at two fields. Each case
runs the script on synthetic files:

  1. agree          both callers agree: HLA-A and HLA-B reach PharmCAT
  2. disagree       OptiType differs at HLA-B: HLA-B stays 'indeterminate',
                    both calls are printed, it is not in the outside-call
                    file; HLA-A still passes
  3. homozygous     T1K's one allele and OptiType's A1 == A2 agree
  4. no call        a file step 08 wrote instead of a result ('# OptiType
                    not run: ...'): both genes withheld
  5. asked, missing --require-optitype without --optitype: both withheld
  6. T1K alone      no OptiType asked for: T1K's call passes, as before
  7. quality 0      T1K's quality-0 rule still wins over an agreement
  8. concordance    the hla-concordance subcommand: HLA-A, -B, -C, yes/no
  9. red control    the agreement test disabled in a copy of the script:
                    the disagreeing HLA-B reaches the outside-call file and
                    case 2's check reports it, so case 2 can fail

Run: python3 tests/test_hla_second_caller.py
"""
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "pgx_outside_calls.py")
FAILS = []

T1K = ("HLA-A\t2\tHLA-A*01:01:01\t110.5\t60\tHLA-A*26:01:01\t95.2\t60\t\n"
       "HLA-B\t2\tHLA-B*57:01:01\t101.0\t60\tHLA-B*38:01:01\t88.7\t48\t\n"
       "HLA-C\t1\tHLA-C*06:02:01\t150.3\t60\t.\t0\t-1\t\n")
T1K_HOM = ("HLA-A\t1\tHLA-A*02:01:01\t200.0\t60\t.\t0\t-1\t\n"
           "HLA-B\t2\tHLA-B*57:01:01\t101.0\t60\tHLA-B*38:01:01\t88.7\t48\t\n")
T1K_Q0 = ("HLA-A\t2\tHLA-A*01:01:01\t110.5\t60\tHLA-A*26:01:01\t95.2\t60\t\n"
          "HLA-B\t2\tHLA-B*57:01:01\t101.0\t0\tHLA-B*38:01:01\t88.7\t48\t\n")
OT_HEAD = "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
OT_AGREE = OT_HEAD + "0\tA*26:01\tA*01:01\tB*38:01\tB*57:01\tC*06:02\tC*06:02\t812\t790.4\n"
OT_DISAGREE = OT_HEAD + "0\tA*01:01\tA*26:01\tB*58:01\tB*38:01\tC*06:02\tC*07:01\t812\t790.4\n"
OT_HOM = OT_HEAD + "0\tA*02:01\tA*02:01\tB*57:01\tB*38:01\t\t\t500\t480.0\n"
OT_NONE = "# OptiType not run: no reads at HLA-A, HLA-B or HLA-C in this BAM\n"


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ("" if ok else f": {detail}"))
    if not ok:
        FAILS.append(name)


def run(tmp, t1k, optitype=None, require=False, script=SCRIPT):
    """(outside-call text, {gene: consensus row}) of one run."""
    hla = os.path.join(tmp, "t1k.tsv")
    with open(hla, "w") as f:
        f.write(t1k)
    calls, cons = os.path.join(tmp, "calls.tsv"), os.path.join(tmp, "cons.tsv")
    cmd = [sys.executable, script, "--calls", calls, "--consensus", cons, "--hla", hla]
    if optitype is not None:
        ot = os.path.join(tmp, "optitype.tsv")
        with open(ot, "w") as f:
            f.write(optitype)
        cmd += ["--optitype", ot]
    if require:
        cmd.append("--require-optitype")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return None, {"error": r.stderr}
    rows = {}
    with open(cons) as f:
        for line in f.read().splitlines()[1:]:
            c = line.split("\t")
            rows[c[0]] = c
    with open(calls) as f:
        return f.read(), rows


def passed(calls, gene, dip):
    return calls is not None and f"{gene}\t{dip}\n" in calls


def withheld(calls, rows, gene, reason_re):
    r = rows.get(gene) or []
    return (calls is not None and not re.search(rf"^{gene}\t", calls, re.M)
            and r[1:3] == ["indeterminate", "no"] and re.search(reason_re, r[3]) is not None)


def disagree_checks(tag, tmp, script=SCRIPT):
    """Case 2's checks; True when every one passed."""
    calls, rows = run(tmp, T1K, OT_DISAGREE, script=script)
    before = len(FAILS)
    check(f"{tag} HLA-B withheld when OptiType disagrees",
          withheld(calls, rows, "HLA-B", r"^T1K and OptiType disagree at two fields$"), f"{calls!r} {rows.get('HLA-B')}")
    check(f"{tag} both HLA-B calls are printed",
          (rows.get("HLA-B") or [""] * 5)[4] == "T1K: HLA-B*57:01:01 (quality 60); HLA-B*38:01:01 (quality 48); "
                                                  "OptiType: *58:01/*38:01", rows.get("HLA-B"))
    check(f"{tag} HLA-A, where they agree, still passes", passed(calls, "HLA-A", "*01:01/*26:01"), calls)
    return len(FAILS) == before


with tempfile.TemporaryDirectory() as tmp:
    # 1. agree
    calls, rows = run(tmp, T1K, OT_AGREE)
    check("1 HLA-A passes when both agree (either order)", passed(calls, "HLA-A", "*01:01/*26:01"), calls)
    check("1 HLA-B passes when both agree", passed(calls, "HLA-B", "*57:01/*38:01"), calls)
    check("1 the reason names both callers",
          (rows.get("HLA-B") or [""] * 4)[3] == "T1K and OptiType (step 08), which agree at two fields", rows.get("HLA-B"))
    check("1 HLA-C is never an outside call", calls is not None and "HLA-C" not in calls, calls)

    # 2. disagree
    disagree_checks("2", tmp)

    # 3. homozygous
    calls, rows = run(tmp, T1K_HOM, OT_HOM)
    check("3 T1K's one allele agrees with OptiType's A1 == A2", passed(calls, "HLA-A", "*02:01/*02:01"), calls)
    check("3 the homozygous HLA-A row says both callers",
          (rows.get("HLA-A") or [""] * 4)[3] == "T1K and OptiType (step 08), which agree at two fields", rows.get("HLA-A"))

    # 4. OptiType did not run (the file step 08 writes instead)
    calls, rows = run(tmp, T1K, OT_NONE)
    for g in ("HLA-A", "HLA-B"):
        check(f"4 {g} withheld: OptiType made no call",
              withheld(calls, rows, g, r"^OptiType, the second caller, made no call$"), f"{calls!r} {rows.get(g)}")

    # 5. asked for, no result
    calls, rows = run(tmp, T1K, None, require=True)
    for g in ("HLA-A", "HLA-B"):
        check(f"5 {g} withheld: OptiType asked for, no result",
              withheld(calls, rows, g, r"^the second caller was asked for \(HLA_QC\)"), f"{calls!r} {rows.get(g)}")

    # 6. T1K alone, as before the second caller
    calls, rows = run(tmp, T1K)
    check("6 without OptiType, T1K's HLA-B passes", passed(calls, "HLA-B", "*57:01/*38:01"), calls)
    check("6 the reason is T1K's", (rows.get("HLA-B") or [""] * 4)[3] == "T1K (step 08)", rows.get("HLA-B"))

    # 7. T1K's quality 0 still withholds
    calls, rows = run(tmp, T1K_Q0, OT_AGREE)
    check("7 a T1K quality-0 allele withholds HLA-B even when OptiType agrees",
          withheld(calls, rows, "HLA-B", r"quality 0 or below"), f"{calls!r} {rows.get('HLA-B')}")

    # 8. the concordance table of step 08
    t1k, ot, out = (os.path.join(tmp, n) for n in ("c_t1k.tsv", "c_ot.tsv", "c_out.tsv"))
    with open(t1k, "w") as f:
        f.write(T1K)
    with open(ot, "w") as f:
        f.write(OT_DISAGREE)
    r = subprocess.run([sys.executable, SCRIPT, "hla-concordance", "--t1k", t1k, "--optitype", ot, "--out", out],
                       capture_output=True, text=True)
    want = ("Gene\tT1K\tOptiType\tAgree\tNote\n"
            "HLA-A\t*01:01/*26:01\t*01:01/*26:01\tyes\t-\n"
            "HLA-B\t*57:01/*38:01\t*58:01/*38:01\tno\tthe callers disagree at two fields; step 36 does not pass it to PharmCAT\n"
            "HLA-C\t*06:02/*06:02\t*06:02/*07:01\tno\tthe callers disagree at two fields\n")
    got = open(out).read() if r.returncode == 0 and os.path.exists(out) else r.stderr
    check("8 the concordance table, byte for byte", got == want, f"\n{got}")

    # 9. red control: case 2 must fail on a script whose agreement test is off
    src = open(SCRIPT).read()
    broken = src.replace("if not same_pair(names, ot):", "if False:")
    check("9 the mutation applies", broken != src)
    bad = os.path.join(tmp, "pgx_outside_calls_broken.py")
    with open(bad, "w") as f:
        f.write(broken)
    saved = list(FAILS)
    caught = not disagree_checks("9 (control, must FAIL)", tmp, script=bad)
    FAILS[:] = saved
    check("9 with the agreement rule disabled, case 2 rejects the disagreeing HLA-B", caught)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
