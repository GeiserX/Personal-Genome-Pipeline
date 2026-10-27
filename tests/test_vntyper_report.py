#!/usr/bin/env python3
"""bin/vntyper_report.py: step 52's reading of VNtyper's MUC1 VNTR result,
and its section in both reports.

On VNtyper run folders made up for each case, in the layout VNtyper 2.0.47
writes (kestrel/kestrel_result.tsv after its ## banner lines,
coverage/coverage_summary.tsv as one header and one data line):
  1. a High_Precision* call, not flagged, at 1,075x: positive, with the
     variant's name, and the sensitivity line of the 100x-or-more band
     (160 of 161 runs);
  2. a Low_Precision call, and a High_Precision call VNtyper flagged, are
     possible, never positive; a call over a VNTR more than half without
     reads is possible;
  3. no call with coverage PASS (148x) or REDUCED (70x) is negative, with the
     band of its depth (160 of 161; 95 of 108); no call with coverage FAIL
     (20x) is not assessed, never negative, with the below-50x band (48 of
     109); no measured depth (NA) says the sensitivity is not known;
  4. VNtyper's non-zero exit code, or a run folder without its files, is not
     assessed; the exit on a BAM with no read over MUC1 says so;
  5. the CLI writes every metric row and exits 0, and a table a reader cannot
     parse exits 2;
  6. both reports: the card and the text section carry the result, the
     variant and the sensitivity line; the list of what was not assessed
     names the MUC1 VNTR when step 52 did not run or gave no reading, and not
     when it read positive or negative; no card without step 52's table;
  7. control: a copy of the script that calls every Kestrel call positive
     fails case 2, so these checks can fail.

Run: python3 tests/test_vntyper_report.py
"""
import importlib.util
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "vntyper_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import vntyper_report as vr  # noqa: E402

FAILS = 0

BANNER = "##VNtyper Kestrel result\n##source=test\n"
# The columns VNtyper 2.0.47 writes for a call, in its order, trimmed.
POS_HEADER = ("Motifs\tVariant\tPOS\tREF\tALT\tEstimated_Depth_AlternateVariant\tEstimated_Depth_Variant_ActiveRegion\t"
              "Depth_Score\tConfidence\tFlag\thaplo_count\tNomenclature\tMotif_sequence\n")
NEG = ("Motif\tVariant\tPOS\tREF\tALT\tMotif_sequence\tEstimated_Depth_AlternateVariant\t"
       "Estimated_Depth_Variant_ActiveRegion\tDepth_Score\tConfidence\n"
       "None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNegative\n")
COV_HEADER = "mean\tmedian\tstdev\tmin\tmax\tregion_length\tuncovered_bases\tpercent_uncovered\tcoverage_qc\n"


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS += 1


def run_dir(root, name, kestrel=None, mean="1075.32", qc="PASS", uncovered="0.00"):
    d = os.path.join(root, name)
    os.makedirs(os.path.join(d, "kestrel"))
    os.makedirs(os.path.join(d, "coverage"))
    if kestrel is not None:
        with open(os.path.join(d, "kestrel", "kestrel_result.tsv"), "w") as f:
            f.write(BANNER + kestrel)
    with open(os.path.join(d, "coverage", "coverage_summary.tsv"), "w") as f:
        f.write(COV_HEADER + f"{mean}\t1000.00\t300.00\t95\t1942\t4500\t0\t{uncovered}\t{qc}\n")
    return d


def call(conf="High_Precision*", flag="Not flagged", name="59dupC"):
    return POS_HEADER + f"2\tInsertion\t67\tG\tGC\t416\t7110\t0.0585\t{conf}\t{flag}\t3\t{name}\tCCCCACC\n"


def cases(mod, t):
    """(description, ok, detail) for cases 1 to 4 run with module mod."""
    out = []
    r = mod.read(run_dir(t, "pos", call()), 0)
    out.append(("1. High_Precision*, not flagged, 1,075x: positive with the name",
                r["status"] == "positive" and r["variant"] == "59dupC" and r["confidence"] == "High_Precision*"
                and r["alt_depth"] == "416", r))
    out.append(("1. the sensitivity line names the depth and the 100x band (160 of 161)",
                "1075x" in r["sensitivity"] and "160 of 161 runs (99%)" in r["sensitivity"], r["sensitivity"]))
    r = mod.read(run_dir(t, "low", call(conf="Low_Precision")), 0)
    out.append(("2. a Low_Precision call is possible, not positive",
                r["status"] == "possible" and "low precision" in r["message"], r))
    r = mod.read(run_dir(t, "flag", call(flag="False_Positive_4bp_Insertion")), 0)
    out.append(("2. a flagged High_Precision call is possible, with the flag",
                r["status"] == "possible" and "False_Positive_4bp_Insertion" in r["message"], r))
    r = mod.read(run_dir(t, "holes", call(), uncovered="62.00", qc="FAIL"), 0)
    out.append(("2. a call over a VNTR 62% without reads is possible",
                r["status"] == "possible" and "62% of the VNTR has no read" in r["message"], r))
    r = mod.read(run_dir(t, "neg", NEG, mean="147.94"), 0)
    out.append(("3. no call at 148x, coverage PASS: negative, 100x band",
                r["status"] == "negative" and r["variant"] == "." and "160 of 161" in r["sensitivity"]
                and r["advice"] == mod.ADVICE_NONE, r))
    r = mod.read(run_dir(t, "red", NEG, mean="70.00", qc="REDUCED"), 0)
    out.append(("3. no call at 70x, coverage REDUCED: negative, 50x-100x band (95 of 108)",
                r["status"] == "negative" and "95 of 108 runs (88%)" in r["sensitivity"]
                and "between 50x and 100x" in r["sensitivity"], r))
    r = mod.read(run_dir(t, "thin", NEG, mean="20.00", qc="FAIL"), 0)
    out.append(("3. no call at 20x, coverage FAIL: not assessed, never negative, below-50x band (48 of 109)",
                r["status"] == "not assessed" and "cannot be ruled out" in r["message"]
                and "48 of 109 runs (44%)" in r["sensitivity"], r))
    r = mod.read(run_dir(t, "na", NEG, mean="NA", qc="NOT_EVALUATED"), 0)
    out.append(("3. no measured depth: not assessed, sensitivity not known",
                r["status"] == "not assessed" and r["sensitivity"].startswith("Sensitivity not known"), r))
    r = mod.read(run_dir(t, "rc", call()), 1)
    out.append(("4. VNtyper's exit code 1: not assessed, whatever it wrote",
                r["status"] == "not assessed" and "exit code 1" in r["message"], r))
    r = mod.read(run_dir(t, "empty"), 1, "ValueError: FASTQ layout 'empty' cannot be consumed without dropping reads.")
    out.append(("4. exit 1 on a BAM with no read over MUC1: not assessed, saying so",
                r["status"] == "not assessed" and "no read in the MUC1 region" in r["message"], r))
    r = mod.read(run_dir(t, "nok"), 0)
    out.append(("4. no kestrel_result.tsv: not assessed, naming the file",
                r["status"] == "not assessed" and "kestrel/kestrel_result.tsv" in r["message"], r))
    return out


def main():
    with tempfile.TemporaryDirectory() as t:
        for desc, ok, detail in cases(vr, os.path.join(t, "real")):
            check(desc, ok, detail)

        # 5. the CLI
        d = run_dir(t, "cli", call())
        out_tsv = os.path.join(t, "cli.tsv")
        p = subprocess.run([sys.executable, SCRIPT, "--run-dir", d, "--out", out_tsv], capture_output=True, text=True)
        rows = [line.rstrip("\n").split("\t") for line in open(out_tsv)] if os.path.exists(out_tsv) else []
        check("5. the CLI exits 0 and writes metric, value and every field once",
              p.returncode == 0 and rows[:1] == [["metric", "value"]]
              and [r[0] for r in rows[1:]] == list(vr.FIELDS) and all(len(r) == 2 for r in rows), (p, rows))
        bad = run_dir(t, "bad", "Confidence\n\"x\n")
        with open(os.path.join(bad, "kestrel", "kestrel_result.tsv"), "wb") as f:
            f.write(b"\xff\xfe\x00bad")
        p = subprocess.run([sys.executable, SCRIPT, "--run-dir", bad, "--out", os.path.join(t, "bad.tsv")],
                           capture_output=True, text=True)
        check("5. an unreadable result file exits 2", p.returncode == 2 and "cannot read" in p.stderr, p)

        # 6. the reports
        import collect_summary
        import render_report

        def sample(name, kestrel, **kw):
            sd = os.path.join(t, "s", name, "S")
            os.makedirs(os.path.join(sd, "vntyper"))
            vr.main(["--run-dir", run_dir(t, "r" + name, kestrel, **kw),
                     "--out", os.path.join(sd, "vntyper", "S_vntyper.tsv")])
            s = collect_summary.collect("S", sd)
            return s, render_report.text_report(s), render_report.html_report(s)

        line = collect_summary.VNTYPER_NOT_ASSESSED
        s, txt, html = sample("pos", call())
        sec = s["sections"].get("vntyper", {})
        check("6. the section is read", sec.get("state") == "ok" and sec["data"]["status"] == "positive", sec)
        check("6. the text section has the result, the variant and the sensitivity line",
              "## MUC1 VNTR Frameshifts (ADTKD-MUC1, VNtyper)" in txt and "Result: positive" in txt
              and "Variant: 59dupC (Kestrel confidence High_Precision*" in txt and "160 of 161 runs" in txt, txt)
        check("6. the HTML card has the result badge and the sensitivity line",
              "<h2>MUC1 VNTR Frameshifts (VNtyper)</h2>" in html and 'badge-red">positive<' in html
              and "160 of 161 runs" in html, html[-3000:])
        check("6. a positive reading takes the MUC1 VNTR off the not-assessed list", line not in s["not_assessed"],
              s["not_assessed"])
        s, txt, html = sample("neg", NEG, mean="147.94")
        check("6. a negative reading takes it off too, and the card says negative",
              line not in s["not_assessed"] and 'badge-green">negative<' in html, s["not_assessed"])
        s, txt, html = sample("thin", NEG, mean="20.00", qc="FAIL")
        check("6. not assessed keeps the MUC1 VNTR on the not-assessed list", line in s["not_assessed"],
              s["not_assessed"])
        os.makedirs(os.path.join(t, "none", "S"))
        s = collect_summary.collect("S", os.path.join(t, "none", "S"))
        check("6. without step 52 the MUC1 VNTR is listed as not assessed, with no card and no text section",
              line in s["not_assessed"] and "MUC1 VNTR Frameshifts" not in render_report.html_report(s)
              and "## MUC1 VNTR" not in render_report.text_report(s), s["not_assessed"])
        check("6. and both reports list it under what was not assessed",
              line in render_report.text_report(s) and render_report.E(line) in render_report.html_report(s))

        # 7. control: a copy that calls every Kestrel call positive
        src = open(SCRIPT).read()
        cond = 'conf.startswith("High_Precision") and not flagged'
        check("7. control: the condition to break is in the script", cond in src)
        broken = os.path.join(t, "vntyper_report_broken.py")
        with open(broken, "w") as f:
            f.write(src.replace(cond, "True"))
        spec = importlib.util.spec_from_file_location("vntyper_report_broken", broken)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        failed = [desc for desc, ok, _ in cases(mod, os.path.join(t, "broken")) if not ok]
        check("7. control: the broken copy fails the possible-not-positive cases",
              any(x.startswith("2.") for x in failed), failed)

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
