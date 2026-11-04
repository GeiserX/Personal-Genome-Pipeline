#!/usr/bin/env python3
"""bin/locityper_report.py: step 53's table from what Locityper genotype
wrote, and its section in both reports.

On Locityper output folders made up for each case, in the layout Locityper
1.7.14 writes (gt/loci/<locus>/res.json.gz, a `success` file beside each
finished locus):
  1. a call with no warning, weighted distance 0.2 and 12 of 900 reads
     unexplained is called and passes the authors' filter, with both
     haplotypes in their columns;
  2. the filter fails on a weighted distance above 30, on 1000 or more and
     20% or more unexplained reads, on a Locityper warning, and on a missing
     weighted distance (a database without `locityper augment`); 1500 of
     10000 unexplained still passes (under 20%);
  3. a locus without `success` or without res.json.gz is not called, with
     Locityper's exit code in the reason when it stopped; a genotype without
     two haplotypes is not called; an unreadable res.json.gz is not called;
  4. the CLI writes one row per BED locus in the BED's order with every
     column, skips '#' lines, and exits 1 on a BED it cannot read;
  5. both reports: the card and the text section carry each locus's call,
     filter and the primary caller's call; no card without step 53's table;
  6. control: a copy of the script whose filter passes everything fails
     case 2, so these checks can fail;
  7. scripts/run-all.sh, under the fakes of scripts/ci/fake-docker (no
     Docker, no Nextflow): a default run neither lists step 53 nor passes
     locityper to Nextflow; TOOLS naming locityper with the database
     installed lists it as running and passes --tools locityper,
     --locityper_data and --locityper_bg_region; without the database it is
     skipped, "data not installed". Control: the opted-in checks fail on the
     default run's output.

Run: python3 tests/test_locityper_report.py
"""
import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "locityper_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import locityper_report as lr  # noqa: E402

FAILS = 0
BED = ("# a comment line\n"
       "chr6\t32652603\t32669077\tHLA-DQB1\n"
       "chr22\t42122284\t42152960\tCYP2D6\n"
       "chr1\t161506761\t161681047\tFCGR\n")


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}" + ("" if ok else f" -- {str(detail)[:800]}"))
    if not ok:
        FAILS += 1


def res(genotype="HG00733.1,NA19240.2", quality=27.8, total=900, unexpl=12, wd=0.2, warnings=None):
    r = {"total_reads": total, "quality": quality, "unexpl_reads": unexpl, "genotype": genotype,
         "options": 120, "dist_type": "edit"}
    if wd is not None:
        r["weight_dist"] = wd
    if warnings is not None:
        r["warnings"] = warnings
    return r


def gt_dir(root, locus, result=None, success=True, raw=None):
    d = os.path.join(root, "loci", locus)
    os.makedirs(d, exist_ok=True)
    if result is not None or raw is not None:
        with gzip.open(os.path.join(d, "res.json.gz"), "wb") as f:
            f.write(raw if raw is not None else json.dumps(result).encode())
    if success:
        open(os.path.join(d, "success"), "w").close()
    return root


def cases(mod, t):
    out = []
    n = [0]

    def row(result=None, success=True, raw=None, exit_code=0):
        n[0] += 1
        root = os.path.join(t, f"c{n[0]}")
        gt_dir(root, "CYP2D6", result, success, raw)
        return mod.row_for("CYP2D6", "chr22:42122285-42152960", root, exit_code)

    r = row(res())
    out.append(("1. a clean call is called and passes, both haplotypes in their columns",
                r["status"] == "called" and r["filter"] == "pass" and r["haplotype_1"] == "HG00733.1"
                and r["haplotype_2"] == "NA19240.2" and r["quality"] == "27.8" and r["weighted_distance"] == "0.200"
                and r["reason"] == ".", r))
    r = row(res(wd=31.5))
    out.append(("2. a weighted distance above 30 fails, saying so",
                r["filter"] == "fail" and "weighted distance 31.50" in r["reason"], r))
    r = row(res(total=3000, unexpl=1200))
    out.append(("2. 1200 of 3000 unexplained reads fails", r["filter"] == "fail" and "1200 of 3000" in r["reason"], r))
    r = row(res(total=10000, unexpl=1500))
    out.append(("2. 1500 of 10000 unexplained reads passes (under 20%)", r["filter"] == "pass", r))
    r = row(res(warnings=["UnexplainedReads"]))
    out.append(("2. a Locityper warning fails, naming it",
                r["filter"] == "fail" and r["warnings"] == "UnexplainedReads" and "UnexplainedReads" in r["reason"], r))
    r = row(res(wd=None))
    out.append(("2. no weighted distance fails, pointing at locityper augment",
                r["filter"] == "fail" and "augment" in r["reason"] and r["weighted_distance"] == ".", r))
    r = row(res(), success=False, exit_code=101)
    out.append(("3. no success file: not called, with Locityper's exit code",
                r["status"] == "not called" and "exit code 101" in r["reason"] and r["genotype"] == "."
                and r["filter"] == ".", r))
    r = row(None)
    out.append(("3. no res.json.gz: not called", r["status"] == "not called" and "no result" in r["reason"], r))
    r = row(res(genotype="HG00733.1"))
    out.append(("3. one haplotype only: not called", r["status"] == "not called" and "no diploid" in r["reason"], r))
    r = row(raw=b"{not json")
    out.append(("3. an unreadable res.json.gz: not called", r["status"] == "not called"
                and "cannot be read" in r["reason"], r))
    return out


def run_all(w, tools, done=True):
    """scripts/run-all.sh on a placeholder BAM and VCF under the fake docker
    and nextflow; returns its output and the arguments it gave nextflow."""
    import shutil
    g = os.path.join(w, "g")
    shutil.rmtree(w, ignore_errors=True)
    for f in ("s1/aligned/s1_sorted.bam", "s1/aligned/s1_sorted.bam.bai", "s1/vcf/s1.vcf.gz", "s1/vcf/s1.vcf.gz.tbi",
              "reference/GRCh38_no_alt_analysis_set.fasta"):
        os.makedirs(os.path.dirname(os.path.join(g, f)), exist_ok=True)
        with open(os.path.join(g, f), "w") as fh:
            fh.write("placeholder\n")
    if done:
        d = os.path.join(g, "reference", "locityper-" + lr_panel_version())
        os.makedirs(d)
        with open(os.path.join(d, "DONE"), "w") as fh:
            fh.write("ok\n")
    env = {"PATH": os.path.join(REPO, "scripts", "ci", "fake-docker") + os.pathsep + os.environ["PATH"], "HOME": w,
           "GENOME_DIR": g, "SKIP_VALIDATION": "true", "FAKE_DOCKER_LOG": os.path.join(w, "docker.log"),
           "LOCITYPER_BG_REGION": "chr20:10000001-10500000"}
    if tools:
        env["TOOLS"] = tools
    p = subprocess.run(["bash", os.path.join(REPO, "scripts", "run-all.sh"), "s1", "male"], env=env, stdin=subprocess.DEVNULL,
                       capture_output=True, text=True)
    log = open(os.path.join(w, "docker.log")).read() if os.path.exists(os.path.join(w, "docker.log")) else ""
    nf = [x for x in log.splitlines() if x.startswith("nextflow ::") and "/main.nf " in x]
    return p.stdout + p.stderr, (nf[-1] if nf else "")


def lr_panel_version():
    for line in open(os.path.join(REPO, "versions.env")):
        if line.startswith("LOCITYPER_PANEL_VERSION="):
            return line.split("=", 1)[1].split("#")[0].strip().strip('"')
    raise SystemExit("LOCITYPER_PANEL_VERSION is not in versions.env")


def opted_in(out, nf):
    import re
    return [("7. TOOLS=pharmcat,locityper: run-all.sh lists step 53 as running",
             re.search(r"^  53 Hard loci \(Locityper\) +runs$", out, re.M) is not None, out[-1500:]),
            ("7. and passes locityper in --tools, with --locityper_data and --locityper_bg_region",
             re.search(r"--tools [^ ]*locityper", nf) is not None and "--locityper_data " in nf
             and "--locityper_bg_region chr20:10000001-10500000" in nf, nf)]


def plan_checks(w):
    import re
    out_off, nf_off = run_all(os.path.join(w, "default"), "")
    check("7. a default run does not list step 53", "53 Hard loci" not in out_off, out_off[-1500:])
    check("7. and passes no locityper to Nextflow", bool(nf_off) and re.search(r"--tools [^ ]*locityper", nf_off) is None
          and "--locityper_" not in nf_off, nf_off)
    out_on, nf_on = run_all(os.path.join(w, "on"), "pharmcat,locityper")
    for desc, ok, detail in opted_in(out_on, nf_on):
        check(desc, ok, detail)
    out, nf = run_all(os.path.join(w, "nodata"), "pharmcat,locityper", done=False)
    check("7. without the database step 53 is skipped, data not installed",
          "53 Hard loci (Locityper)     skipped    (data not installed: reference/locityper-" in out
          and "--locityper_data" not in nf, out[-1500:])
    failed = [d for d, ok, _ in opted_in(out_off, nf_off) if not ok]
    check("7. control: the opted-in checks fail on the default run's output", len(failed) == 2, failed)


def main():
    with tempfile.TemporaryDirectory() as t:
        for desc, ok, detail in cases(lr, os.path.join(t, "real")):
            check(desc, ok, detail)

        # 4. the CLI
        bed = os.path.join(t, "loci.bed")
        with open(bed, "w") as f:
            f.write(BED)
        root = os.path.join(t, "cli")
        gt_dir(root, "CYP2D6", res())
        gt_dir(root, "HLA-DQB1", res(wd=50))
        out_tsv = os.path.join(t, "cli.tsv")
        p = subprocess.run([sys.executable, SCRIPT, "--gt-dir", root, "--loci", bed, "--exit-code", "0",
                            "--out", out_tsv], capture_output=True, text=True)
        rows = [line.rstrip("\n").split("\t") for line in open(out_tsv)] if os.path.exists(out_tsv) else []
        check("4. the CLI exits 0, one row per BED locus in its order, every column",
              p.returncode == 0 and rows[:1] == [lr.COLUMNS] and [r[0] for r in rows[1:]] == ["HLA-DQB1", "CYP2D6", "FCGR"]
              and all(len(r) == len(lr.COLUMNS) for r in rows) and rows[1][1] == "chr6:32652604-32669077"
              and [r[2] for r in rows[1:]] == ["called", "called", "not called"]
              and "2 of 3 loci called, 1 pass" in p.stdout, (p, rows))
        badbed = os.path.join(t, "bad.bed")
        with open(badbed, "w") as f:
            f.write("chr6\tx\t10\tA\n")
        p = subprocess.run([sys.executable, SCRIPT, "--gt-dir", root, "--loci", badbed, "--out",
                            os.path.join(t, "bad.tsv")], capture_output=True, text=True)
        check("4. a BED it cannot read exits 1", p.returncode == 1 and "ERROR" in p.stderr, p)

        # 5. the reports
        import collect_summary
        import render_report

        sd = os.path.join(t, "s", "S")
        os.makedirs(os.path.join(sd, "locityper"))
        lr.main(["--gt-dir", root, "--loci", bed, "--out", os.path.join(sd, "locityper", "S_locityper.tsv")])
        s = collect_summary.collect("S", sd)
        txt, html = render_report.text_report(s), render_report.html_report(s)
        sec = s["sections"].get("locityper", {})
        check("5. the section is read, with the counts",
              sec.get("state") == "ok" and sec["data"]["called"] == 2 and sec["data"]["passed"] == 1, sec)
        check("5. each locus carries the primary caller's call",
              [x.get("primary", "")[:5] for x in sec["data"]["loci"]] == ["T1K (", "pypgx", "no ot"], sec["data"])
        check("5. the text section has the calls, the filter and the primary caller",
              "## Locityper Second Opinion at Hard Loci (step 53)" in txt and "HG00733.1 + NA19240.2" in txt
              and "filter fail: weighted distance 50.00" in txt and "primary:   no other step calls FCGR" in txt
              and "not called: Locityper wrote no result" in txt, txt)
        check("5. the HTML card has the table and the filter badges",
              "Locityper Second Opinion at Hard Loci</h2>" in html and 'badge-green">pass<' in html
              and 'badge-yellow">fail<' in html and "<th>Primary caller</th>" in html, html[-3000:])
        os.makedirs(os.path.join(t, "none", "S"))
        s = collect_summary.collect("S", os.path.join(t, "none", "S"))
        check("5. without step 53 there is no card and no text section",
              "Locityper Second Opinion" not in render_report.html_report(s)
              and "Locityper Second Opinion" not in render_report.text_report(s))

        # 6. control: a copy whose filter passes everything
        src = open(SCRIPT).read()
        cond = 'return ("fail" if why else "pass"), "; ".join(why)'
        check("6. control: the line to break is in the script", cond in src)
        broken = os.path.join(t, "locityper_report_broken.py")
        with open(broken, "w") as f:
            f.write(src.replace(cond, 'return "pass", ""'))
        spec = importlib.util.spec_from_file_location("locityper_report_broken", broken)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        failed = [desc for desc, ok, _ in cases(mod, os.path.join(t, "broken")) if not ok]
        check("6. control: the broken copy fails the filter cases", any(x.startswith("2.") for x in failed), failed)

        # 7. run-all.sh plans step 53 only when TOOLS names it
        plan_checks(os.path.join(t, "plan"))

    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
