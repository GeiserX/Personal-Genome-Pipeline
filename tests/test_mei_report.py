#!/usr/bin/env python3
"""bin/mei_report.py: the Alu, L1 and SVA insertions MEGAnE calls, and the known disease insertions (step 50).

Each case writes a synthetic VCF in the format MEGAnE writes (the header and
INFO fields of its MEI_final_*_genotyped.vcf), the catalog of
assets/mei_sites.tsv and a samtools bedcov -c file with chosen depths, runs
the script, and checks what it writes:

  1. an SVA planted inside the FKTN window is seen, with its position,
     subfamily and copies; the same VCF without it says not seen; the same
     SVA just outside the window is not seen; an Alu at the FKTN site is not
     the catalog's SVA
  2. too few reads: a window under 10x mean depth, or a contig the BAM does
     not have (no bedcov row), is not callable; a call is still seen at low
     depth
  3. the calls table keeps Alu, L1 and SVA with MEPRED=PASS and one or two
     copies, with MEGAnE's FILTER and copy-number confidence, in karyotype
     order; it leaves out MEPRED=FAILED, zero copies and other classes (DNA
     transposons, LTRs); the subfamily comes from MEI, else from a breakpoint
  4. damaged input stops the script: a VCF without its #CHROM line, a VCF
     that is not MEGAnE's (no SVTYPE or MEPRED), a depth file that is not
     samtools bedcov -c output, a catalog without its columns
  5. collect_summary.py's mei section and render_report.py's text block and
     HTML card read the tables; without them the summary's "not assessed"
     list names mobile element insertions, with them it does not
  6. the real catalog: every row is an SVA, its window is the site -/+ 500,
     and the script reads it

Run: python3 tests/test_mei_report.py
"""
import csv
import os
import subprocess
import sys
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "mei_report.py")
CATALOG = os.path.join(REPO, "assets", "mei_sites.tsv")
sys.path.insert(0, os.path.join(REPO, "bin"))

FAILS = []
HEAD = """##fileformat=VCFv4.1
##source=MEI search version "v1.0.1 2022/04/12"
##ALT=<ID=INS:ME,Description="Insertion of mobile element">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MEPRED,Number=1,Type=String,Description="ME prediction status">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
"""
SITES = """# a test catalog
gene\tchrom\tstart\tend\tsite\telement\tinheritance\tcondition\tclinvar\tpmid
FKTN\tchr9\t105639156\t105640156\t105639656\tSVA\tautosomal recessive\tFukuyama congenital muscular dystrophy\t3199\t9690476
TAF1\tchrX\t71439993\t71440993\t71440493\tSVA\tX-linked\tX-linked dystonia-parkinsonism\t9862\t29474918
"""


def rec(chrom, pos, svtype, mei="MEI=SVA_F,1/1500,+/+", filt="PASS", mepred="PASS", cn="1", conf="H"):
    """One record as MEGAnE writes it."""
    info = f"SVTYPE={svtype};MEPRED={mepred};HOMLEN=9;MEI_rpos={pos + 8};{mei};CN_conf={conf}"
    return f"{chrom}\t{pos}\t.\tN\t<INS:ME>\t.\t{filt}\t{info}\tCN\t{cn}\n"


FKTN_SVA = rec("chr9", 105639660, "Retroposon/SVA")
# chrom, bedcov start (0-based), end, gene, summed depth over the 1001 bases, reads
DEPTH_30X = "chr9\t105639155\t105640156\tFKTN\t30030\t200\nchrX\t71439992\t71440993\tTAF1\t30030\t200\n"


def check(name, got, want):
    if got == want:
        print(f"PASS {name}")
    else:
        print(f"FAIL {name}: got {got!r}, want {want!r}")
        FAILS.append(name)


def run(records, depth=DEPTH_30X, sites=SITES, vcf=None):
    """(calls, sites) as lists of dicts, or the error text when the script fails."""
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/m.vcf", "w") as f:
            f.write(vcf if vcf is not None else HEAD + "".join(records))
        with open(f"{d}/sites.tsv", "w") as f:
            f.write(sites)
        with open(f"{d}/depth.tsv", "w") as f:
            f.write(depth)
        p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--vcf", f"{d}/m.vcf", "--sites",
                            f"{d}/sites.tsv", "--depth", f"{d}/depth.tsv", "--out-calls", f"{d}/c.tsv",
                            "--out-sites", f"{d}/s.tsv"], capture_output=True, text=True)
        if p.returncode:
            return p.stderr
        return [list(csv.DictReader(open(f"{d}/{n}.tsv"), delimiter="\t")) for n in ("c", "s")]


def status(result, gene):
    return next(r for r in result[1] if r["gene"] == gene)


# 1. the planted FKTN SVA
r = run([rec("chr20", 5040478, "SINE/Alu"), FKTN_SVA])
s = status(r, "FKTN")
check("1 the planted FKTN SVA is seen", s["status"], "seen")
check("1 its call names position, subfamily, copies and filter", s["call"], "chr9:105639660 SVA_F 1 copy PASS")
check("1 its depth is the bedcov sum over the window length", s["mean_depth"], "30.0")
check("1 without it, FKTN is not seen", status(run([rec("chr20", 5040478, "SINE/Alu")]), "FKTN")["status"],
      "not seen")
check("1 the same SVA 1 kb past the window is not seen",
      status(run([rec("chr9", 105641200, "Retroposon/SVA")]), "FKTN")["status"], "not seen")
check("1 an Alu at the FKTN site is not the catalog's SVA",
      status(run([rec("chr9", 105639660, "SINE/Alu", mei="MEI=AluYa5,0/288,+/+")]), "FKTN")["status"], "not seen")
check("1 two copies read as two", status(run([rec("chr9", 105639660, "Retroposon/SVA", cn="2")]), "FKTN")["call"],
      "chr9:105639660 SVA_F 2 copies PASS")
check("1 an SVA MEGAnE could not identify (MEPRED=FAILED) is not seen",
      status(run([rec("chr9", 105639660, "Retroposon/SVA", mepred="FAILED")]), "FKTN")["status"], "not seen")

# 2. too few reads
thin = "chr9\t105639155\t105640156\tFKTN\t9009\t60\n"
check("2 under 10x over the window: not callable", status(run([], depth=thin), "FKTN")["status"], "not callable")
check("2 and its depth is shown", status(run([], depth=thin), "FKTN")["mean_depth"], "9.0")
check("2 a contig the BAM does not have: not callable", status(run([], depth=thin), "TAF1")["status"],
      "not callable")
check("2 with no depth to show", status(run([], depth=thin), "TAF1")["mean_depth"], ".")
check("2 a call at low depth is still seen", status(run([FKTN_SVA], depth=thin), "FKTN")["status"], "seen")
check("2 an empty depth file (no catalog contig in the BAM): every row not callable",
      [x["status"] for x in run([FKTN_SVA], depth="")[1]], ["not callable", "not callable"])

# 3. the calls table
r = run([rec("chr20", 5040478, "SINE/Alu", mei="MEI=AluYb8,0/298,+/+|AluYb9,0/298,+/+", filt="LC", conf="L"),
         rec("chr2", 1000, "LINE/L1", mei="MEI_lbp=L1HS,5581,+;MEI_rbp=pA", cn="2"),
         rec("chr2", 2000, "LINE/L1", mei="MEI_lbp=NA,NA,+;MEI_rbp=pA"),
         rec("chrX", 3000, "Retroposon/SVA", mei="MEI=SVA_E,1/900,+/+"),
         rec("chr20", 6665193, "DNA/hAT-Charlie", mei="MEI=Charlie5,2546/2611,+/+"),
         rec("chr20", 16134625, "LTR/ERV1", mei="MEI=MER41A,91/60,-/-"),
         rec("chr20", 28576054, "LINE/L1", mepred="FAILED"),
         rec("chr20", 30000000, "SINE/Alu", cn="0")])
check("3 Alu, L1 and SVA with MEPRED=PASS and copies, in karyotype order",
      [(c["chrom"], c["pos"], c["element"]) for c in r[0]],
      [("chr2", "1000", "L1"), ("chr2", "2000", "L1"), ("chr20", "5040478", "Alu"), ("chrX", "3000", "SVA")])
check("3 the subfamily: the first MEI name", r[0][2]["subfamily"], "AluYb8")
check("3 the subfamily from a breakpoint when there is no MEI", r[0][0]["subfamily"], "L1HS")
check("3 no name anywhere: '.'", r[0][1]["subfamily"], ".")
check("3 copies, filter, confidence and end are MEGAnE's",
      {k: r[0][2][k] for k in ("copies", "filter", "confidence", "end")},
      {"copies": "1", "filter": "LC", "confidence": "L", "end": "5040486"})
check("3 the columns", list(r[0][0].keys()),
      ["chrom", "pos", "end", "element", "subfamily", "copies", "filter", "confidence"])
check("3 a VCF with no record: an empty table and every site not seen",
      (run([])[0], [x["status"] for x in run([])[1]]), ([], ["not seen", "not seen"]))

# 4. damaged input
err = run([], vcf=FKTN_SVA)
check("4 a VCF without #CHROM stops the script", isinstance(err, str) and "not a VCF" in err, True)
err = run([], vcf=HEAD + "chr9\t1\t.\tN\t<INS:ME>\t.\tPASS\tEND=5\tCN\t1\n")
check("4 a VCF that is not MEGAnE's stops it", isinstance(err, str) and "not MEGAnE's VCF" in err, True)
err = run([], depth="chr9\t105639155\t105640156\n")
check("4 a depth file that is not bedcov -c output stops it", isinstance(err, str) and "samtools bedcov" in err, True)
err = run([], sites="gene\tchrom\nFKTN\tchr9\n")
check("4 a catalog without its columns stops it", isinstance(err, str) and "column" in err, True)

# 5. the summary and the reports
import collect_summary as cs  # noqa: E402
import mei_report as mr  # noqa: E402
import render_report as rr  # noqa: E402
with tempfile.TemporaryDirectory() as d:
    os.makedirs(f"{d}/S")
    summ = cs.collect("S", f"{d}/S")
    check("5 without step 50 the not-assessed list names mobile element insertions",
          cs.MEI_NOT_ASSESSED in summ["not_assessed"] and cs.MEI_NOT_ASSESSED.startswith("Mobile element insertions"),
          True)
    check("5 no table: no section", cs.sec_mei(f"{d}/S", "S"), (None, {}))
    # A step 50 that ran and failed was looked for: it is listed as failed, not as opt-in.
    os.makedirs(f"{d}/F/logs")
    with open(f"{d}/F/logs/run_status.tsv", "w") as f:
        f.write(f"meta\tstarted_epoch\t{time.time() - 60}\nstep\t50\tfailed (exit 1)\n")
    summ = cs.collect("F", f"{d}/F")
    check("5 a failed step 50 is listed under not run, not as not looked for",
          (summ["sections"]["mei"]["state"], cs.MEI_NOT_ASSESSED in summ["not_assessed"],
           any("step 50 failed" in x for x in summ["not_run"])), ("failed", False, True))
    os.makedirs(f"{d}/S/mei")
    for name, text in (("m.vcf", HEAD + FKTN_SVA + rec("chr20", 5040478, "SINE/Alu")), ("sites.tsv", SITES),
                       ("depth.tsv", DEPTH_30X)):
        with open(f"{d}/{name}", "w") as f:
            f.write(text)
    mr.main(["--sample", "S", "--vcf", f"{d}/m.vcf", "--sites", f"{d}/sites.tsv", "--depth", f"{d}/depth.tsv",
             "--out-calls", f"{d}/S/mei/S_mei.tsv", "--out-sites", f"{d}/S/mei/S_mei_sites.tsv"])
    path, data = cs.sec_mei(f"{d}/S", "S")
    check("5 the summary finds the table", path.endswith("mei/S_mei_sites.tsv"), True)
    check("5 the summary counts the calls per element", (data["calls"], data["by_element"]), (2, {"Alu": 1, "SVA": 1}))
    check("5 the summary has each catalog row", [(x["gene"], x["status"]) for x in data["sites"]],
          [("FKTN", "seen"), ("TAF1", "not seen")])
    summ = cs.collect("S", f"{d}/S")
    check("5 with step 50's tables, mobile element insertions are not listed as not assessed",
          cs.MEI_NOT_ASSESSED in summ["not_assessed"], False)
    lines = rr.mei_lines(data)
    check("5 the text block opens with the counts", lines[0], "Insertions called: 2 (Alu 1, SVA 1)")
    check("5 the FKTN line", lines[1], "FKTN SVA (Fukuyama congenital muscular dystrophy; autosomal recessive): "
          "seen (depth 30.0, chr9:105639660 SVA_F 1 copy PASS)")
    check("5 and it ends with the note", lines[-1], rr.MEI_NOTE)
    html = "".join(rr.mei_body(data))
    check("5 the HTML card has the FKTN row and the call", "<td>FKTN</td>" in html and "SVA_F" in html, True)
    with open(f"{d}/S/mei/S_mei_sites.tsv", "w") as f:
        f.write("gene\tstatus\nFKTN\tseen\n")
    try:
        cs.sec_mei(f"{d}/S", "S")
        check("5 a sites table without its columns is unreadable", "read", "Unreadable")
    except cs.Unreadable:
        check("5 a sites table without its columns is unreadable", "Unreadable", "Unreadable")

# 6. the real catalog
rows = mr.read_sites(CATALOG)
check("6 the catalog: FKTN, TAF1 and BBS1", [x["gene"] for x in rows], ["FKTN", "TAF1", "BBS1"])
check("6 every row is an SVA", {x["element"] for x in rows}, {"SVA"})
check("6 every window is the site -/+ 500",
      all(int(x["start"]) == int(x["site"]) - 500 and int(x["end"]) == int(x["site"]) + 500 for x in rows), True)
check("6 the script reads it", [x["status"] for x in run([FKTN_SVA], sites=open(CATALOG).read())[1]],
      ["seen", "not seen", "not callable"])

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
