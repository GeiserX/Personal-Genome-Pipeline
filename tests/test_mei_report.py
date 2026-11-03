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
     HTML card read the tables, and the full HTML report renders the card;
     an empty calls table, or one without its columns, is unreadable;
     without the tables the summary's "not assessed" list names mobile
     element insertions, with them it does not
  6. the real catalog: every row is an SVA, its window is the site -/+ 500,
     and the script reads it
  7. AnnotSV: --annotsv-input writes the listed calls, and only them, as
     <INS:ME> records with END=POS, an ID per call and the VCF's contigs;
     --annotsv joins a table in AnnotSV's columns back by that ID: genes,
     the location in each gene, the benign insertion frequency and the
     class; a call AnnotSV did not annotate, and every call without
     --annotsv, has "." there; a table that is not AnnotSV's stops it; an
     empty one (AnnotSV had nothing to annotate) leaves every call "."

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
      ["chrom", "pos", "end", "element", "subfamily", "copies", "filter", "confidence",
       "genes", "location", "benign_af", "acmg_class"])
check("3 without --annotsv the AnnotSV columns are '.'",
      {r[0][0][k] for k in ("genes", "location", "benign_af", "acmg_class")}, {"."})
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
    try:
        page = rr.html_report(summ)
    except Exception as e:  # the whole step 24 page is lost if the card raises
        page = f"html_report raised {e!r}"
    check("5 the full HTML report renders the card", "Mobile Element Insertions" in page and "SVA_F" in page, True)
    lines = rr.mei_lines(data)
    check("5 the text block opens with the counts", lines[0], "Insertions called: 2 (Alu 1, SVA 1)")
    check("5 the FKTN line", lines[1], "FKTN SVA (Fukuyama congenital muscular dystrophy; autosomal recessive): "
          "seen (depth 30.0, chr9:105639660 SVA_F 1 copy PASS)")
    check("5 and it ends with the note", lines[-1], rr.MEI_NOTE)
    html = "".join(rr.mei_body(data))
    check("5 the HTML card has the FKTN row and the call", "<td>FKTN</td>" in html and "SVA_F" in html, True)
    calls_path = f"{d}/S/mei/S_mei.tsv"
    good_calls = open(calls_path).read()
    with open(calls_path, "w") as f:
        f.write(good_calls.splitlines()[0] + "\n")
    check("5 a calls table with its header and no row is zero calls", cs.sec_mei(f"{d}/S", "S")[1]["calls"], 0)
    for name, text in (("an empty calls table", ""), ("a calls table without its columns", "chrom\tpos\n")):
        with open(calls_path, "w") as f:
            f.write(text)
        try:
            cs.sec_mei(f"{d}/S", "S")
            check(f"5 {name} is unreadable", "read", "Unreadable")
        except cs.Unreadable:
            check(f"5 {name} is unreadable", "Unreadable", "Unreadable")
    with open(calls_path, "w") as f:
        f.write(good_calls)
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

# 7. AnnotSV
def annotsv_run(records, table):
    """(AnnotSV input records, calls) for the records and an AnnotSV table text, or the error text."""
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/m.vcf", "w") as f:
            f.write(HEAD.replace("##ALT", "##contig=<ID=chr20,length=64444167>\n##ALT", 1) + "".join(records))
        p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--vcf", f"{d}/m.vcf",
                            "--annotsv-input", f"{d}/in.vcf"], capture_output=True, text=True)
        if p.returncode:
            return p.stderr
        vcf = open(f"{d}/in.vcf").read()
        for name, text in (("sites.tsv", SITES), ("depth.tsv", DEPTH_30X), ("a.tsv", table)):
            with open(f"{d}/{name}", "w") as f:
                f.write(text)
        p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--vcf", f"{d}/m.vcf", "--sites",
                            f"{d}/sites.tsv", "--depth", f"{d}/depth.tsv", "--annotsv", f"{d}/a.tsv",
                            "--out-calls", f"{d}/c.tsv", "--out-sites", f"{d}/s.tsv"], capture_output=True, text=True)
        if p.returncode:
            return p.stderr
        return vcf, list(csv.DictReader(open(f"{d}/c.tsv"), delimiter="\t"))


AHEAD = "AnnotSV_ID\tSV_chrom\tSV_start\tSV_end\tSV_type\tID\tAnnotation_mode\tGene_name\tLocation\tB_ins_AFmax\tACMG_class\n"
TABLE = (AHEAD
         + "20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tfull\tFIXG1;FIXG2\t\t0.0312\tNA\n"
         + "20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tsplit\tFIXG1\tintron2-intron2\t\tfull=NA\n"
         + "20_5040477_5040478_INS_1\t20\t5040477\t5040478\tINS\tchr20_5040478_Alu\tsplit\tFIXG2\texon1-exon1\t\tfull=NA\n"
         + "9_105639659_105639660_INS_1\t9\t105639659\t105639660\tINS\tchr9_105639660_SVA\tfull\t\t\t\tNA\n")
RECS = [rec("chr20", 5040478, "SINE/Alu", mei="MEI=AluYa5,0/288,+/+"), FKTN_SVA,
        rec("chr20", 7000000, "LINE/L1", cn="2"), rec("chr20", 8000000, "SINE/Alu", mepred="FAILED"),
        rec("chr20", 9000000, "DNA/hAT-Charlie")]
r = annotsv_run(RECS, TABLE)
body = [x.split("\t") for x in r[0].splitlines() if not x.startswith("#")]
check("7 the AnnotSV input holds the listed calls only (not MEPRED=FAILED, not a DNA transposon)",
      [x[2] for x in body], ["chr20_5040478_Alu", "chr9_105639660_SVA", "chr20_7000000_L1"])
check("7 each is an <INS:ME> with END=POS and its FILTER, one copy 0/1, two 1/1",
      [(x[4], x[6], x[7], x[9]) for x in body],
      [("<INS:ME>", "PASS", "SVTYPE=INS;END=5040478", "0/1"), ("<INS:ME>", "PASS", "SVTYPE=INS;END=105639660", "0/1"),
       ("<INS:ME>", "PASS", "SVTYPE=INS;END=7000000", "1/1")])
check("7 the input keeps the VCF's contig lines", "##contig=<ID=chr20,length=64444167>" in r[0], True)
calls = {c["pos"]: c for c in r[1]}
check("7 the Alu: its genes, the location in each, the benign frequency and the class",
      {k: calls["5040478"][k] for k in ("genes", "location", "benign_af", "acmg_class")},
      {"genes": "FIXG1,FIXG2", "location": "FIXG1:intron2-intron2,FIXG2:exon1-exon1", "benign_af": "0.0312",
       "acmg_class": "NA"})
check("7 the SVA: annotated, in no gene: class NA, the rest '.'",
      {k: calls["105639660"][k] for k in ("genes", "location", "benign_af", "acmg_class")},
      {"genes": ".", "location": ".", "benign_af": ".", "acmg_class": "NA"})
check("7 the L1 AnnotSV did not annotate: '.'",
      {calls["7000000"][k] for k in ("genes", "location", "benign_af", "acmg_class")}, {"."})
check("7 red control: with an ID that matches no call, the Alu is '.'",
      annotsv_run(RECS, TABLE.replace("chr20_5040478_Alu", "chr20_5040479_Alu"))[1][1]["genes"], ".")
check("7 an empty AnnotSV table: every call '.'",
      {c["acmg_class"] for c in annotsv_run(RECS, "")[1]}, {"."})
quoted = annotsv_run(RECS, TABLE.replace("\tsplit\tFIXG1\t", "\tsplit\t\"FIXG1\t"))[1]
check("7 a field that starts with a double quote is read as it is, not as a quoted field",
      ({c["pos"]: c for c in quoted}["105639660"]["acmg_class"], {c["pos"]: c for c in quoted}["5040478"]["genes"]),
      ("NA", "FIXG1,FIXG2"))
err = annotsv_run(RECS, "SV_chrom\tSV_start\n20\t1\n")
check("7 a table that is not AnnotSV's stops it", isinstance(err, str) and "not AnnotSV's table" in err, True)
dup = annotsv_run([rec("chr20", 100, "SINE/Alu"), rec("chr20", 100, "SINE/Alu")], AHEAD)
check("7 two calls at one place get two IDs", [x.split("\t")[2] for x in dup[0].splitlines() if x[0] != "#"],
      ["chr20_100_Alu", "chr20_100_Alu_2"])
with tempfile.TemporaryDirectory() as d:
    with open(f"{d}/m.vcf", "w") as f:
        f.write(HEAD)
    p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--vcf", f"{d}/m.vcf"], capture_output=True, text=True)
check("7 without --annotsv-input the report options are required", p.returncode == 2 and "--sites" in p.stderr, True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
