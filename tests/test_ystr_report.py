#!/usr/bin/env python3
"""bin/ystr_report.py: Y-STR repeat counts from HipSTR's calls (step 45).

Each case writes a synthetic HipSTR VCF (one sample, only the records it
plants), runs the script on assets/ystr_panel.tsv, and checks the table:

  1. the panel: 111 markers in FamilyTreeDNA's order (the Y-12, Y-25 and
     Y-37 heads checked by name), every multi-copy marker kind multicopy,
     the BED holds every region of the single and parts markers once, and
     one region over the whole locus for each span marker
  2. counts: reference copies plus GB / motif length; a partial repeat as
     N.r; a negative GB; the sum of a marker's parts (DYS449); a span
     marker's parts' copies plus its one GB (DYS510, DYS389II at 250 bp);
     an offset from a panel that gives one
  3. confidence: high; low for a posterior under 0.9, for over 15 % of
     reads with stutter or with a flank indel, and for under 10 reads; no
     call for a region without a record or without a genotype, for two
     lengths on the haploid Y, and for a marker whose other part is missing;
     a region without a record carries the reason HipSTR's log gives for it
  4. not callable: multi-copy markers (DYS464, DYF395S1), markers the
     marker set has no region for (DYS439), and a region longer than the
     reads can span (DYS452 at 150 bp reads, callable at 250 bp; the spans
     of DYS389II and DYS448)
  5. the CSV: one row per marker in the panel's order, the count and the
     confidence
  6. a damaged input: a VCF with two samples, a panel out of order
  7. the reports: bin/collect_summary.py reads the table (counts per
     confidence) and refuses one with an unknown confidence; the HTML card of
     bin/render_report.py lists the called markers only, with the note, and
     the HTML report shows an unreadable table as Unreadable

Run: python3 tests/test_ystr_report.py
"""
import atexit
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "ystr_report.py")
PANEL = os.path.join(REPO, "assets", "ystr_panel.tsv")
sys.path.insert(0, os.path.join(REPO, "bin"))
import ystr_report as yr  # noqa: E402

FAILS = []


def check(name, got, want):
    if got == want:
        print(f"[PASS] {name}")
    else:
        print(f"[FAIL] {name}: got {got!r}, want {want!r}")
        FAILS.append(name)


TMP = tempfile.mkdtemp(prefix="ystr-test-")
atexit.register(shutil.rmtree, TMP, ignore_errors=True)
FORMAT = "GT:GB:Q:PQ:DP:DSNP:PSNP:PDP:GLDIFF:DSTUTTER:DFLANKINDEL"


def vcf(records, samples=("S1",), gz=True):
    """records: {HipSTR name: (GT, GB, Q, DP, DSTUTTER, DFLANKINDEL)}."""
    panel = yr.read_panel(PANEL)
    where = {n: (c, s) for r in panel for c, s, _e, _p, _k, n in r["regions"]}
    lines = ["##fileformat=VCFv4.1", "##source=HipSTR",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)]
    for name, (gt, gb, q, dp, st, fl) in records.items():
        c, s = where.get(name, ("chrY", 1))
        cell = f"{gt}:{gb}:{q}:{q}:{dp}:0:0|0:{dp}|0:.:{st}:{fl}"
        lines.append("\t".join([c, str(s - 1), name, "A", ".", ".", "PASS", "PERIOD=4", FORMAT]
                               + [cell] * len(samples)))
    p = os.path.join(TMP, f"calls{len(os.listdir(TMP))}.vcf" + (".gz" if gz else ""))
    data = ("\n".join(lines) + "\n").encode()
    with (gzip.open(p, "wb") if gz else open(p, "wb")) as f:
        f.write(data)
    return p


def run(records=None, panel=PANEL, read_length=None, samples=("S1",), want_rc=0, log=None):
    v = vcf(records or {}, samples)
    tsv, out_csv = v + ".tsv", v + ".csv"
    cmd = [sys.executable, SCRIPT, "report", "--panel", panel, "--vcf", v, "--tsv", tsv, "--csv", out_csv]
    if read_length:
        cmd += ["--read-length", str(read_length)]
    if log is not None:
        lp = v + ".log"
        with open(lp, "w") as f:
            f.write(log)
        cmd += ["--log", lp]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != want_rc:
        print(p.stdout, p.stderr)
    check(f"exit code {want_rc} ({' '.join(cmd[2:3])})", p.returncode, want_rc)
    if want_rc:
        return p
    with open(tsv, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(out_csv, encoding="utf-8") as f:
        crow = list(csv.reader(f))
    return {r["marker"]: r for r in rows}, rows, crow, p.stdout


GOOD = ("0", "0", "0.99", "20", "0", "0")


def good(gb="0", dp="20"):
    return ("0" if gb == "0" else "1", gb, "0.99", dp, "0", "0")


# --- 1. the panel --------------------------------------------------------------------
panel = yr.read_panel(PANEL)
names = [r["marker"] for r in panel]
check("1 111 markers", len(panel), 111)
check("1 Y-12 in FamilyTreeDNA's order", names[:12],
      ["DYS393", "DYS390", "DYS19", "DYS391", "DYS385a", "DYS385b", "DYS426", "DYS388", "DYS439",
       "DYS389I", "DYS392", "DYS389II"])
check("1 Y-25 (13-25)", names[12:25],
      ["DYS458", "DYS459a", "DYS459b", "DYS455", "DYS454", "DYS447", "DYS437", "DYS448", "DYS449",
       "DYS464a", "DYS464b", "DYS464c", "DYS464d"])
check("1 Y-37 (26-37)", names[25:37],
      ["DYS460", "Y-GATA-H4", "YCAIIa", "YCAIIb", "DYS456", "DYS607", "DYS576", "DYS570", "CDYa", "CDYb",
       "DYS442", "DYS438"])
check("1 the last marker is DYS435", names[-1], "DYS435")
check("1 panel sizes", [sum(1 for r in panel if r["panel"] == p) for p in ("Y-12", "Y-25", "Y-37", "Y-67", "Y-111")],
      [12, 13, 12, 30, 44])
multi = {r["marker"] for r in panel if r["kind"] == "multicopy"}
check("1 the multi-copy markers", multi,
      {"DYS385a", "DYS385b", "DYS459a", "DYS459b", "DYS464a", "DYS464b", "DYS464c", "DYS464d", "YCAIIa", "YCAIIb",
       "CDYa", "CDYb", "DYF395S1a", "DYF395S1b", "DYS413a", "DYS413b"})
bed = os.path.join(TMP, "regions.bed")
check("1 bed exits 0", yr.main(["bed", "--panel", PANEL, "--out", bed]), 0)
with open(bed) as f:
    bl = [line.rstrip("\n").split("\t") for line in f]
want_regions = {n for r in panel for *_x, n in r["regions"]}
check("1 the BED has every region once", sorted(x[5] for x in bl), sorted(want_regions))
check("1 the BED has six columns, 1-based start below the end", all(len(x) == 6 and 0 < int(x[1]) < int(x[2]) for x in bl),
      True)
check("1 the BED is sorted", [int(x[1]) for x in bl], sorted(int(x[1]) for x in bl))
check("1 DYS393's region", [x for x in bl if x[5] == "DYS393"], [["chrY", "3263111", "3263158", "4", "12", "DYS393"]])
spans = {r["marker"] for r in panel if r["kind"] == "span"}
check("1 the span markers", spans, {"DYS389II", "DYS448", "DYS450", "DYS504", "DYS510", "DYS532"})
check("1 the BED: 86 regions (76 single, 6 spans, the 4 parts of DYS449 and DYS552)", len(bl), 86)
check("1 DYS450 is one region from its first part's start to its last part's end, 9 + 5 copies",
      [x for x in bl if x[5].startswith("DYS450")], [["chrY", "8258259", "8258341", "5", "14", "DYS450"]])
check("1 DYS510's three parts are one region, 17 + 5 + 3 copies (4.5 halves up)",
      [x for x in bl if x[5].startswith("DYS510")], [["chrY", "15188021", "15188138", "4", "25", "DYS510"]])
check("1 no part of a span marker is in the BED", [x[5] for x in bl if x[5].split(".")[0] in spans and "." in x[5]], [])
check("1 DYS449 keeps its two parts", sorted(x[5] for x in bl if x[5].startswith("DYS449")), ["DYS449.1", "DYS449.2"])

# --- 2. counts -------------------------------------------------------------------------
by, rows, crow, out = run({
    "DYS393": good("4"),                     # one repeat longer than the reference's 12
    "DYS390": good("-8"),                    # two shorter than 24
    "DYS391": good("2"),                     # a partial repeat: 11.2
    "DYS392": good("-1"),                    # three-base motif, one base short: 12.2
    "DYS449.1": good("4"), "DYS449.2": good("-4"),      # 15 (15.25 rounded) + 1, plus 21 - 1
    "DYS510": good("4"),                                # one region: 17 + 5 (4.5 halves up) + 3, plus 1
    "DYS389II": good("4"), "DYS448": good("6"),         # spans longer than 150 bp reads can cover
    "DYS19/DYS394": GOOD,
})
check("2 a reference call keeps the reference copies (DYS19: 16)", by["DYS19"]["count"], "16")
check("2 GB +4 at a 4-base motif: one repeat more (DYS393 13)", by["DYS393"]["count"], "13")
check("2 GB -8: two repeats fewer (DYS390 22)", by["DYS390"]["count"], "22")
check("2 GB +2: a partial repeat (DYS391 11.2)", by["DYS391"]["count"], "11.2")
check("2 GB -1 at a 3-base motif (DYS392 12.2)", by["DYS392"]["count"], "12.2")
check("2 DYS449 is the sum of its two parts (16 + 20)", by["DYS449"]["count"], "36")
check("2 DYS510: its parts' copies (17 + 5 + 3, 4.5 halves up) plus its one GB", by["DYS510"]["count"], "26")
check("2 DYS510 reads one region: one read count, no part names in the note",
      (by["DYS510"]["confidence"], by["DYS510"]["reads"], by["DYS510"]["note"]), ("high", "20", "."))
for m, bp in (("DYS389II", 164), ("DYS448", 156)):
    check(f"2 {m}: the span is {bp} bp, not callable from 150 bp reads",
          (by[m]["confidence"], by[m]["count"], f"a {bp} bp repeat" in by[m]["note"]), ("not callable", "", True))
by250 = run({"DYS389II": good("4")}, read_length=250)[0]
check("2 DYS389II at 250 bp reads: 17 + 12 copies, plus 1", by250["DYS389II"]["count"], "30")
check("2 the DYS389II note says the count covers DYS389I", "covers DYS389I" in by250["DYS389II"]["note"], True)
check("2 the summary line counts the confidences", out.strip().startswith("Y-STR markers: "), True)
off = os.path.join(TMP, "offset_panel.tsv")
with open(PANEL) as f, open(off, "w") as g:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) == 8 and c[2] == "DYS393":
            c[6] = "-1"
        g.write("\t".join(c) + "\n")
oby = run({"DYS393": good("4")}, panel=off)[0]
check("2 a panel offset is added (13 - 1)", oby["DYS393"]["count"], "12")

# --- 3. confidence ---------------------------------------------------------------------
by = run({
    "DYS393": GOOD,
    "DYS390": ("0", "0", "0.6", "20", "0", "0"),     # posterior below 0.9
    "DYS391": ("0", "0", "0.99", "20", "4", "0"),    # 20 % stutter
    "DYS426": ("0", "0", "0.99", "20", "0", "4"),    # 20 % flank indels
    "DYS388": ("0", "0", "0.99", "6", "0", "0"),     # six reads
    "DYS392": (".", ".", ".", "0", "0", "0"),        # no genotype
    "DYS458": ("0/1", "0/4", "0.99", "20", "0", "0"),  # two lengths
    "DYS449.1": GOOD,                                # DYS449.2 missing
    "DYS456": ("0", "0", "0.99", "20", "3", "3"),    # exactly 15 %: still high
})[0]
check("3 a clean call is high", (by["DYS393"]["confidence"], by["DYS393"]["reads"], by["DYS393"]["q"]),
      ("high", "20", "0.99"))
check("3 posterior 0.6: low", by["DYS390"]["confidence"], "low")
check("3 and the note says why", "posterior 0.6 below 0.9" in by["DYS390"]["note"], True)
check("3 a low call keeps its count", by["DYS390"]["count"], "24")
check("3 4 of 20 reads with stutter: low", (by["DYS391"]["confidence"], by["DYS391"]["stutter_reads"]), ("low", "4"))
check("3 4 of 20 reads with a flank indel: low", by["DYS426"]["confidence"], "low")
check("3 6 reads: low", (by["DYS388"]["confidence"], "6 reads, fewer than 10" in by["DYS388"]["note"]), ("low", True))
check("3 3 of 20 (15 %) is not above the limit: high", by["DYS456"]["confidence"], "high")
check("3 no genotype: no call, no count", (by["DYS392"]["confidence"], by["DYS392"]["count"]), ("no call", ""))
check("3 two lengths on the Y: no call", (by["DYS458"]["confidence"], "two lengths" in by["DYS458"]["note"]),
      ("no call", True))
check("3 a missing part: no call", (by["DYS449"]["confidence"], "DYS449.2" in by["DYS449"]["note"]),
      ("no call", True))
check("3 a region without a record: no call", by["DYS19"]["confidence"], "no call")
check("3 without HipSTR's log, the note points to it", "see HipSTR's log" in by["DYS19"]["note"], True)
# HipSTR's log: a "Processing region" line per region (0-based start), then why it stopped.
LOG = """Processing region chrY 9684379 9684443
Aligning reads to each candidate haplotype
Aborting genotyping of the locus as the sequence downstream of the repeat is too repetitive for accurate genotyping
\tFlanking sequence = AAACACTATATATATATAACACTATATATATAATAC
Processing region chrY 18041621 18041651
Skipping locus with too few reads: TOTAL=3, MIN=5
Processing region chrY 3263110 3263158
Locus timing:
Processing region chrY 8258258 8258341
Skipping locus with too few informative reads for stutter training: TOTAL=2, MIN=5
Processing region chrY 12500447 12500611
Skipping region as the reference allele length exceeds the threshold (164 vs 150)
"""
by = run({"DYS390": GOOD}, log=LOG, read_length=250)[0]
check("3 the log's reason for DYS19 (its flank)", by["DYS19"]["note"],
      "HipSTR wrote no record for DYS19/DYS394; its log: Aborting genotyping of the locus as the sequence "
      "downstream of the repeat is too repetitive for accurate genotyping")
check("3 the log's reason for DYS425 (too few reads)", by["DYS425"]["note"],
      "HipSTR wrote no record for DYS425; its log: Skipping locus with too few reads: TOTAL=3, MIN=5")
check("3 a span's reason, found by its one region (DYS450)",
      "too few informative reads for stutter training" in by["DYS450"]["note"], True)
check("3 a span longer than --max-str-len (DYS389II at 250 bp reads)",
      "exceeds the threshold (164 vs 150)" in by["DYS389II"]["note"], True)
check("3 a region the log names without a reason: the note says so", by["DYS393"]["note"],
      "HipSTR wrote no record for DYS393, and its log gives no reason")
check("3 a region the log does not name: the note says so", by["DYS391"]["note"],
      "HipSTR wrote no record for DYS391, and its log gives no reason")
check("3 a called marker carries no log line", by["DYS390"]["note"], ".")

# --- 4. not callable ---------------------------------------------------------------------
by = run({"DYS452": GOOD})[0]
for m in ("DYS464a", "DYS464d", "DYF395S1a", "DYS385b", "CDYa"):
    check(f"4 {m}: not callable, multi-copy", (by[m]["confidence"], "multi-copy" in by[m]["note"], by[m]["count"]),
          ("not callable", True, ""))
check("4 DYS439 (no region in the marker set): not callable", (by["DYS439"]["confidence"], "no region" in by["DYS439"]["note"]),
      ("not callable", True))
check("4 DYS452 (155 bp) at 150 bp reads: not callable even with a call",
      (by["DYS452"]["confidence"], by["DYS452"]["count"]), ("not callable", ""))
by = run({"DYS452": GOOD}, read_length=250)[0]
check("4 DYS452 at 250 bp reads: called", (by["DYS452"]["confidence"], by["DYS452"]["count"]), ("high", "31"))

# --- 5. the CSV ----------------------------------------------------------------------------
_by, rows, crow, _o = run({"DYS393": good("4")})
check("5 CSV header", crow[0], ["marker", "value", "confidence"])
check("5 CSV: a row per marker in the panel's order", [r[0] for r in crow[1:]], names)
check("5 CSV: DYS393's row", crow[1], ["DYS393", "13", "high"])
check("5 CSV: a multi-copy row has no value", crow[5], ["DYS385a", "", "not callable"])
check("5 TSV: rows in the panel's order", [r["marker"] for r in rows], names)

# --- 6. damaged input --------------------------------------------------------------------
p = run({"DYS393": GOOD}, samples=("S1", "S2"), want_rc=1)
check("6 two samples: the error says so", "2 samples" in p.stderr, True)
bad = os.path.join(TMP, "bad_panel.tsv")
with open(PANEL) as f, open(bad, "w") as g:
    lines = f.readlines()
    head = [x for x in lines if x.startswith("#") or x.startswith("order\t")]
    body = [x for x in lines if x not in head]
    g.writelines(head + [body[1], body[0]] + body[2:])
p = run({}, panel=bad, want_rc=1)
check("6 a panel out of order stops", "order column" in p.stderr, True)

# --- 7. the reports -----------------------------------------------------------------------
import collect_summary as cs  # noqa: E402
import render_report as rr  # noqa: E402

sd = os.path.join(TMP, "sample")
os.makedirs(os.path.join(sd, "y_str"))
v = vcf({"DYS393": good("4"), "DYS390": ("0", "0", "0.6", "20", "0", "0")})
yr.main(["report", "--panel", PANEL, "--vcf", v, "--tsv", os.path.join(sd, "y_str", "S1_y_str.tsv"),
         "--csv", os.path.join(sd, "y_str", "S1_y_str.csv")])
path, data = cs.sec_y_str(sd, "S1")
check("7 collect_summary finds the table", os.path.basename(path or ""), "S1_y_str.tsv")
check("7 counts per confidence", (data["high"], data["low"], data["not_callable"], data["no_call"]),
      (1, 1, sum(1 for r in panel if r["kind"] in ("multicopy", "absent")) + 3,  # + DYS452, DYS389II, DYS448: too long
       111 - 2 - sum(1 for r in panel if r["kind"] in ("multicopy", "absent")) - 3))
check("7 the markers keep their order", [x["marker"] for x in data["markers"]], names)
check("7 a sample without the step has no section", cs.sec_y_str(TMP, "nobody"), (None, {}))
html = "\n".join(rr.ystr_body(data))
check("7 the card lists DYS393 at 13", "<td>DYS393</td><td>13</td><td>high</td>" in html, True)
check("7 the card lists the low call too", "<td>DYS390</td><td>24</td><td>low</td>" in html, True)
check("7 and leaves out markers without a count", "DYS464a" in html, False)
check("7 the card carries the convention note", "fixed number of repeats" in html, True)
with open(os.path.join(sd, "y_str", "S1_y_str.tsv")) as f:
    t = f.read().replace("\thigh\t", "\tgreat\t", 1)
with open(os.path.join(sd, "y_str", "S1_y_str.tsv"), "w") as f:
    f.write(t)
try:
    cs.sec_y_str(sd, "S1")
    check("7 an unknown confidence is unreadable", "read", "Unreadable")
except cs.Unreadable:
    check("7 an unknown confidence is unreadable", "Unreadable", "Unreadable")
try:
    page = rr.html_report(cs.collect("S1", sd))
except (KeyError, TypeError) as e:
    page = f"{type(e).__name__}: {e}"
check("7 the HTML report renders an unreadable Y-STR table as such",
      "<h2>Y-STR Markers (HipSTR)</h2>" in page and "Unreadable" in page, True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
