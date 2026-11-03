#!/usr/bin/env python3
"""bin/bloodtype_report.py: the blood group antigens a genome predicts (step 41).

Each case writes a synthetic one-sample VCF (only the records it plants) and,
when the case has a BAM, a samtools bedcov file over assets/bloodtype_markers.tsv
with chosen depths, runs the script, and checks rows of the table it writes:

  1. the table and the code agree: every site has a rule, every rule a site,
     every phenotype template names both antigens of its site
  2. nothing planted, reads everywhere: the reference alleles (O, k, Fy(a+b-) ...)
     at confidence high; without the BAM the same calls at medium, and RhD
     not determined
  3. ABO: A1, A2, B, AB, O.02, and markers that do not add up
  4. RhD from depth: two copies, one, none (the RHD deletion), the RHD
     pseudogene, a hybrid gene (exons 3 to 7 without reads), low flank depth
  5. Duffy and Kell flip with planted variants: Fyb, the GATA box, FY*X; K
  6. the VCF: a site with too few reads, a FILTER-failed call, a RefCall 0/0,
     contigs without chr, a multi-allelic record, a record whose REF is not the
     marker's REF, a second sample
  7. a BAM without some contigs (no depth rows there): those systems and RhD
     are not determined, the others are read as usual

Run: python3 tests/test_bloodtype_report.py
"""
import csv
import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "bloodtype_report.py")
MARKERS = os.path.join(REPO, "assets", "bloodtype_markers.tsv")
sys.path.insert(0, os.path.join(REPO, "bin"))
import bloodtype_report as bt  # noqa: E402

FAILS = []
SITES, REGIONS = bt.read_markers(MARKERS)


def check(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def vcf_line(name, gt, filt="PASS", chrom=None, alt=None, ref=None):
    r = SITES[name]
    return "\t".join([chrom or r["chrom"], str(r["pos"]), ".", ref or r["ref"], alt or r["alt"], "50", filt, ".",
                      "GT", gt]) + "\n"


def run(records=(), depth=None, gz=False, samples=("S",), raw=None, drop=()):
    """{system: row} of the table the script writes, or the error text."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf = os.path.join(tmp, "s.vcf.gz" if gz else "s.vcf")
        text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(samples) + "\n" + "".join(records) + (raw or ""))
        if gz:
            with gzip.open(vcf, "wt") as f:
                f.write(text)
        else:
            with open(vcf, "w") as f:
                f.write(text)
        cmd = [sys.executable, SCRIPT, "--sample", "S", "--vcf", vcf, "--markers", MARKERS,
               "--out", os.path.join(tmp, "t.tsv"), "--evidence", os.path.join(tmp, "e.tsv")]
        if depth is not None:
            dp = os.path.join(tmp, "d.tsv")
            with open(dp, "w") as f:
                for r in list(SITES.values()) + list(REGIONS.values()):
                    if r["chrom"] in drop:
                        continue   # a contig the BAM has no @SQ line for: step 41 leaves its rows out
                    mean = depth.get(r["name"], depth.get("*", 30.0))
                    cols = [r["chrom"], str(r["start"]), str(r["end"]), r["kind"], r["system"], r["name"]]
                    f.write("\t".join(cols + ["", "", "", "", "", "", "",
                                              str(int(mean * (r["end"] - r["start"])))]) + "\n")
            cmd += ["--depth", dp]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            return p.stderr.strip() or f"exit {p.returncode}"
        with open(os.path.join(tmp, "t.tsv")) as f:
            return {r["system"]: r for r in csv.DictReader(f, delimiter="\t")}


def rhd(copies, flank=30.0, exons=None):
    """Depth means for an RHD copy number (2 = the flanks' depth)."""
    d = {"*": 30.0, "RHD_FLANK_UP": flank, "RHD_FLANK_DOWN": flank, "RHD_GENE": flank * copies / 2}
    for n in REGIONS:
        if n.startswith("RHD_EX"):
            d[n] = flank * (exons or {}).get(int(n[6:]), copies) / 2
    return d


def ph(t, system):
    return t[system]["phenotype"] if isinstance(t, dict) else t


def cf(t, system):
    return t[system]["confidence"] if isinstance(t, dict) else t


# 1. table and code
names = set(SITES)
NULLS = {n for n, r in SITES.items() if r["alt_means"] == "null"}
ruled = set(bt.PAIRS) | bt.SPECIAL_SITES | NULLS
check("1 every site of the table has a rule", sorted(names - ruled), [])
check("1 every rule has a site in the table", sorted(ruled - names), [])
for n, (label, tmpl) in bt.PAIRS.items():
    r = SITES.get(n)
    if r:
        check(f"1 {n}: the template names both antigens",
              ("{" + r["ref_means"] + "}" in tmpl, "{" + r["alt_means"] + "}" in tmpl), (True, True))
check("1 every system has a site", sorted({s for s, _ in bt.SYSTEMS} - {r["system"] for r in SITES.values()}), [])
check("1 positions are unique", len({(r["chrom"], r["pos"]) for r in SITES.values()}), len(SITES))
check("1 the RhD depth regions are there",
      sorted(n for n in REGIONS if not n.startswith("RHD_EX")), ["RHD_FLANK_DOWN", "RHD_FLANK_UP", "RHD_GENE"])
check("1 RHD has ten exons", sorted(int(n[6:]) for n in REGIONS if n.startswith("RHD_EX")), list(range(1, 11)))
check("1 every null row is in a system that reads null rows",
      sorted(n for n in NULLS if SITES[n]["system"] not in set(bt.NULL_SYSTEMS) | bt.NULL_HOSTS), [])
check("1 every system read from null rows has some",
      sorted(s for s in bt.NULL_SYSTEMS if not any(SITES[n]["system"] == s for n in NULLS)), [])
check("1 a null row of Kell or Kidd names the alleles it sits with",
      sorted(n for n in NULLS if SITES[n]["system"] in bt.NULL_HOSTS and not SITES[n]["null_on"]), [])
check("1 the high_prevalence antigen is one the row gives",
      sorted(n for n, r in SITES.items() if r["high_prevalence"] and r["high_prevalence"] not in (r["ref_means"], r["alt_means"])), [])

# 2. nothing planted
ref = run(depth=rhd(2))
want = {"ABO": "O", "RhD": "D+", "RhCE": "C-c+ E-e+", "Kell": "K-k+ Kp(a-b+) Js(a-b+)", "Duffy": "Fy(a+b-)",
        "Kidd": "Jk(a+b-)", "MNS": "M-N+ S-s+", "Lutheran": "Lu(a-b+)", "Diego": "Di(a-b+)", "Dombrock": "Do(a-b+)",
        "Colton": "Co(a+b-)", "HPA-1": "HPA-1(a+b-)", "HPA-2": "HPA-2(a+b-)", "HPA-3": "HPA-3(a+b-)",
        "HPA-4": "HPA-4(a+b-)", "HPA-5": "HPA-5(a+b-)", "HPA-15": "HPA-15(a-b+)"}
check("2 reference alleles, reads everywhere: every system", {k: ph(ref, k) for k in want}, want)
check("2 ABO genotype", ref["ABO"]["genotype"] if isinstance(ref, dict) else ref, "O.01/O.01")
check("2 ABO confidence high with reads", cf(ref, "ABO"), "high")
check("2 RhCE low: c.307 sits in RHCE's RHD-like exon 2", cf(ref, "RhCE"), "low")
check("2 Kell high", cf(ref, "Kell"), "high")
nobam = run()
check("2 without the BAM, ABO is still O", ph(nobam, "ABO"), "O")
check("2 without the BAM, confidence medium", cf(nobam, "ABO"), "medium")
check("2 without the BAM, RhD is not determined", ph(nobam, "RhD"), "not determined")
check("2 without the BAM, the note says why", "deletes the whole RHD gene" in nobam["RhD"]["note"], True)

# 3. ABO
check("3 A1/O: one c.261 insertion", ph(run([vcf_line("ABO_261", "0/1")], rhd(2)), "ABO"), "A1")
a2 = run([vcf_line("ABO_261", "1/1"), vcf_line("ABO_1061", "1/1")], rhd(2))
check("3 A2/A2: c.1061delC on both", (ph(a2, "ABO"), a2["ABO"]["genotype"]), ("A2", "A2/A2"))
check("3 A1/A2 reads A1", ph(run([vcf_line("ABO_261", "1/1"), vcf_line("ABO_1061", "0/1")], rhd(2)), "ABO"), "A1")
b = run([vcf_line("ABO_261", "0/1"), vcf_line("ABO_803", "0/1"), vcf_line("ABO_796", "0/1")], rhd(2))
check("3 B/O", (ph(b, "ABO"), b["ABO"]["genotype"]), ("B", "B/O.01"))
ab = run([vcf_line("ABO_261", "1/1"), vcf_line("ABO_803", "0/1"), vcf_line("ABO_796", "0/1")], rhd(2))
check("3 A1B", (ph(ab, "ABO"), cf(ab, "ABO")), ("A1B", "high"))
check("3 the A1 or A2 caveat is in the note", "c.1061delC, a deletion in a run of C" in ab["ABO"]["note"], True)
o2 = run([vcf_line("ABO_261", "0/1"), vcf_line("ABO_802", "0/1")], rhd(2))
check("3 O.01/O.02 is O", (ph(o2, "ABO"), o2["ABO"]["genotype"]), ("O", "O.01/O.02"))
bad = run([vcf_line("ABO_261", "0/1"), vcf_line("ABO_803", "1/1")], rhd(2))
check("3 two B markers on one functional allele: not determined", ph(bad, "ABO"), "not determined")
odd = run([vcf_line("ABO_261", "0/1"), vcf_line("ABO_803", "0/1")], rhd(2))
check("3 c.803 without c.796: B, confidence low", (ph(odd, "ABO"), cf(odd, "ABO")), ("B", "low"))

# 4. RhD from depth
check("4 two copies: D+", ph(run(depth=rhd(2)), "RhD"), "D+")
het = run(depth=rhd(1))
check("4 one copy (one RHD deleted): D+", (ph(het, "RhD"), het["RhD"]["genotype"].split(" (")[0]),
      ("D+", "RHD/RHD deleted"))
gone = run(depth=rhd(0))
check("4 no copy (RHD deleted on both): D-", (ph(gone, "RhD"), cf(gone, "RhD")), ("D-", "high"))
check("4 the deleted genotype", gone["RhD"]["genotype"].split(" (")[0] if isinstance(gone, dict) else gone,
      "RHD deleted/RHD deleted")
psi = run([vcf_line("RHD_PSI_DUP", "0/1"), vcf_line("RHD_PSI_609", "0/1")], rhd(1))
check("4 one RHD copy and it is the pseudogene: D-", ph(psi, "RhD"), "D-")
psi2 = run([vcf_line("RHD_PSI_DUP", "0/1"), vcf_line("RHD_PSI_609", "0/1")], rhd(2))
check("4 two copies, one of them the pseudogene: D+", ph(psi2, "RhD"), "D+")
hyb = run(depth=rhd(1.6, exons={3: 0, 4: 0, 5: 0, 6: 0, 7: 0}))
check("4 exons 3 to 7 without reads: a hybrid, RhD not determined", ph(hyb, "RhD"), "not determined")
check("4 the hybrid note names the exons", "exons 3, 4, 5, 6, 7" in hyb["RhD"]["note"], True)
check("4 flanks under 10x: not determined", ph(run(depth=rhd(2, flank=6.0)), "RhD"), "not determined")
check("4 near a boundary (1.45 copies): confidence low", cf(run(depth=rhd(1.45)), "RhD"), "low")

# 5. Duffy and Kell flip
check("5 Fyb het: Fy(a+b+)", ph(run([vcf_line("FY_125", "0/1")], rhd(2)), "Duffy"), "Fy(a+b+)")
gata = run([vcf_line("FY_125", "0/1"), vcf_line("FY_GATA", "0/1")], rhd(2))
check("5 Fyb het with the GATA box: Fy(a+b-)", (ph(gata, "Duffy"), gata["Duffy"]["genotype"]), ("Fy(a+b-)", "FY*A/FY*B GATA"))
check("5 Fyb hom, GATA hom: Fy(a-b-)",
      ph(run([vcf_line("FY_125", "1/1"), vcf_line("FY_GATA", "1/1")], rhd(2)), "Duffy"), "Fy(a-b-)")
check("5 FY*X: Fy(a-b+w)",
      ph(run([vcf_line("FY_125", "1/1"), vcf_line("FY_265", "1/1")], rhd(2)), "Duffy"), "Fy(a-b+w)")
k = run([vcf_line("KEL_578", "0/1")], rhd(2))
check("5 K het: K+k+", (ph(k, "Kell"), k["Kell"]["genotype"]), ("K+k+ Kp(a-b+) Js(a-b+)", "K/k Kpb/Kpb Jsb/Jsb"))
check("5 K hom: K+k-", ph(run([vcf_line("KEL_578", "1/1")], rhd(2)), "Kell").split()[0], "K+k-")
check("5 HPA-15: the reference is 15b, the alt 15a",
      ph(run([vcf_line("HPA15_2108", "1/1")], rhd(2)), "HPA-15"), "HPA-15(a+b-)")

# 6. the VCF
low = run(depth=dict(rhd(2), KEL_578=3.0))
check("6 a missing site with 3 reads: that pair is not determined",
      (ph(low, "Kell"), cf(low, "Kell")), ("K/k not determined Kp(a-b+) Js(a-b+)", "low"))
check("6 a FILTER-failed call is no call",
      ph(run([vcf_line("JK_838", "1/1", filt="LowQual")], rhd(2)), "Kidd"), "not determined")
check("6 a RefCall 0/0 is a reference call",
      (ph(run([vcf_line("JK_838", "0/0", filt="RefCall")], rhd(2, flank=30.0)), "Kidd")), "Jk(a+b-)")
check("6 contigs without chr, gzip",
      ph(run([vcf_line("JK_838", "1/1", chrom="18")], rhd(2), gz=True), "Kidd"), "Jk(a-b+)")
multi = run([vcf_line("KEL_578", "1/2", alt="C,A")], rhd(2))
check("6 a multi-allelic record: the second ALT is K", ph(multi, "Kell").split()[0], "K+k+")
check("6 another allele at the site lowers the confidence", cf(multi, "Kell"), "low")
mism = run([vcf_line("KEL_578", "1/1", ref="C", alt="T")], rhd(2))
check("6 a record whose REF is not the marker's REF: that pair is not determined",
      (ph(mism, "Kell").split()[0:3], cf(mism, "Kell")), (["K/k", "not", "determined"], "low"))
check("6 and the note names both REFs", "REF C at this position, the marker's REF is G" in mism["Kell"]["note"], True)
span = run([vcf_line("KEL_578", "0/0", ref="GT", alt="G")], rhd(2))
check("6 a reference call over a longer REF that starts with the marker's REF is reference",
      ph(span, "Kell").split()[0], "K-k+")
check("6 a no-call genotype", ph(run([vcf_line("JK_838", "./.")], rhd(2)), "Kidd"), "not determined")
two = run(samples=("S", "T"))
check("6 two samples stop the script", isinstance(two, str) and "2 samples" in two, True)

# 7. a BAM without chr1 and chr9: no depth rows for RH, Duffy and ABO
part = run([vcf_line("JK_838", "1/1")], rhd(2), drop=("chr1", "chr9"))
check("7 no depth over RHD: RhD not determined", ph(part, "RhD"), "not determined")
check("7 the RhD note names the missing contig", "the BAM has no chr1 contig" in part["RhD"]["note"], True)
check("7 ABO with no record and no depth row: not determined", ph(part, "ABO"), "not determined")
check("7 the ABO note says why", "no depth row" in part["ABO"]["note"], True)
check("7 Duffy (chr1) not determined", ph(part, "Duffy"), "not determined")
check("7 Kidd (chr18, called) is read as usual", (ph(part, "Kidd"), cf(part, "Kidd")), ("Jk(a-b+)", "high"))
check("7 Kell (chr7, covered) is read as usual", (ph(part, "Kell"), cf(part, "Kell")),
      ("K-k+ Kp(a-b+) Js(a-b+)", "high"))
empty = run(depth={}, drop=tuple({r["chrom"] for r in SITES.values()}))
check("7 an empty depth file (no marker contig in the BAM): RhD not determined", ph(empty, "RhD"), "not determined")
check("7 and a system with no record is not determined", ph(empty, "Kidd"), "not determined")

# 8. the orientation of the null and HPA-7 to -35 rows: each row's REF and ALT are
# its c. change on the gene's strand (Ensembl GRCh38, 2026-11-03). A row with REF
# and ALT swapped, or on the wrong strand, fails here.
STRAND = {"SMIM1": 1, "ABCB6": -1, "ABCG2": -1, "FUT1": -1, "FUT2": 1, "KEL": -1, "SLC14A1": 1,
          "ITGB3": 1, "ITGA2B": -1, "GP1BB": 1, "ITGA2": 1, "GP9": 1}
COMP = str.maketrans("ACGT", "TGCA")


def oriented(gene, ref, alt, change):
    """The VCF REF/ALT agree with one c. change on the gene's strand: the bases
    of a substitution or a delins, the length of a deletion, duplication or insertion."""
    import re
    flip = (lambda s: s.translate(COMP)[::-1]) if STRAND[gene] == -1 else (lambda s: s)
    m = re.fullmatch(r"c\.[-*0-9+_]+([ACGT])>([ACGT])", change)
    if m:
        return (flip(ref), flip(alt)) == m.groups()
    m = re.fullmatch(r"c\.([-*0-9+]+)(?:_([-*0-9+]+))?delins([ACGT]+)", change)
    if m:
        return len(ref) == span(m.group(1), m.group(2)) and flip(alt) == m.group(3)
    m = re.fullmatch(r"c\.([-*0-9+]+)(?:_([-*0-9+]+))?(del|dup)([ACGT]*)", change)
    if m:
        n = len(m.group(4)) or span(m.group(1), m.group(2))
        d = len(ref) - len(alt) if m.group(3) == "del" else len(alt) - len(ref)
        return d == n and ref[0] == alt[0]
    m = re.fullmatch(r"c\.[0-9]+_[0-9]+ins([ACGT]+)", change)
    if m:
        return len(alt) - len(ref) == len(m.group(1)) and ref[0] == alt[0]
    return False


def span(a, b):
    def pos(x):
        base, _, off = x.partition("+")
        return int(base) + (int(off) if off else 0)
    return 1 if b is None else pos(b) - pos(a) + 1


NEW = sorted(NULLS | set(bt.HPA_RARE))
for n in NEW:
    r = SITES[n]
    gene, _, changes = r["hgvs"].partition(" ")
    changes = changes.split(", ")
    alts = r["alt"].split(",")
    check(f"8 {n}: one c. change per ALT", len(changes), len(alts))
    check(f"8 {n}: REF and ALT read the c. change on the {gene} strand",
          [c for a, c in zip(alts, changes) if not oriented(gene, r["ref"], a, c)], [])
check("8 control: a swapped REF and ALT fails the orientation", oriented("KEL", "A", "G", "c.578C>T"), False)
check("8 control: the right one passes", oriented("KEL", "G", "A", "c.578C>T"), True)
check("8 control: a deletion one base short fails", oriented("FUT2", "CG", "C", "c.719_720delGT"), False)


# 9. every new row, positive and negative: a reference genome carries the antigen
# (the a antigen for HPA-7 to -35), the row's ALT on both chromosomes lacks it.
def run_inproc(records):
    """Like run(), in this process (bin/bloodtype_report.py's main), with reads everywhere."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf, dp, out = os.path.join(tmp, "s.vcf"), os.path.join(tmp, "d.tsv"), os.path.join(tmp, "t.tsv")
        with open(vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + "".join(records))
        with open(dp, "w") as f:
            for r in list(SITES.values()) + list(REGIONS.values()):
                f.write("\t".join([r["chrom"], str(r["start"]), str(r["end"]), r["kind"], r["system"], r["name"],
                                   str(30 * (r["end"] - r["start"]))]) + "\n")
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            rc = bt.main(["--sample", "S", "--vcf", vcf, "--markers", MARKERS, "--depth", dp, "--out", out,
                          "--evidence", os.path.join(tmp, "e.tsv")])
        if rc:
            return f"exit {rc}"
        with open(out) as f:
            return [dict(r) for r in csv.DictReader(f, delimiter="\t")]


NEG = {"Vel": "Vel-", "Langereis": "Lan-", "Junior": "Jr(a-)", "H": "H- (Bombay, Oh)",
       "Secretor": "non-secretor", "Kell": "K-k- Kp(a-b-) Js(a-b-)", "Kidd": "Jk(a-b-)"}
POS = {"Vel": "Vel+", "Langereis": "Lan+", "Junior": "Jr(a+)", "H": "H+", "Secretor": "secretor",
       "Kell": "K-k+ Kp(a-b+) Js(a-b+)", "Kidd": "Jk(a+b-)"}
base = {r["system"]: r for r in run_inproc([])}
for n in NEW:
    r = SITES[n]
    s = r["system"]
    want_pos = POS.get(s) or f"{s}(a+b-)"
    check(f"9 {n}: the reference genome is {want_pos}", base[s]["phenotype"], want_pos)
    for a in r["alt"].split(","):
        recs = [vcf_line(n, "1/1", alt=a)]
        if s == "H":
            recs.append(vcf_line("FUT2_01N02", "1/1"))   # FUT1 null with a non-secretor: Bombay
        t = run_inproc(recs)
        got = {x["system"]: x for x in t} if isinstance(t, list) else {}
        want_neg = NEG.get(s) or f"{s}(a-b+)"
        check(f"9 {n} {r['ref']}>{a} on both chromosomes: {want_neg}",
              got.get(s, {}).get("phenotype"), want_neg)
check("9 a reference genome lacks no high-prevalence antigen: ABO is the first row",
      [x["system"] for x in run_inproc([])][:3], ["ABO", "RhD", "RhCE"])

# 10. synthetic gVCF records: Jk(a-b-) from SLC14A1 c.342-1G>A and HPA-21b/21b,
# each with the <*> allele DeepVariant's gVCF adds; a reference block at the
# sites gives the common phenotypes.
jk, hpa = SITES["JK_01N06"], SITES["HPA21_1960"]
gv = [f"{jk['chrom']}\t{jk['pos']}\t.\tG\tA,<*>\t40\tPASS\t.\tGT:GQ\t1/1:40\n",
      f"{hpa['chrom']}\t{hpa['pos']}\t.\t{hpa['ref']}\t{hpa['alt']},<*>\t40\tPASS\t.\tGT:GQ\t1/1:40\n"]
g = {x["system"]: x for x in run_inproc(gv)}
check("10 gVCF: c.342-1G>A 1/1 is Jk(a-b-), lacking Jk3",
      (g["Kidd"]["phenotype"], g["Kidd"]["genotype"], "lacks Jk3, carried by most people" in g["Kidd"]["note"]),
      ("Jk(a-b-)", "null/null", True))
check("10 gVCF: HPA-21b/21b", (g["HPA-21"]["phenotype"], g["HPA-21"]["genotype"]), ("HPA-21(a-b+)", "21b/21b"))
order = [x["system"] for x in run_inproc(gv)]
check("10 the systems that lack a high-prevalence antigen come first", order[:2], ["Kidd", "HPA-21"])
blocks = [f"{r['chrom']}\t{r['pos']}\t.\t{r['ref'][0]}\t<*>\t0\t.\tEND={r['pos'] + 5}\tGT:GQ\t0/0:50\n" for r in (jk, hpa)]
gb = {x["system"]: x for x in run_inproc(blocks)}
check("10 gVCF reference blocks: Jk(a+b-) and HPA-21(a+b-), confidence high",
      (gb["Kidd"]["phenotype"], gb["HPA-21"]["phenotype"], gb["Kidd"]["confidence"]),
      ("Jk(a+b-)", "HPA-21(a+b-)", "high"))

# 11. null alleles in detail
one = {x["system"]: x for x in run_inproc([vcf_line("VEL_01N01", "0/1")])}
check("11 one SMIM1 null copy: Vel+, Vel/null, the note names it",
      (one["Vel"]["phenotype"], one["Vel"]["genotype"], "SMIM1 c.64_80del" in one["Vel"]["note"]), ("Vel+", "Vel/null", True))
cpd = {x["system"]: x for x in run_inproc([vcf_line("LAN_01N02", "0/1"), vcf_line("LAN_01N13", "0/1")])}
check("11 two different ABCB6 nulls: Lan-, confidence low (no phase), the note says so",
      (cpd["Langereis"]["phenotype"], cpd["Langereis"]["confidence"], "phase is not read" in cpd["Langereis"]["note"]),
      ("Lan-", "low", True))
kk = {x["system"]: x for x in run_inproc([vcf_line("KEL_578", "0/1"), vcf_line("KEL_02N03", "0/1")])}
check("11 K/k with one KEL null: the null goes on the k chromosome, K+k-",
      (kk["Kell"]["phenotype"], kk["Kell"]["genotype"]), ("K+k- Kp(a-b+) Js(a-b+)", "K/null Kpb/null Jsb/null"))
check("11 and that lacks k, a high-prevalence antigen, so Kell comes first",
      (run_inproc([vcf_line("KEL_578", "0/1"), vcf_line("KEL_02N03", "0/1")])[0]["system"],
       "lacks k, carried by most people (high prevalence)" in kk["Kell"]["note"]), ("Kell", True))
jb = {x["system"]: x for x in run_inproc([vcf_line("JK_838", "0/1"), vcf_line("JK_02N06", "0/1")])}
check("11 Jka/Jkb with JK*02N.06 (c.871T>C): the null goes on Jkb, Jk(a+b-)",
      (jb["Kidd"]["phenotype"], jb["Kidd"]["genotype"]), ("Jk(a+b-)", "Jka/null"))
para = {x["system"]: x for x in run_inproc([vcf_line("FUT1_01N08", "1/1"), vcf_line("ABO_261", "0/1")])}
check("11 FUT1 null with a working FUT2: para-Bombay", para["H"]["phenotype"], "H- (para-Bombay)")
check("11 and the ABO row says the red cells group as O", "the cells usually group as O" in para["ABO"]["note"], True)
check("11 a non-secretor alone leaves H and ABO as they are",
      [(x["system"], x["phenotype"]) for x in run_inproc([vcf_line("FUT2_01N02", "1/1")])
       if x["system"] in ("H", "ABO", "Secretor")], [("ABO", "O"), ("H", "H+"), ("Secretor", "non-secretor")])
lowd = run(depth=dict(rhd(2), VEL_01N01=2.0))
check("11 the only Vel null row with 2 reads: not determined", ph(lowd, "Vel"), "not determined")
lowk = run(depth=dict(rhd(2), KEL_02N03=2.0))
check("11 one KEL null row with 2 reads: Kell still read, confidence medium",
      (ph(lowk, "Kell"), cf(lowk, "Kell")), ("K-k+ Kp(a-b+) Js(a-b+)", "medium"))
multi_null = {x["system"]: x for x in run_inproc([vcf_line("KEL_02N01", "1/1", alt="A")])}
check("11 a null row with several ALTs counts each (c.223+1G>T, the third)",
      multi_null["Kell"]["phenotype"], "K-k- Kp(a-b-) Js(a-b-)")


# 12. common variants are not null alleles. FUT2 c.772G>A (rs602662, about half of
# all chromosomes in gnomAD) and c.335C>T (rs200157007, 0.14 in South Asians) ride on
# the c.461G>A haplotype, and SLC14A1 c.810G>A (rs17675299, 0.08 in Finns) is far
# too common for a Jk null. As null rows they turned Se/se carriers into
# non-secretors and Jk(a+b+) into Jk(a+b-) at confidence high.
def raw_line(chrom, pos, ref, alt, gt):
    return "\t".join([chrom, str(pos), ".", ref, alt, "50", "PASS", ".", "GT", gt]) + "\n"


for desc, recs, system, want in [
        ("FUT2 c.461G>A 0/1 with c.772G>A 0/1: secretor",
         [vcf_line("FUT2_01N02", "0/1"), raw_line("chr19", 48703728, "G", "A", "0/1")], "Secretor", "secretor"),
        ("FUT2 c.772G>A 1/1 alone: secretor",
         [raw_line("chr19", 48703728, "G", "A", "1/1")], "Secretor", "secretor"),
        ("FUT2 c.461G>A 0/1 with c.335C>T 0/1: secretor",
         [vcf_line("FUT2_01N02", "0/1"), raw_line("chr19", 48703291, "C", "T", "0/1")], "Secretor", "secretor"),
        ("JK_838 0/1 with SLC14A1 c.810G>A 0/1: Jk(a+b+)",
         [vcf_line("JK_838", "0/1"), raw_line("chr18", 45739309, "G", "A", "0/1")], "Kidd", "Jk(a+b+)")]:
    got = {x["system"]: x for x in run_inproc(recs)}
    check(f"12 {desc}", got[system]["phenotype"], want)
check("12 those three positions are no null row",
      sorted(n for n in NULLS if (SITES[n]["chrom"], SITES[n]["pos"]) in
             {("chr19", 48703728), ("chr19", 48703291), ("chr18", 45739309)}), [])

# 13. another ALT at a null row's position (KEL c.1546C>T's position, G>T) is not that
# null allele: the antigens stay fully read at confidence high, and the note names it.
kel = SITES["KEL_02N17"]
oth = {x["system"]: x for x in run_inproc([raw_line(kel["chrom"], kel["pos"], "G", "T", "0/1")])}
check("13 another ALT at a KEL null row: K-k+, confidence high, the note names it",
      (oth["Kell"]["phenotype"], oth["Kell"]["confidence"], "another allele at this position" in oth["Kell"]["note"]),
      ("K-k+ Kp(a-b+) Js(a-b+)", "high", True))

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
