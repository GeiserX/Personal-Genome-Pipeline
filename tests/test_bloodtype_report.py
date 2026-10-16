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
ruled = set(bt.PAIRS) | bt.SPECIAL_SITES
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

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
