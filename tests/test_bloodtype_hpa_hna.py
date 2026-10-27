#!/usr/bin/env python3
"""bin/bloodtype_report.py: the platelet (HPA) and neutrophil (HNA) antigens of step 41.

Each case writes a synthetic one-sample VCF with planted genotypes at the HPA and
HNA sites of assets/bloodtype_markers.tsv, a samtools bedcov file with 30 reads
everywhere (so a site with no record is the reference allele, confidence high),
runs the script and checks the rows it writes:

  1. the table: every HPA and HNA system has its site(s), at the GRCh38 position
     and REF base Ensembl gave for its dbSNP id (FCGR3B's on FCGR3B, not FCGR3A)
  2. each two-allele antigen (HPA-1 to -6, -15, HNA-3, -4, -5) for 0/0, 0/1 and
     1/1 of its planted variant: phenotype and genotype
  3. the meaning: a person who lacks an antigen gets the table's meaning in the
     note (HPA-1b/1b: FNAIT and anti-HPA-1a), a person with both antigens none
  4. HNA-1: FCGR3B*01 and *02 from the four markers, *03 (HNA-1c, which carries
     HNA-1b too), c.233 with no call (*02 or *03), markers that disagree (not
     determined), c.194 that disagrees (a note), and the confidence is always
     low (FCGR3A, copy number)

Run: python3 tests/test_bloodtype_hpa_hna.py
"""
import csv
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


def rec(name, gt):
    r = SITES[name]
    return "\t".join([r["chrom"], str(r["pos"]), ".", r["ref"], r["alt"], "50", "PASS", ".", "GT", gt]) + "\n"


def run(records=()):
    """{system: row} of the table the script writes, or the error text."""
    with tempfile.TemporaryDirectory() as tmp:
        vcf, dp = os.path.join(tmp, "s.vcf"), os.path.join(tmp, "d.tsv")
        with open(vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                    + "".join(records))
        with open(dp, "w") as f:
            for r in list(SITES.values()) + list(REGIONS.values()):
                cols = [r["chrom"], str(r["start"]), str(r["end"]), r["kind"], r["system"], r["name"]]
                f.write("\t".join(cols + [str(30 * (r["end"] - r["start"]))]) + "\n")
        p = subprocess.run([sys.executable, SCRIPT, "--sample", "S", "--vcf", vcf, "--markers", MARKERS,
                            "--depth", dp, "--out", os.path.join(tmp, "t.tsv"),
                            "--evidence", os.path.join(tmp, "e.tsv")], capture_output=True, text=True)
        if p.returncode != 0:
            return p.stderr.strip() or f"exit {p.returncode}"
        with open(os.path.join(tmp, "t.tsv")) as f:
            return {r["system"]: r for r in csv.DictReader(f, delimiter="\t")}


def cell(t, system, col):
    if not isinstance(t, dict):
        return t
    return t[system][col] if system in t else f"no {system} row"


# 1. the table
want_sites = {
    # name: (system, chrom, GRCh38 position, REF, ALT, rsid)
    "HPA1_176": ("HPA-1", "chr17", 47283364, "T", "C", "rs5918"),
    "HPA2_482": ("HPA-2", "chr17", 4933086, "C", "T", "rs6065"),
    "HPA3_2621": ("HPA-3", "chr17", 44375697, "A", "C", "rs5911"),
    "HPA4_506": ("HPA-4", "chr17", 47284587, "G", "A", "rs5917"),
    "HPA5_1600": ("HPA-5", "chr5", 53062927, "G", "A", "rs1801106"),
    "HPA6_1544": ("HPA-6", "chr17", 47292422, "G", "A", "rs13306487"),
    "HPA15_2108": ("HPA-15", "chr6", 73783709, "A", "C", "rs10455097"),
    "FCGR3B_108": ("HNA-1", "chr1", 161629989, "G", "C", "rs200688856"),
    "FCGR3B_114": ("HNA-1", "chr1", 161629983, "A", "G", "rs527909462"),
    "FCGR3B_194": ("HNA-1", "chr1", 161629903, "T", "C", "rs448740"),
    "FCGR3B_233": ("HNA-1", "chr1", 161629864, "G", "T", "rs5030738"),
    "FCGR3B_244": ("HNA-1", "chr1", 161629853, "T", "C", "rs147574249"),
    "FCGR3B_316": ("HNA-1", "chr1", 161629781, "T", "C", "rs2290834"),
    "HNA3_461": ("HNA-3", "chr19", 10631494, "A", "G", "rs2288904"),
    "HNA4_230": ("HNA-4", "chr16", 31265490, "G", "A", "rs1143679"),
    "HNA5_2372": ("HNA-5", "chr16", 30506720, "G", "C", "rs2230433"),
}
for n, w in want_sites.items():
    r = SITES.get(n)
    got = (r["system"], r["chrom"], r["pos"], r["ref"], r["alt"], r["rsid"]) if r else None
    check(f"1 {n} at its GRCh38 position", got, w)
systems = [s for s, _ in bt.SYSTEMS]
for s in ("HPA-6", "HNA-1", "HNA-3", "HNA-4", "HNA-5"):
    check(f"1 {s} is a system of the table", s in systems, True)
check("1 every HPA and HNA site has a meaning",
      sorted(n for n, r in SITES.items() if r["system"][:4] in ("HPA-", "HNA-") and not r.get("meaning")), [])
check("1 no red cell site has a meaning",
      sorted(n for n, r in SITES.items() if r["system"][:4] not in ("HPA-", "HNA-") and r.get("meaning")), [])
check("1 every FCGR3B site carries the FCGR3A and copy-number caveat",
      sorted(n for n in want_sites if n.startswith("FCGR3B") and "copy number" not in SITES.get(n, {}).get("caveat", "")), [])

# 2. the two-allele antigens, 0/0, 0/1, 1/1
pairs = {
    # name: (system, phenotype at 0/0, 0/1, 1/1, genotype at 0/1)
    "HPA1_176": ("HPA-1", "HPA-1(a+b-)", "HPA-1(a+b+)", "HPA-1(a-b+)", "1a/1b"),
    "HPA2_482": ("HPA-2", "HPA-2(a+b-)", "HPA-2(a+b+)", "HPA-2(a-b+)", "2a/2b"),
    "HPA3_2621": ("HPA-3", "HPA-3(a+b-)", "HPA-3(a+b+)", "HPA-3(a-b+)", "3a/3b"),
    "HPA4_506": ("HPA-4", "HPA-4(a+b-)", "HPA-4(a+b+)", "HPA-4(a-b+)", "4a/4b"),
    "HPA5_1600": ("HPA-5", "HPA-5(a+b-)", "HPA-5(a+b+)", "HPA-5(a-b+)", "5a/5b"),
    "HPA6_1544": ("HPA-6", "HPA-6(a+b-)", "HPA-6(a+b+)", "HPA-6(a-b+)", "6a/6b"),
    "HPA15_2108": ("HPA-15", "HPA-15(a-b+)", "HPA-15(a+b+)", "HPA-15(a+b-)", "15a/15b"),
    "HNA3_461": ("HNA-3", "HNA-3(a-b+)", "HNA-3(a+b+)", "HNA-3(a+b-)", "3a/3b"),
    "HNA4_230": ("HNA-4", "HNA-4(a+b-)", "HNA-4(a+b+)", "HNA-4(a-b+)", "4a/4b"),
    "HNA5_2372": ("HNA-5", "HNA-5(a+b-)", "HNA-5(a+b+)", "HNA-5(a-b+)", "5a/5b"),
}
none = run()
for n, (system, p0, p1, p2, g1) in pairs.items():
    het, hom = run([rec(n, "0/1")]), run([rec(n, "1/1")])
    check(f"2 {system}: 0/0, 0/1, 1/1", (cell(none, system, "phenotype"), cell(het, system, "phenotype"),
                                         cell(hom, system, "phenotype")), (p0, p1, p2))
    check(f"2 {system}: the 0/1 genotype at confidence high",
          (cell(het, system, "genotype"), cell(het, system, "confidence")), (g1, "high"))
    check(f"2 {system}: both antigens, no meaning in the note", "lacks " in cell(het, system, "note"), False)

# 3. the meaning
mum = run([rec("HPA1_176", "1/1")])
note = cell(mum, "HPA-1", "note")
check("3 HPA-1b/1b: the note says the person lacks HPA-1a", note.startswith("lacks HPA-1a: "), True)
check("3 and names FNAIT and anti-HPA-1a", "FNAIT" in note and "anti-HPA-1a" in note, True)
check("3 HPA-1a/1a lacks HPA-1b", cell(none, "HPA-1", "note").startswith("lacks HPA-1b: "), True)
check("3 GRCh38 carries HNA-3b: no record lacks HNA-3a, with the neutrophil meaning",
      (cell(none, "HNA-3", "note").startswith("lacks HNA-3a: "), "neonatal alloimmune neutropenia"
       in cell(none, "HNA-3", "note")), (True, True))
check("3 a red cell system never gets a meaning", "lacks " in cell(mum, "Kidd", "note"), False)

# 4. HNA-1
A4 = ["FCGR3B_108", "FCGR3B_114", "FCGR3B_244", "FCGR3B_316"]
b = none
check("4 the GRCh38 bases (FCGR3B*02 at four markers): HNA-1(a-b+c-), FCGR3B*02 twice",
      (cell(b, "HNA-1", "phenotype"), cell(b, "HNA-1", "genotype")), ("HNA-1(a-b+c-)", "FCGR3B*02/FCGR3B*02"))
check("4 confidence is always low (FCGR3A, copy number)", cell(b, "HNA-1", "confidence"), "low")
check("4 c.194 that does not match is named in the note",
      "FCGR3B c.194A>G (x0) does not match" in cell(b, "HNA-1", "note"), True)
check("4 the note lacks HNA-1a and HNA-1c", "lacks HNA-1a, HNA-1c: " in cell(b, "HNA-1", "note"), True)
ab = run([rec(n, "0/1") for n in A4] + [rec("FCGR3B_194", "0/1")])
check("4 HNA-1a/1b", (cell(ab, "HNA-1", "phenotype"), cell(ab, "HNA-1", "genotype")),
      ("HNA-1(a+b+c-)", "FCGR3B*01/FCGR3B*02"))
check("4 c.194 that matches: no note about it", "does not match" in cell(ab, "HNA-1", "note"), False)
aa = run([rec(n, "1/1") for n in A4] + [rec("FCGR3B_194", "1/1")])
check("4 HNA-1a/1a", (cell(aa, "HNA-1", "phenotype"), cell(aa, "HNA-1", "genotype")),
      ("HNA-1(a+b-c-)", "FCGR3B*01/FCGR3B*01"))
check("4 HNA-1a/1a lacks HNA-1b", "lacks HNA-1b, HNA-1c: " in cell(aa, "HNA-1", "note"), True)
c = run([rec("FCGR3B_233", "0/1")] + [rec("FCGR3B_194", "1/1")])
check("4 HNA-1c on an FCGR3B*02: HNA-1(a-b+c+), FCGR3B*02/FCGR3B*03",
      (cell(c, "HNA-1", "phenotype"), cell(c, "HNA-1", "genotype")), ("HNA-1(a-b+c+)", "FCGR3B*02/FCGR3B*03"))
ac = run([rec(n, "1/1") for n in A4] + [rec("FCGR3B_194", "1/1"), rec("FCGR3B_233", "0/1")])
check("4 HNA-1c with two FCGR3B*01: b+ (FCGR3B*03 carries HNA-1b) and the extra-copy note",
      (cell(ac, "HNA-1", "phenotype"), "extra FCGR3B copy" in cell(ac, "HNA-1", "note")), ("HNA-1(a+b+c+)", True))
nc = run([rec("FCGR3B_233", "./."), rec("FCGR3B_194", "1/1")])
check("4 c.233 with no call: the HNA-1b alleles read *02 or *03, HNA-1c not determined",
      (cell(nc, "HNA-1", "phenotype"), cell(nc, "HNA-1", "genotype"), "HNA-1c not determined" in cell(nc, "HNA-1", "note")),
      ("HNA-1(a-b+)", "FCGR3B*02 or *03/FCGR3B*02 or *03", True))
mixed = run([rec("FCGR3B_108", "1/1"), rec("FCGR3B_114", "0/1")])
check("4 markers that disagree: not determined", cell(mixed, "HNA-1", "phenotype"), "not determined")
check("4 and the note says why", "do not agree on HNA-1a or HNA-1b" in cell(mixed, "HNA-1", "note"), True)

if FAILS:
    print(f"{len(FAILS)} check(s) failed")
    sys.exit(1)
print("all checks passed")
