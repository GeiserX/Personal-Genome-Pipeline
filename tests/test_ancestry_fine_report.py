#!/usr/bin/env python3
"""bin/ancestry_fine_report.py, the Python half of step 26b, on made-up inputs:

  panel-map    groups from the population, related and grouped-out people
               are no reference, an excluded person is gone, the painting
               references are capped per group and taken in turn from each
               population, a person without metadata stops it
  order        the column of each group in a supervised ADMIXTURE .Q, and the
               errors for a K that is not the number of groups or two columns
               of one group
  sample-bed   the .bed code of each genotype (ALT is the .bim's A1), a site
               without a genotype or with another allele is left out of the
               .bed, the .bim and every level's .P.in alike, chr22 in the VCF
               meets 22 in the .bim, and the .ped allele pair of each code
  chroms       a chromosome is painted only with 100 SNVs spanning 5 cM
  report       fractions sorted with the bootstrap spread, neighbouring RFMix
               windows of one group merged into one BED row per haplotype,
               the painted shares, a real PNG of the right size, and the
               honest-limits text; and a level without a bootstrap has no range

Run: python3 tests/test_ancestry_fine_report.py
"""
import contextlib
import io
import os
import struct
import sys
import tempfile
import zlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import ancestry_fine_report as af  # noqa: E402

FAILS = []


def expect(desc, got, want):
    if got == want:
        print(f"PASS {desc}")
    else:
        print(f"FAIL {desc}: got {got!r}, want {want!r}")
        FAILS.append(desc)


def raises(desc, fn):
    try:
        with contextlib.redirect_stderr(io.StringIO()):
            fn()
    except (ValueError, SystemExit):
        print(f"PASS {desc}")
        return
    print(f"FAIL {desc}: no error")
    FAILS.append(desc)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


# --- panel-map ---------------------------------------------------------------------
def meta(s, pop, related="false"):
    return {"s": s, "population": pop, "sample_filters.all_samples_related": related}


rows = af.panel_map(
    [meta("y1", "YRI"), meta("y2", "YRI"), meta("y3", "YRI", "true"), meta("m1", "Mandenka"),
     meta("i1", "IBS"), meta("b1", "Basque"), meta("t1", "TSI"), meta("a1", "ASW"), meta("x1", "CHB")],
    ["y1", "y2", "y3", "m1", "i1", "b1", "t1", "a1", "x1"], exclude={"x1"}, paint_per_group=2)
by = {r[0]: r for r in rows}
expect("panel-map: an excluded person is not in the table", "x1" in by, False)
expect("panel-map: groups from the population (IBS: West Eurasia, Iberia)", by["i1"][2:4], ["WEA", "IBE"])
expect("panel-map: a Basque is in Iberia too", by["b1"][3], "IBE")
expect("panel-map: a related person is no reference", by["y3"][4], "0")
expect("panel-map: an admixed population (ASW) has no group and is no reference", by["a1"][2:5], ["-", "-", "0"])
expect("panel-map: at most 2 painting references per continental group",
       sum(r[5] == "1" for r in rows if r[2] == "AFR"), 2)
expect("panel-map: taken in turn from each population (Mandenka before a second YRI)",
       sorted(r[0] for r in rows if r[2] == "AFR" and r[5] == "1"), ["m1", "y1"])
raises("panel-map: a panel person without metadata stops it", lambda: af.panel_map([meta("y1", "YRI")], ["y1", "zz"]))

# --- order -------------------------------------------------------------------------
pops = ["EUR", "EUR", "AFR", "EAS", "AFR"]
q = [[1, 0, 0], [0.99, 0.01, 0], [0, 0, 1], [0, 1, 0], [0, 0.02, 0.98]]
expect("order: each column is the group whose people fill it", af.column_order(pops, q), ["EUR", "EAS", "AFR"])
raises("order: K other than the number of groups", lambda: af.column_order(pops, [r[:2] for r in q]))
raises("order: two columns of one group", lambda: af.column_order(["A", "A", "B"], [[1, 0], [1, 0], [1, 0]]))

# --- sample-bed --------------------------------------------------------------------
tmp = tempfile.mkdtemp()
bim = os.path.join(tmp, "ref.bim")
write(bim, "".join(f"22\t22:{p}\t0\t{p}\t{alt}\t{ref}\n" for p, ref, alt in
                   [(100, "A", "G"), (200, "C", "T"), (300, "G", "A"), (400, "T", "C"), (500, "A", "C")]))
vcf = os.path.join(tmp, "g.vcf")
write(vcf, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
      "chr22\t100\t.\tA\tG\t.\t.\t.\tGT\t1/1\n"     # hom ALT: code 0
      "chr22\t200\t.\tC\tT\t.\t.\t.\tGT\t0|1\n"     # het: code 2
      "chr22\t300\t.\tG\tA\t.\t.\t.\tGT\t0/0\n"     # hom REF: code 3
      "chr22\t400\t.\tT\tG\t.\t.\t.\tGT\t0/1\n")    # another ALT: left out; 500 missing: left out
keep, data = af.sample_bed(af.read_bim(bim), af.read_vcf_genotypes(vcf))
expect("sample-bed: the called SNVs (rows 0 to 2), chr22 meeting 22", keep, [0, 1, 2])
expect("sample-bed: .bed magic, then hom ALT 0, het 2, hom REF 3", data, b"\x6c\x1b\x01\x00\x02\x03")
p_file = os.path.join(tmp, "continent.P")
write(p_file, "".join(f"0.{i} 0.{9 - i}\n" for i in range(5)))
out = os.path.join(tmp, "adm")
with contextlib.redirect_stdout(io.StringIO()):
    af.main(["sample-bed", "--bim", bim, "--vcf", vcf, "--sample", "S", "--level", f"continent:2:{p_file}",
             "--bootstrap", "7", "--min-sites", "3", "--out-dir", out])
expect("sample-bed: the .P.in rows of the kept SNVs only",
       open(os.path.join(out, "S_continent.2.P.in")).read(), "0.0 0.9\n0.1 0.8\n0.2 0.7\n")
expect("sample-bed: the .bim rows of the kept SNVs only", len(open(os.path.join(out, "S_continent.bim")).readlines()), 3)
expect("sample-bed: .ped allele pairs (2 is A1): hom ALT, het, hom REF",
       open(os.path.join(out, "S.geno")).read(), "2 2\n1 2\n1 1\n")
sh = open(os.path.join(out, "run_admixture.sh")).read()
expect("sample-bed: run_admixture.sh runs 7 replicates of sample S", ("-le 7 ]" in sh, "S_${level}.bed" in sh), (True, True))
raises("sample-bed: fewer genotyped SNVs than --min-sites stops it", lambda: af.main(
    ["sample-bed", "--bim", bim, "--vcf", vcf, "--sample", "S", "--level", f"continent:2:{p_file}", "--out-dir", out]))

# --- chroms ------------------------------------------------------------------------
maps = os.path.join(tmp, "maps")
os.makedirs(maps)
write(os.path.join(maps, "chr20.map"), "chr20 . 0.0 1000000\nchr20 . 10.0 11000000\n")   # 1 cM per Mb
write(os.path.join(maps, "chr22.map"), "chr22 . 0.0 1000000\nchr22 . 10.0 11000000\n")
write(os.path.join(maps, "chr19.map"), "chr19 . 0.0 1000000\nchr19 . 10.0 11000000\n")
counts = [("chr20", 500, 2000000, 8000000), ("chr22", 500, 2000000, 4000000), ("chr21", 999, 1, 40000000),
          ("chr19", 50, 1000000, 11000000)]
expect("chroms: 6 cM painted; 2 cM, no map and 50 SNVs not", af.paint_chroms(counts, maps), ["chr20"])
expect("chroms: genetic position between two markers", af.cm_at(af.read_map(os.path.join(maps, "chr20.map")), 6000000), 5.0)

# --- report ------------------------------------------------------------------------
model = os.path.join(tmp, "model")
os.makedirs(model)
write(os.path.join(model, "levels.tsv"), "continent\t3\tEAS,WEA,AFR\nregion\t2\tIBE,ITA\n")
adm = os.path.join(tmp, "admixture")
os.makedirs(adm)
write(os.path.join(adm, "S_continent.3.Q"), "0.05 0.9 0.05\n")
write(os.path.join(adm, "S_continent.boot.Q"), "".join(f"0.05 {0.8 + i / 100:.2f} {0.15 - i / 100:.2f}\n" for i in range(11)))
write(os.path.join(adm, "S_region.2.Q"), "0.6 0.4\n")
paint = os.path.join(tmp, "paint")
os.makedirs(paint)
write(os.path.join(paint, "chr20.msp.tsv"),
      "#Subpopulation order/codes: AFR=0\tEAS=1\tWEA=2\n#chm\tspos\tepos\tsgpos\tegpos\tn snps\tS.0\tS.1\n"
      "chr20\t1000\t5000\t1\t2\t10\t2\t2\nchr20\t5000\t9000\t2\t3\t10\t2\t0\nchr20\t9000\t11000\t3\t4\t10\t2\t0\n")
rep = os.path.join(tmp, "out")
with contextlib.redirect_stdout(io.StringIO()):
    af.main(["report", "--sample", "S", "--model", model, "--admixture", adm, "--paint", paint,
             "--panel", "test", "--out", rep])
table = [line.rstrip("\n").split("\t") for line in open(os.path.join(rep, "S_ancestry_fine.tsv"))]
expect("report: the header", table[0], ["level", "code", "group", "fraction", "se", "low", "high"])
cont = [r for r in table if r[0] == "continent"]
expect("report: continental groups, largest first", [r[1] for r in cont], ["WEA", "EAS", "AFR"])
expect("report: West Eurasia's fraction and bootstrap range (2.5 to 97.5 percentile)", cont[0][3:], ["0.9000", "0.0332", "0.8025", "0.8975"])
expect("report: a level without bootstrap replicates has no range",
       [r[4:] for r in table if r[0] == "region"][0], [".", ".", "."])
bed = [line.rstrip("\n").split("\t") for line in open(os.path.join(rep, "S_painting.bed"))][1:]
expect("report: haplotype 1 one West Eurasian row, haplotype 2 split at 5000",
       bed, [["chr20", "999", "11000", "1", "WEA"], ["chr20", "999", "5000", "2", "WEA"],
             ["chr20", "4999", "11000", "2", "AFR"]])
paint_rows = {r[1]: r[3] for r in table if r[0] == "painting"}
expect("report: painted shares by length", paint_rows, {"WEA": "0.7000", "AFR": "0.3000"})
png = open(os.path.join(rep, "S_painting.png"), "rb").read()
w, h = struct.unpack(">II", png[16:24])
expect("report: a PNG of 1000 pixels by the 22 chromosomes", (png[:8], w, h), (b"\x89PNG\r\n\x1a\n", 1000, 6 + 22 * 19))
idat = png[png.index(b"IDAT") + 4:png.index(b"IEND") - 8]
raw = zlib.decompress(idat)
expect("report: the PNG's pixel data is whole", len(raw), h * (1 + 3 * w))
text = open(os.path.join(rep, "S_ancestry_fine.txt")).read()
expect("report: the text names the group, the fraction and its range",
       "West Eurasia (Europe, Middle East, North Africa)     90%  (80% to 90%)" in text, True)
expect("report: the text says a name is ours, not a nationality", "not a nationality" in text, True)
expect("report: and that Greek and Southern Italian overlap", "Greek and" in text and "Southern Italian ancestry" in text, True)
os.remove(os.path.join(paint, "chr20.msp.tsv"))
with contextlib.redirect_stdout(io.StringIO()):
    af.main(["report", "--sample", "S", "--model", model, "--admixture", adm, "--paint", paint,
             "--panel", "test", "--out", rep])
expect("report: without a painting the text says none",
       "Chromosome painting: none." in open(os.path.join(rep, "S_ancestry_fine.txt")).read(), True)

print("")
if FAILS:
    print(f"RESULT: {len(FAILS)} FAILED")
    sys.exit(1)
print("RESULT: ALL PASS")
