#!/usr/bin/env python3
"""bin/pgx_mito_merge.py: the MT-RNR1 outside call PharmCAT gets from the chrM
calls of step 20 / MITO_VARIANTS (step 07 and PHARMCAT run it).

The test builds a stand-in for pharmcat.jar, a zip holding
org/pharmgkb/pharmcat/phenotype/MT_RNR1.json with some of PharmCAT 3.4.0's
MT-RNR1 alleles and their function values, writes invented Mutect2-style chrM
VCFs (GT 0/1 at every AF, as Mutect2 writes them), and checks:
  1. a homoplasmic PASS m.1555A>G (AF 0.99) is the call; step 36's lines are
     copied before it and an MT-RNR1 line there is dropped;
  2. AF 0.95 counts as homoplasmic, 0.949 does not; --min-af moves the bound;
  3. a heteroplasmic, a filtered or an AF-less m.1555A>G passes nothing and
     says why; so do a homoplasmic PASS change at a named position that
     PharmCAT does not name and a PASS deletion near m.961; a filtered one,
     or a low-AF unnamed change, does not;
  4. PASS chrM records with no named allele give Reference when samtools
     depth shows 20 reads or more at every named position; without a depth
     file, with one position below 20 (--min-depth moves the bound) or absent
     from it, no call; a homoplasmic named allele needs no depth; a file with
     no PASS chrM record gives no call; records on other contigs are ignored;
  5. two homoplasmic named alleles: the increased-risk one is the call, the
     other is listed;
  6. INFO AF when FORMAT has none, the AF of the second ALT of a multi-allelic
     record, a bgzip/gzip input and the contig name MT are read;
  7. a jar without the MT-RNR1 definitions, or an unreadable VCF, exits 2.

Run: python3 tests/test_pgx_mito_merge.py
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "pgx_mito_merge.py")
FAILS = []

INC = "Increased risk of aminoglycoside-induced hearing loss"
UNC = "Uncertain risk of aminoglycoside-induced hearing loss"
NOR = "Normal risk of aminoglycoside-induced hearing loss"
ALLELES = [("Reference", NOR), ("m.827A>G", NOR), ("m.960C>del", UNC), ("m.961T>del", UNC),
           ("m.961T>del+Cn", UNC), ("m.1095T>C", INC), ("m.1189T>C", UNC), ("m.1555A>G", INC)]
HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


def check(desc, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS.append(desc)


def rec(pos, ref, alt, af, filt="PASS", chrom="chrM", info="."):
    fmt, val = ("GT:AF:DP", f"0/1:{af}:3000") if af is not None else ("GT:DP", "0/1:3000")
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filt}\t{info}\t{fmt}\t{val}\n"


BACKGROUND = rec(263, "A", "G", 0.998) + rec(750, "A", "G", 0.999)


class Case:
    def __init__(self, tmp):
        self.tmp = tmp
        self.jar = os.path.join(tmp, "pharmcat.jar")
        with zipfile.ZipFile(self.jar, "w") as z:
            z.writestr("org/pharmgkb/pharmcat/phenotype/MT_RNR1.json", json.dumps({
                "gene": "MT-RNR1",
                "namedAlleles": [{"name": n, "functionValue": f, "activityValue": None} for n, f in ALLELES]}))
        self.n = 0

    def run(self, records, outside=None, gz=False, extra=(), jar=None, depth=None):
        self.n += 1
        if depth is not None:
            path = os.path.join(self.tmp, f"d{self.n}.depth")
            with open(path, "w") as f:
                f.writelines(f"chrM\t{p}\t{depth.get(p, 500)}\n" for p in range(648, 1602) if depth.get(p, 500) is not None)
            extra = (*extra, "--depth", path)
        vcf = os.path.join(self.tmp, f"m{self.n}.vcf" + (".gz" if gz else ""))
        text = HEADER + records
        if gz:
            with gzip.open(vcf, "wt") as f:
                f.write(text)
        else:
            with open(vcf, "w") as f:
                f.write(text)
        out = os.path.join(self.tmp, f"o{self.n}.tsv")
        cmd = [sys.executable, SCRIPT, "--jar", jar or self.jar, "--mito", vcf, "--out", out, *extra]
        if outside is not None:
            path = os.path.join(self.tmp, f"c{self.n}.tsv")
            with open(path, "w") as f:
                f.write(outside)
            cmd += ["--outside", path]
        p = subprocess.run(cmd, capture_output=True, text=True)
        lines = open(out).read().splitlines() if os.path.exists(out) else None
        return p.returncode, p.stdout + p.stderr, lines

    def mt(self, lines):
        return [x.split("\t", 1)[1] for x in (lines or []) if x.startswith("MT-RNR1\t")]


def main():
    with tempfile.TemporaryDirectory() as tmp:
        c = Case(tmp)

        # 1. homoplasmic m.1555A>G
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.99),
                               outside="HLA-A\t*01:01/*26:01\nMT-RNR1\tm.1095T>C\n")
        check("1: exit 0", rc == 0, log)
        check("1: MT-RNR1 m.1555A>G is the call", c.mt(lines) == ["m.1555A>G"], lines)
        check("1: step 36's HLA-A line is copied first", (lines or [""])[0] == "HLA-A\t*01:01/*26:01", lines)
        check("1: the MT-RNR1 line of step 36's file is dropped", len(c.mt(lines)) == 1, lines)
        check("1: the log says it passed m.1555A>G, homoplasmic, with the risk",
              "MT-RNR1 m.1555A>G passed to PharmCAT as an outside call (homoplasmic, chrM:1555 A>G AF 0.990; "
              + INC in log, log)
        check("1: the log names the dropped line", "dropped the MT-RNR1 line" in log, log)

        # 2. the AF bound
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.95))
        check("2: AF 0.95 is homoplasmic", c.mt(lines) == ["m.1555A>G"], log)
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.949))
        check("2: AF 0.949 passes nothing", c.mt(lines) == [], log)
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.9), extra=("--min-af", "0.85"))
        check("2: --min-af 0.85 takes AF 0.9", c.mt(lines) == ["m.1555A>G"], log)

        # 3. what passes nothing
        for desc, records, want in [
            ("heteroplasmic", rec(1555, "A", "G", 0.4), "is heteroplasmic (AF below 0.95)"),
            ("filtered", rec(1555, "A", "G", 0.99, filt="possible_numt"), "is filtered (possible_numt)"),
            ("without an allele fraction", rec(1555, "A", "G", None), "has no allele fraction"),
            ("a change PharmCAT does not name", rec(1555, "A", "C", 0.99),
             "a change at an MT-RNR1 position that PharmCAT does not name"),
            ("a deletion near m.961", rec(955, "AC", "A", 0.99), "at or near m.960"),
            ("an insertion at m.961", rec(961, "T", "TC", 0.99), "(chrM:961 T>TC AF 0.990) at or near m.960"),
        ]:
            rc, log, lines = c.run(BACKGROUND + records, outside="HLA-B\t*35:01/*44:02\n")
            check(f"3: {desc}: exit 0", rc == 0, log)
            check(f"3: {desc}: no MT-RNR1 line", c.mt(lines) == [], lines)
            check(f"3: {desc}: step 36's line is still copied", lines == ["HLA-B\t*35:01/*44:02"], lines)
            check(f"3: {desc}: the log says why ({want})", want in log, log)
            check(f"3: {desc}: the log says the gene stays No Result", "it stays No Result" in log, log)
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "C", 0.3, filt="weak_evidence")
                               + rec(955, "AC", "A", 0.2, filt="weak_evidence"), depth={})
        check("3: a filtered unnamed change and a filtered deletion near m.961 hold nothing back",
              c.mt(lines) == ["Reference"], log)

        # 4. Reference, and no PASS record
        rc, log, lines = c.run(BACKGROUND + rec(2706, "A", "G", 0.999) + rec(1555, "A", "G", 0.99, chrom="chr1"),
                               depth={})
        check("4: no named allele on chrM, depth 500: Reference", c.mt(lines) == ["Reference"], log)
        check("4: the log says all 6 named positions are covered and counts the PASS records",
              "all 6 named positions have depth 20 or more (lowest 500); 3 PASS chrM records read" in log, log)
        rc, log, lines = c.run(BACKGROUND, depth={1555: 20, 961: 20})
        check("4: depth exactly 20 counts", c.mt(lines) == ["Reference"], log)
        rc, log, lines = c.run(BACKGROUND)
        check("4: no depth file: no call", c.mt(lines) == [], log)
        check("4: no depth file: the log says so", "there is no read depth" in log and "stays No Result" in log, log)
        rc, log, lines = c.run(BACKGROUND, depth={1555: 5})
        check("4: chrM:1555 at depth 5: no call", c.mt(lines) == [], log)
        check("4: chrM:1555 at depth 5: the log names it", "1 of its 6 positions have depth below 20 (chrM:1555 depth 5)" in log, log)
        rc, log, lines = c.run(BACKGROUND, depth={1555: 5}, extra=("--min-depth", "5"))
        check("4: --min-depth 5 takes depth 5", c.mt(lines) == ["Reference"], log)
        rc, log, lines = c.run(BACKGROUND, depth={960: None})
        check("4: a named position absent from the depth file counts as depth 0",
              c.mt(lines) == [] and "chrM:960 depth 0" in log, log)
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.99))
        check("4: a homoplasmic named allele needs no depth file", c.mt(lines) == ["m.1555A>G"], log)
        rc, log, lines = c.run(rec(263, "A", "G", 0.99, filt="weak_evidence"))
        check("4: no PASS chrM record: exit 0", rc == 0, log)
        check("4: no PASS chrM record: no call", c.mt(lines) == [], lines)
        check("4: no PASS chrM record: the log says so", "the chrM call set has no PASS record" in log, log)
        check("4: an empty output file is still written", lines == [], lines)

        # 5. two homoplasmic named alleles
        rc, log, lines = c.run(BACKGROUND + rec(1189, "T", "C", 0.99) + rec(1555, "A", "G", 0.98))
        check("5: the increased-risk allele is the call", c.mt(lines) == ["m.1555A>G"], log)
        check("5: the other is listed", "also homoplasmic, not passed" in log and "m.1189T>C" in log, log)

        # 6. AF sources, multi-allelic, gzip, MT
        rc, log, lines = c.run(BACKGROUND + "chrM\t1555\t.\tA\tG\t.\tPASS\tAF=0.97\tGT\t0/1\n")
        check("6: INFO AF is read when FORMAT has none", c.mt(lines) == ["m.1555A>G"], log)
        rc, log, lines = c.run(BACKGROUND + "chrM\t1555\t.\tA\tC,G\t.\tPASS\t.\tGT:AF\t0/1/2:0.01,0.98\n")
        check("6: the AF of the second ALT", c.mt(lines) == ["m.1555A>G"], log)
        rc, log, lines = c.run(BACKGROUND + rec(1555, "A", "G", 0.99), gz=True)
        check("6: a gzip VCF", c.mt(lines) == ["m.1555A>G"], log)
        rc, log, lines = c.run(rec(263, "A", "G", 0.99, chrom="MT") + rec(1555, "A", "G", 0.99, chrom="MT"))
        check("6: contig MT", c.mt(lines) == ["m.1555A>G"], log)

        # 7. errors
        bad = os.path.join(tmp, "empty.jar")
        with zipfile.ZipFile(bad, "w") as z:
            z.writestr("logback.xml", "<configuration/>")
        rc, log, lines = c.run(BACKGROUND, jar=bad)
        check("7: a jar without MT_RNR1.json exits 2", rc == 2, log)
        check("7: and names the resource", "org/pharmgkb/pharmcat/phenotype/MT_RNR1.json" in log, log)
        rc, log, lines = c.run("chrM\t1555\tA\n")
        check("7: a record with too few columns exits 2", rc == 2, log)

    if FAILS:
        print(f"\n{len(FAILS)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
