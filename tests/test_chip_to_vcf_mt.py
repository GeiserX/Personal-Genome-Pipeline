#!/usr/bin/env python3
"""scripts/chip-to-vcf.sh: MyHeritage input and the mitochondrion.

Runs the script with the real bcftools and Picard images on a made-up
reference: 3 kb contigs 1, X and Y, and a 16,569-base mitochondrion that is
the same sequence in the "hg19" FASTA (MT, the rCRS of human_g1k_v37) and in
the "GRCh38" one (chrM). The chain maps 1, X and Y one to one, and chrM the
way UCSC's hg19ToHg38 chain does: hg19's chrM is 16,571 bases, so past an
early gap every position moves two bases (1555 to 1553 here).

Checks:
  1. tests/fixtures/chip/myheritage.csv (quoted CSV, a column header row,
     CRLF line endings): the script exits 0, the TSV it hands bcftools has
     ten rows and no carriage return, and the VCF has ten records, the
     no-call as a missing genotype, every genotype a plain call, on chr1,
     chrM, chrX and chrY. The script before this fix left the carriage
     return on the genotype column of that TSV.
  2. chrM keeps its rCRS positions: m.73, m.1555 and m.3243 sit at 73, 1555
     and 3243 with the REF of the rCRS there, and nothing sits at 1553 or
     3241. The script before this fix lifted chrM with the chain: m.1555
     landed on 1553 (the same base there, so Picard kept it) and m.3243 went
     to the reject file.
  3. The same rows as a 23andMe file, its MT sites numbered on hg19's chrM
     as 23andMe numbers them (1557 and 3245 under this chain, plus 301
     inside the chain's gap), give the same chrM records: each MT position
     moves to the rCRS through the chrM chain and the gap site is dropped
     with its count printed. Without that move 23andMe's
     sites would stay at chrM:1557 and 3245, with no call at 1555 or 3243.

Needs docker. Run: python3 tests/test_chip_to_vcf_mt.py
"""
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "scripts", "chip-to-vcf.sh")
FIXTURE = os.path.join(REPO, "tests", "fixtures", "chip", "myheritage.csv")
MT_LEN = 16569
SHORT = 3000
# rCRS-like bases the checks rely on: 1553 and 1555 the same base (so the
# shifted record passes Picard's REF check), 3241 and 3243 different ones.
FIXED = {73: "A", 1553: "A", 1555: "A", 3241: "G", 3243: "A"}
# 23andMe's MT numbering under the test chain: rCRS position to hg19 chrM
# (two bases more past the gap at 301-302, as m.3243 is 3244 or 3245 there).
HG19_MT = {1555: 1557, 3243: 3245}

FAILS = []


def check(desc, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + desc + ("" if ok else "  " + detail))
    if not ok:
        FAILS.append(desc)


def image(name):
    out = subprocess.run(["bash", "-c", '. "$1/versions.env"; printf %s "${!2}"', "-", REPO, name],
                         capture_output=True, text=True, check=True).stdout
    if not out:
        sys.exit("FAIL: %s is not set in versions.env" % name)
    return out


def sequence(n, seed, fixed=None):
    s, bases = seed, "ACGT"
    out = []
    for _ in range(n):
        s = (s * 75 + 74) % 65537
        out.append(bases[s % 4])
    for pos, b in (fixed or {}).items():
        out[pos - 1] = b
    return "".join(out)


def write_fasta(path, contigs):
    off = 0
    with open(path, "w") as fa, open(path + ".fai", "w") as fai:
        for name, seq in contigs:
            hdr = ">%s\n" % name
            fa.write(hdr)
            off += len(hdr)
            lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
            fa.write("".join(l + "\n" for l in lines))
            fai.write("%s\t%d\t%d\t60\t61\n" % (name, len(seq), off))
            off += len(seq) + len(lines)


def setup_genome(g):
    for sub in ("reference_hg19", "reference", "liftover"):
        os.makedirs(os.path.join(g, sub))
    nuclear = [sequence(SHORT, 7 + i) for i in range(3)]
    mt = sequence(MT_LEN, 11, FIXED)
    write_fasta(os.path.join(g, "reference_hg19", "human_g1k_v37.fasta"),
                list(zip(["1", "X", "Y"], nuclear)) + [("MT", mt)])
    write_fasta(os.path.join(g, "reference", "test38.fasta"),
                list(zip(["chr1", "chrX", "chrY"], nuclear)) + [("chrM", mt)])
    with open(os.path.join(g, "reference", "test38.dict"), "w") as d:
        d.write("@HD\tVN:1.6\tSO:unsorted\n")
        for c, n in (("chr1", SHORT), ("chrX", SHORT), ("chrY", SHORT), ("chrM", MT_LEN)):
            d.write("@SQ\tSN:%s\tLN:%d\n" % (c, n))
    chain = ""
    for c in ("chr1", "chrX", "chrY"):
        chain += "chain 3000 %s 3000 + 0 3000 %s 3000 + 0 3000 1\n3000\n\n" % (c, c)
    # hg19's chrM has two more bases than the rCRS; the gap after 300 bases
    # moves every later position by two, as UCSC's chain moves m.1555 to 1553.
    chain += ("chain 16569 chrM 16571 + 0 16571 chrM 16569 + 0 16569 2\n"
              "300\t2\t0\n16269\n\n")
    with gzip.open(os.path.join(g, "liftover", "hg19ToHg38.over.chain.gz"), "wt") as f:
        f.write(chain)
    return mt


def run(g, fmt):
    env = dict(os.environ, GENOME_DIR=g, REF_FASTA=os.path.join(g, "reference", "test38.fasta"))
    p = subprocess.run([SCRIPT, "S1", fmt], env=env, capture_output=True, text=True)
    log = p.stdout + p.stderr
    print("\n".join("    " + l for l in log.splitlines()))
    return p.returncode, log


def records(g, img):
    p = subprocess.run(["docker", "run", "--rm", "-v", "%s:/g:ro" % g, img, "bcftools", "query",
                        "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n", "/g/S1/vcf/S1.vcf.gz"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        return None, p.stderr
    return [l.split("\t") for l in p.stdout.splitlines()], ""


def check_mt(name, rows, mt):
    mrows = {int(r[1]): r for r in rows if r[0] == "chrM"}
    for pos in (73, 1555, 3243):
        r = mrows.get(pos)
        check("%s: m.%d stays at chrM:%d with the rCRS base %s as REF" % (name, pos, pos, mt[pos - 1]),
              r is not None and r[2] == mt[pos - 1], str(sorted(mrows)))
    check("%s: nothing moved to chrM:1553 or chrM:3241" % name,
          not ({1553, 3241} & set(mrows)), str(sorted(mrows)))
    r = mrows.get(1555)
    check("%s: m.1555 A>G is the alt call" % name,
          r is not None and (r[3], r[4]) == ("G", "1"), str(r))


def main():
    if shutil.which("docker") is None:
        sys.exit("FAIL: docker is needed to run this test")
    img = image("BCFTOOLS_IMAGE")
    with tempfile.TemporaryDirectory() as work:
        # 1 and 2: the MyHeritage fixture
        g = os.path.join(work, "myheritage")
        mt = setup_genome(g)
        os.makedirs(os.path.join(g, "S1", "raw"))
        shutil.copy(FIXTURE, os.path.join(g, "S1", "raw", "MyHeritage_raw_dna_data.csv"))
        with open(FIXTURE, "rb") as f:
            check("the MyHeritage fixture has CRLF line endings", b"\r\n" in f.read())
        print("--- myheritage: chip-to-vcf.sh S1 auto")
        rc, log = run(g, "auto")
        check("myheritage: auto-detected", "Auto-detected format: myheritage" in log, log[-400:])
        check("myheritage: chip-to-vcf.sh exits 0 (exit %d)" % rc, rc == 0, log[-400:])
        tsv = os.path.join(g, "S1", "raw", "S1_tsv2vcf.tsv")
        with open(tsv, "rb") as f:
            body = f.read() if os.path.exists(tsv) else b""
        check("myheritage: the converted TSV has ten rows and no carriage return",
              body.count(b"\n") == 10 and b"\r" not in body, repr(body[:200]))
        rows, err = records(g, img) if rc == 0 else (None, "no run")
        if rows is None:
            check("myheritage: the VCF reads back", False, err)
        else:
            check("myheritage: 10 records", len(rows) == 10, str(rows))
            gts = sorted({r[4] for r in rows})
            check("myheritage: every genotype is a plain call, no carriage return",
                  all(re.fullmatch(r"(\.|\d+)([/|](\.|\d+))?", gt) for gt in gts), repr(gts))
            nocall = [r[4] for r in rows if (r[0], r[1]) == ("chr1", "404")]
            check("myheritage: the no-call row is a missing genotype", nocall in (["./."], ["."]), str(nocall))
            chroms = sorted({r[0] for r in rows})
            check("myheritage: chromosomes chr1 chrM chrX chrY", chroms == ["chr1", "chrM", "chrX", "chrY"], str(chroms))
            check_mt("myheritage", rows, mt)

        # 3: the same rows as a 23andMe file
        g = os.path.join(work, "23andme")
        mt = setup_genome(g)
        os.makedirs(os.path.join(g, "S1", "raw"))
        with open(FIXTURE) as f, open(os.path.join(g, "S1", "raw", "S1_raw.txt"), "w") as o:
            o.write("# rsid\tchromosome\tposition\tgenotype\n")
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and not line.startswith("RSID"):
                    f4 = [x.strip('"') for x in line.split(",")]
                    if f4[1] == "MT":
                        # rCRS to hg19 chrM numbering under the test chain
                        f4[2] = str(HG19_MT.get(int(f4[2]), int(f4[2])))
                    o.write("\t".join(f4) + "\n")
            # a site inside the chain's gap: no rCRS position, dropped
            o.write("rs9000099\tMT\t301\tA\n")
        print("--- 23andme: chip-to-vcf.sh S1 23andme")
        rc, log = run(g, "23andme")
        check("23andme: chip-to-vcf.sh exits 0 (exit %d)" % rc, rc == 0, log[-400:])
        check("23andme: the MT move and the gap drop are printed",
              "3 site(s) moved from hg19 chrM to rCRS positions, 1 in a chain gap dropped" in log, log[-600:])
        rows, err = records(g, img) if rc == 0 else (None, "no run")
        if rows is None:
            check("23andme: the VCF reads back", False, err)
        else:
            check_mt("23andme", rows, mt)
            mpos = sorted(int(r[1]) for r in rows if r[0] == "chrM")
            check("23andme: chrM holds 73, 1555 and 3243 only (nothing left at 301, 1557 or 3245)",
                  mpos == [73, 1555, 3243], str(mpos))

    if FAILS:
        print("%d check(s) failed" % len(FAILS))
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
