#!/usr/bin/env python3
"""bin/archaic_report.py (step 43) on a synthetic archaic panel, with the real IBDmix image.

Checks:
  1. table: a 1000 Genomes site keeps its frequencies and gets each archaic
     genome's ALT count; an archaic record with another ALT or another REF is
     unknown (9); a site only an archaic genome varies at gets frequency 0; a
     site no archaic genome was called at, or only called as reference with
     no modern variation, is left out; a position given twice by 1000 Genomes
     is left out;
  2. genotypes: the stand-in columns carry the population's ALT count; your
     column is your ALT count, 0 where your VCF has no record, 9 for a
     non-PASS record; a site where your VCF has another ALT is left out; a
     chromosome your VCF has no record on is not analysed;
  3. report: Altai and Vindija segments merge into Neanderthal ones, a
     Denisova segment that overlaps one is not counted as Denisovan, segments
     under 50 kb or LOD 4 are dropped, the share is of the analysed span, and
     the genes of the longest segments are named; bin/collect_summary.py
     reads the summary and stops on a damaged one, and bin/render_report.py
     prints both shares;
  4. with IBDmix itself (IBDMIX_IMAGE, docker): a planted 600 kb stretch where
     the synthetic person carries the Neanderthal haplotype comes out as one
     Neanderthal segment over it, and a planted Denisovan one as Denisovan;
     the same person without the planted stretches gets no segment at all.
     Set ARCHAIC_TEST_NO_DOCKER=1 to skip part 4 on a machine without docker;
     CI never sets it.

Run: python3 tests/test_archaic_report.py
"""
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import archaic_report as ar  # noqa: E402

SCRIPT = os.path.join(REPO, "bin", "archaic_report.py")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def write(path, text):
    with open(path, "w") as f:
        f.write(text)


def run(*args):
    return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)


def image(var):
    for line in open(os.path.join(REPO, "versions.env")):
        if line.startswith(f"{var}="):
            return line.split("=", 1)[1].split("#")[0].strip().strip('"')
    raise SystemExit(f"{var} not in versions.env")


TMP = tempfile.mkdtemp()

# --- 1. table ------------------------------------------------------------------
afs = "0.1\t0.2\t0.05\t0\t0.15\t0.1"
write(f"{TMP}/modern.tsv",
      f"100\tA\tG\t{afs}\n"        # every archaic genome called
      f"200\tC\tT\t{afs}\n"        # Vindija has another ALT, Denisova another REF
      f"300\tG\tA\t{afs}\n"        # no archaic record
      f"400\tT\tC\t{afs}\n400\tT\tG\t{afs}\n")  # a split multi-allelic site
write(f"{TMP}/altai.tsv", "100\tA\tG\t1/1\n200\tC\t.\t0/0\n500\tG\tC\t0/1\n600\tA\t.\t0/0\n")
write(f"{TMP}/vindija.tsv", "100\tA\t.\t0/0\n200\tC\tA\t0/1\n500\tG\t.\t0/0\n600\tA\t.\t0/0\n")
write(f"{TMP}/denisova.tsv", "100\tA\tG\t0|1\n200\tT\tC\t0/1\n400\tT\tC\t1/1\n500\tG\tT\t1/1\n")
r = run("table", "--chrom", "22", "--modern", f"{TMP}/modern.tsv", "--archaic", f"Altai={TMP}/altai.tsv",
        "--archaic", f"Vindija={TMP}/vindija.tsv", "--archaic", f"Denisova={TMP}/denisova.tsv", "--out", f"{TMP}/t.vcf")
check("table exits 0", r.returncode == 0, r.stderr)
recs = {int(l.split("\t")[1]): l.rstrip("\n").split("\t") for l in open(f"{TMP}/t.vcf") if not l.startswith("#")}
info = {p: dict(kv.split("=") for kv in c[7].split(";")) for p, c in recs.items()}
check("rows kept: the modern sites with an archaic call and the archaic-only site", sorted(recs) == [100, 200, 500], sorted(recs))
check("a 1000 Genomes site: each genome's ALT count", (info[100]["Altai"], info[100]["Vindija"], info[100]["Denisova"]) == ("2", "0", "1"), info.get(100))
check("its frequencies are kept", info[100]["AF_ALL"] == "0.1" and info[100]["AF_EUR"] == "0.15", info.get(100))
check("another ALT or another REF is unknown (9)", (info[200]["Altai"], info[200]["Vindija"], info[200]["Denisova"]) == ("0", "9", "9"), info.get(200))
check("an archaic-only site takes the first non-reference ALT and frequency 0",
      recs[500][4] == "C" and info[500]["Altai"] == "1" and info[500]["Vindija"] == "0" and info[500]["Denisova"] == "9"
      and info[500]["AF_ALL"] == "0", (recs.get(500), info.get(500)))
check("chrom gets its chr prefix for the liftover", recs[100][0] == "chr22", recs[100][0])

# --- 2. genotypes ------------------------------------------------------------------
tdir = f"{TMP}/table"
os.makedirs(tdir)
hdr = "\t".join(ar.TABLE_HEADER) + "\n"
rows = ["chr22\t1000\tA\tG\t2\t2\t0\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5",
        "chr22\t2000\tC\tT\t0\t0\t0\t0.01\t0\t0\t0\t0.02\t0",
        "chr22\t3000\tG\tA\t1\t1\t2\t0\t0\t0\t0\t0\t0",
        "chr22\t4000\tT\tC\t2\t9\t0\t0.2\t0\t0\t0\t0\t0",
        "chr22\t5000\tA\tC\t0\t0\t0\t0.2\t0\t0\t0\t0\t0"]
with gzip.open(f"{tdir}/chr22.tsv.gz", "wt") as f:
    f.write(hdr + "\n".join(rows) + "\n")
with gzip.open(f"{tdir}/chr21.tsv.gz", "wt") as f:
    f.write(hdr + "chr21\t1000\tA\tG\t2\t2\t0\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\n")
vcf = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
       "chr22\t1000\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/1:30\n"
       "chr22\t3000\t.\tG\tA\t50\tRefCall\t.\tGT\t1/1\n"
       "chr22\t5000\t.\tA\tT\t50\tPASS\t.\tGT\t1/1\n")
write(f"{TMP}/s.vcf", vcf)
r = run("genotypes", "--table-dir", tdir, "--vcf", f"{TMP}/s.vcf", "--pop", "EUR", "--pseudo", "10", "--out-dir", f"{TMP}/gt")
check("genotypes exits 0", r.returncode == 0, r.stderr)
lines = gzip.open(f"{TMP}/gt/chr22.gt.gz", "rt").read().splitlines()
g = {int(l.split("\t")[1]): l.split("\t") for l in lines[1:]}
check("header: four columns, three archaic, SAMPLE, ten stand-ins",
      lines[0].split("\t")[:8] == ["chrom", "pos", "ref", "alt", "Altai", "Vindija", "Denisova", "SAMPLE"]
      and len(lines[0].split("\t")) == 18, lines[0])
check("the chromosome is an integer for IBDmix", g[1000][0] == "22", g[1000][0])
check("your heterozygous call is 1", g[1000][7] == "1", g[1000])
check("no record is 0 (reference)", g[2000][7] == "0", g[2000])
check("a non-PASS record is 9", g[3000][7] == "9", g[3000])
check("a site where you carry another ALT is left out", 5000 not in g, sorted(g))
check("stand-ins carry EUR's ALT count (0.5 of 20 alleles = 10)", sum(int(x) for x in g[1000][8:]) == 10, g[1000][8:])
check("and 0.02 of 20 rounds to 0", sum(int(x) for x in g[2000][8:]) == 0, g[2000][8:])
check("each stand-in column is one character", all(len(x) == 1 for x in g[1000][4:]), g[1000])
check("a chromosome the VCF has no record on is not analysed", not os.path.exists(f"{TMP}/gt/chr21.gt.gz"))
spans = open(f"{TMP}/gt/spans.tsv").read().splitlines()
check("spans: chr22 from 1000 to 4000, four sites", spans[1].split("\t")[:4] == ["chr22", "1000", "4000", "4"], spans)
samples = open(f"{TMP}/gt/samples.txt").read().split()
check("samples.txt: SAMPLE and the stand-ins", samples[0] == "SAMPLE" and len(samples) == 11, samples)
write(f"{TMP}/two.vcf", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n")
r = run("genotypes", "--table-dir", tdir, "--vcf", f"{TMP}/two.vcf", "--out-dir", f"{TMP}/gt2")
check("a two-sample VCF stops the step", r.returncode != 0 and "one-sample" in r.stderr, r.stderr)

# --- 3. report ---------------------------------------------------------------------
idir = f"{TMP}/ibd"
os.makedirs(idir)
head = "ID\tchrom\tstart\tend\tslod\n"
write(f"{idir}/chr22.Altai.txt", head + "SAMPLE\t22\t1000000\t1300000\t40\nP001\t22\t5000000\t9000000\t99\n"
      "SAMPLE\t22\t2000000\t2040000\t30\nSAMPLE\t22\t3000000\t3100000\t3.5\n")
write(f"{idir}/chr22.Vindija.txt", head + "SAMPLE\t22\t1200000\t1500000\t50\n")
write(f"{idir}/chr22.Denisova.txt", head + "SAMPLE\t22\t1400000\t1450000\t12\nSAMPLE\t22\t6000000\t6200000\t20\n")
write(f"{TMP}/spans.tsv", "chrom\tfirst\tlast\tsites\nchr22\t1\t10000000\t5000\n")
write(f"{TMP}/genes.gtf",
      'chr22\tHAVANA\tgene\t1100000\t1110000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "GENEA";\n'
      'chr22\tHAVANA\tgene\t1120000\t1130000\t.\t+\t.\tgene_id "G2"; gene_type "lncRNA"; gene_name "LNC1";\n'
      'chr22\tHAVANA\tgene\t6100000\t6300000\t.\t-\t.\tgene_id "G3"; gene_type "protein_coding"; gene_name "GENEB";\n')
r = run("report", "--sample", "S1", "--ibdmix-dir", idir, "--spans", f"{TMP}/spans.tsv", "--genes", f"{TMP}/genes.gtf",
        "--out-dir", f"{TMP}/out")
check("report exits 0", r.returncode == 0, r.stderr)
segs = [l.split("\t") for l in open(f"{TMP}/out/S1_archaic_segments.tsv").read().splitlines()]
check("segments header", segs[0] == ar.SEGMENT_HEADER, segs[0])
body = segs[1:]
check("two segments: the merged Neanderthal one and the Denisovan one", len(body) == 2, body)
check("Altai 1.0-1.3 and Vindija 1.2-1.5 Mb merge into Neanderthal 1.0-1.5 Mb",
      body and body[0][:6] == ["chr22", "1000000", "1500000", "500000", "Neanderthal", "Altai,Vindija,Denisova"], body[:1])
check("the merged segment keeps the higher LOD", body and body[0][6] == "50.0", body[:1])
check("a Denisova stretch inside a Neanderthal one is not Denisovan; one apart is", len(body) > 1 and body[1][:5] == ["chr22", "6000000", "6200000", "200000", "Denisovan"], body)
check("protein-coding genes only", body and body[0][7] == "GENEA", body[:1])
check("a stand-in's segment is not yours, nor one under 50 kb or LOD 4", all(r[1] not in ("5000000", "2000000", "3000000") for r in body), body)
summ = dict(l.split("\t", 1) for l in open(f"{TMP}/out/S1_archaic.tsv").read().splitlines()[1:])
check("Neanderthal share of the 10 Mb span: 5.00%", summ.get("neanderthal_percent") == "5.00", summ)
check("Denisovan share: 2.00%", summ.get("denisovan_percent") == "2.00", summ)
check("segment counts", summ.get("neanderthal_segments") == "1" and summ.get("denisovan_segments") == "1", summ)
check("the longest segment is listed first with its genes", summ.get("longest_1", "").startswith("chr22:1000000-1500000 Neanderthal")
      and "GENEA" in summ.get("longest_1", ""), summ.get("longest_1"))
txt = open(f"{TMP}/out/S1_archaic.txt").read()
check("the text names the caveat", "An estimate, not a measurement." in txt, txt[-300:])
r = run("report", "--sample", "S1", "--ibdmix-dir", f"{TMP}/nothing", "--spans", f"{TMP}/spans.tsv", "--out-dir", f"{TMP}/out2")
check("a missing IBDmix output stops the report", r.returncode != 0 and "no IBDmix output" in r.stderr, r.stderr)

# --- 3b. the reports read it ------------------------------------------------------
import collect_summary as cs  # noqa: E402
import render_report as rr  # noqa: E402

sdir = f"{TMP}/sample"
os.makedirs(f"{sdir}/archaic")
shutil.copy(f"{TMP}/out/S1_archaic.tsv", f"{sdir}/archaic/S1_archaic.tsv")
path, data = cs.sec_archaic(sdir, "S1")
check("collect_summary reads the summary", path and data.get("neanderthal_percent") == "5.00"
      and data.get("denisovan_segments") == 1 and len(data.get("longest", [])) == 2, data)
check("the archaic section is in SECTIONS as step 43", any(x[0] == "archaic" and x[2] == "43" for x in cs.SECTIONS))
lines = rr.archaic_lines(data)
check("the report lines give both shares", lines[0].startswith("Neanderthal: 5.00% (0.5 Mb in 1 segments)")
      and lines[1].startswith("Denisovan: 2.00%"), lines)
write(f"{sdir}/archaic/S1_archaic.tsv", "key\tvalue\nneanderthal_percent\tlots\n")
try:
    cs.sec_archaic(sdir, "S1")
    check("a damaged summary is unreadable, not a share", False)
except cs.Unreadable:
    check("a damaged summary is unreadable, not a share", True)
write(f"{sdir}/archaic/S1_archaic.tsv", "key\tvalue\nneanderthal_segments\n")
try:
    cs.sec_archaic(sdir, "S1")
    check("a row cut after its key is unreadable, not a crash", False)
except cs.Unreadable:
    check("a row cut after its key is unreadable, not a crash", True)
genes = {"chr22": [(90, 100, "ENDS_AT_START"), (100, 120, "INSIDE"), (200, 210, "STARTS_AT_END"), (50, 99, "BEFORE")]}
check("genes: a segment [100, 200) takes a gene ending at its first site, not one starting at its end",
      ar.genes_in(genes, "22", 100, 200) == ["ENDS_AT_START", "INSIDE"], ar.genes_in(genes, "22", 100, 200))
check("no file: the section is missing", cs.sec_archaic(f"{TMP}/nowhere", "S1") == (None, {}))

# --- 4. the real IBDmix -------------------------------------------------------------
NEA = (2000000, 2600000)
DEN = (3200000, 3600000)


def synthetic(d, plant):
    """A 4 Mb chr22 table and one person. A third of the sites carry an allele
    the Neanderthals share and modern people rarely do, a third one only the
    Denisovan carries, a third are common. The person draws both copies from
    the frequencies, except inside the planted stretches, where one copy
    carries the archaic allele."""
    rnd = random.Random(11)
    os.makedirs(f"{d}/table")
    rows, calls = [], []
    for i, pos in enumerate(range(1000000, 5000000, 1000)):
        kind = i % 3
        if kind == 0:
            f, al, vi, de = rnd.uniform(0.01, 0.04), 2, 2, 0
        elif kind == 1:
            f, al, vi, de = rnd.uniform(0.01, 0.04), 0, 0, 2
        else:
            f = rnd.uniform(0.1, 0.9)
            al = vi = de = sum(rnd.random() < f for _ in range(2))
        rows.append(f"chr22\t{pos}\tA\tG\t{al}\t{vi}\t{de}\t" + "\t".join([f"{f:.4f}"] * 6))
        g = [rnd.random() < f, rnd.random() < f]
        if plant and NEA[0] <= pos < NEA[1]:
            g[0] = al == 2
        if plant and DEN[0] <= pos < DEN[1]:
            g[0] = de == 2
        if any(g):
            calls.append(f"chr22\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{'1/1' if all(g) else '0/1'}")
    with gzip.open(f"{d}/table/chr22.tsv.gz", "wt") as f:
        f.write("\t".join(ar.TABLE_HEADER) + "\n" + "\n".join(rows) + "\n")
    write(f"{d}/s.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" + "\n".join(calls) + "\n")


def ibdmix_run(d):
    r = run("genotypes", "--table-dir", f"{d}/table", "--vcf", f"{d}/s.vcf", "--out-dir", f"{d}/gt")
    if r.returncode:
        return r
    os.makedirs(f"{d}/ibd")
    # The flags of scripts/43-archaic.sh.
    sh = ("set -e; gzip -dc /w/gt/chr22.gt.gz > /w/gt/chr22.gt; for a in Altai Vindija Denisova; do "
          "ibdmix -g /w/gt/chr22.gt -s /w/gt/samples.txt -n $a -o /w/ibd/chr22.$a.txt "
          "-d 4 -m 1 -a 0.01 -e 0.002 -c 2; done")
    r = subprocess.run(["docker", "run", "--rm", "--network", "none", "-u", f"{os.getuid()}:{os.getgid()}",
                        "-v", f"{d}:/w", image("IBDMIX_IMAGE"), "sh", "-c", sh], capture_output=True, text=True)
    if r.returncode:
        return r
    return run("report", "--sample", "S", "--ibdmix-dir", f"{d}/ibd", "--spans", f"{d}/gt/spans.tsv", "--out-dir", f"{d}/out")


if os.environ.get("ARCHAIC_TEST_NO_DOCKER") == "1":
    print("[SKIP] part 4: ARCHAIC_TEST_NO_DOCKER=1")
else:
    for plant in (True, False):
        d = f"{TMP}/real_{plant}"
        os.makedirs(d)
        os.chmod(d, 0o777)
        synthetic(d, plant)
        r = ibdmix_run(d)
        check(f"IBDmix chain exits 0 ({'planted' if plant else 'control'})", r.returncode == 0, r.stderr + r.stdout)
        if r.returncode:
            continue
        body = [l.split("\t") for l in open(f"{d}/out/S_archaic_segments.tsv").read().splitlines()[1:]]
        print(f"  segments ({'planted' if plant else 'control'}): {[(b[1], b[2], b[4]) for b in body]}")
        if plant:
            nea = [b for b in body if b[4] == "Neanderthal"]
            den = [b for b in body if b[4] == "Denisovan"]
            cover = lambda segs, lo, hi: sum(max(0, min(int(b[2]), hi) - max(int(b[1]), lo)) for b in segs)
            check("the planted Neanderthal stretch: one Neanderthal segment covering 90% of it",
                  len(nea) == 1 and cover(nea, *NEA) >= 0.9 * (NEA[1] - NEA[0]), nea)
            check("and no more than 100 kb beyond it", nea and int(nea[0][3]) <= NEA[1] - NEA[0] + 100000, nea)
            check("the planted Denisovan stretch: one Denisovan segment covering 90% of it",
                  len(den) == 1 and cover(den, *DEN) >= 0.9 * (DEN[1] - DEN[0]), den)
        else:
            check("without the planted stretches: no segment", body == [], body)

shutil.rmtree(TMP, ignore_errors=True)
print(f"{FAILS} check(s) failed" if FAILS else "all checks passed")
sys.exit(1 if FAILS else 0)
