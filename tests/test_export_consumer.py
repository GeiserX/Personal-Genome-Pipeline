#!/usr/bin/env python3
"""bin/export_consumer.py (step 40) on synthetic files.

Checks:
  1. build-sites: an rsID gets its GRCh38 position from dbSNP; a GRCh37 base
     on the other strand gives strand '-'; 23andMe's own i-numbers, rsIDs
     dbSNP lacks, non-SNV rsIDs, rsIDs dbSNP places twice or on another
     chromosome, and a base on neither strand are dropped, each with its
     reason; a pseudoautosomal rsID with a chrX and a chrY record keeps the
     chrX one; an MT site takes dbSNP's rCRS position; the table is sorted
     and its header counts the sites;
  2. export: a site inside a reference block is homozygous reference, a het
     and a hom-alt call keep their bases, a filtered call, a no-call block,
     a gap between records, a site inside a called deletion (also under
     the block DeepVariant starts inside it) and a contig the gVCF lacks are
     '--'; a block over a reference-called deletion is reference; a variant
     record wins over a block at the same position; strand '-' is complemented; a male's chrX outside the PAR
     gets one letter ('--' when het), inside the PAR two; chrY and chrM one;
     a female's chrY is '--'; the summary counts each reason;
  3. both layouts pass the strict reader, the AncestryDNA file has codes
     23, 24, 25 (X PAR on GRCh37) and 26 and '0 0' for no-calls, and reading
     the 23andMe file back gives each site the genotype the gVCF holds;
  4. the strict reader fails on CRLF line ends, a missing column header, a
     bad genotype, one letter on an autosome, rows out of order, a repeated
     rsID and a cut last line;
  5. a planted hole (the gVCF record over a site removed) makes that site
     '--' and leaves the others as they were;
  6. the site table shifted by one position makes the export stop with the
     reference-base message and write nothing;
  7. an unsorted gVCF stops the export.

Run: python3 tests/test_export_consumer.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(REPO, "bin", "export_consumer.py")
PAR = os.path.join(REPO, "assets", "par_grch38.bed")
sys.path.insert(0, os.path.join(REPO, "bin"))
import export_consumer  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def run(*args):
    p = subprocess.run([sys.executable, TOOL, *args], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def write(path, text, gz=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with (gzip.open(path, "wt") if gz else open(path, "w")) as f:
        f.write(text)


def table(path, rows, template_sites=None):
    head = (f"## template_sites: {template_sites or len(rows)}\n## dropped: 0\n"
            "#rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\n")
    write(path, head + "".join("\t".join(str(x) for x in r) + "\n" for r in rows), gz=True)


def unzip(path):
    with zipfile.ZipFile(path) as z:
        return z.read(z.namelist()[0]).decode()


def calls_23andme(text):
    """The genotype of each row, keyed by the number of its rsID."""
    return {int(c[0][2:]): c[3] for c in (l.split("\t") for l in text.splitlines() if l and not l.startswith("#"))}


# --- 1. build-sites ----------------------------------------------------------
TEMPLATE = "".join(f"{c}\t{p}\t{i}\t{r}\n" for c, p, i, r in [
    ("chr1", 1000, "rs10", "A"),      # kept, + strand
    ("chr1", 2000, "rs11", "T"),      # kept, - strand (dbSNP REF A, ALT G)
    ("chr1", 3000, "i7001", "C"),     # no rsID
    ("chr1", 4000, "rs12", "G"),      # not in dbSNP
    ("chr1", 5000, "rs13", "C"),      # only an indel in dbSNP
    ("chr1", 6000, "rs14", "C"),      # two GRCh38 positions
    ("chr2", 7000, "rs15", "C"),      # on chr3 in GRCh38
    ("chr1", 8000, "rs16", "C"),      # neither strand: dbSNP A/T
    ("chrX", 9000, "rs17", "G"),      # PAR: chrX and chrY records
    ("chrM", 3244, "rs18", "G"),      # hg19 chrM numbering; rCRS 3243 A
    ("chr1", 500, "rs19", "G"),       # kept, sorts first
])
DBSNP = "".join("\t".join(r) + "\n" for r in [
    ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
    ("NC_000001.11", "900", "rs10", "A", "G", ".", ".", "RS=10;VC=SNV"),
    ("NC_000001.11", "1900", "rs11", "A", "G", ".", ".", "RS=11;VC=SNV"),
    ("NC_000001.11", "4900", "rs13", "CTG", "C", ".", ".", "RS=13;VC=INDEL"),
    ("NC_000001.11", "5900", "rs14", "C", "T", ".", ".", "RS=14;VC=SNV"),
    ("NC_000001.11", "5950", "rs14", "C", "T", ".", ".", "RS=14;VC=SNV"),
    ("NC_000003.12", "6900", "rs15", "C", "T", ".", ".", "RS=15;VC=SNV"),
    ("NC_000001.11", "7900", "rs16", "A", "T", ".", ".", "RS=16;VC=SNV"),
    ("NC_000023.11", "8900", "rs17", "G", "A", ".", ".", "RS=17;VC=SNV"),
    ("NC_000024.10", "8900", "rs17", "G", "A", ".", ".", "RS=17;VC=SNV"),
    ("NC_012920.1", "3243", "rs18", "A", "G", ".", ".", "RS=18;VC=SNV"),
    ("NT_187361.1", "10", "rs10", "A", "G", ".", ".", "RS=10;VC=SNV"),
    ("NC_000001.11", "400", "rs19", "G", "C", ".", ".", "RS=19;VC=SNV"),
])


def test_build(work):
    tpl, db = os.path.join(work, "tpl.txt.gz"), os.path.join(work, "dbsnp.vcf")
    write(tpl, TEMPLATE, gz=True)
    write(db, DBSNP)
    out, dropped = os.path.join(work, "sites.tsv.gz"), os.path.join(work, "dropped.tsv")
    rc, log = run("build-sites", "--template", tpl, "--dbsnp", db, "--out", out, "--dropped", dropped, "--source", "unit")
    check("build-sites exits 0", rc == 0, log)
    meta, sites = export_consumer.read_sites(out)
    got = {s["rsid"]: (s["chrom"], s["pos37"], s["pos38"], s["ref38"], s["strand"]) for s in sites}
    check("kept: rs10 (+), rs11 (-), rs17 (chrX record), rs18 (rCRS 3243 in both builds), rs19",
          got == {"rs10": ("1", 1000, 900, "A", "+"), "rs11": ("1", 2000, 1900, "A", "-"),
                  "rs17": ("X", 9000, 8900, "G", "+"), "rs18": ("MT", 3243, 3243, "A", "+"),
                  "rs19": ("1", 500, 400, "G", "+")}, got)
    check("sorted by chromosome, then GRCh37 position", [s["rsid"] for s in sites] == ["rs19", "rs10", "rs11", "rs17", "rs18"],
          [s["rsid"] for s in sites])
    check("header counts the template and the dropped sites",
          meta.get("template_sites") == "11" and meta.get("kept") == "5" and meta.get("dropped") == "6", meta)
    reasons = {l.split("\t")[0]: l.rstrip("\n").split("\t")[3] for l in open(dropped) if not l.startswith("rsid")}
    check("each dropped site with its reason", reasons == {
        "i7001": "no rsID (a 23andMe internal id)", "rs12": "not in dbSNP (merged or withdrawn)",
        "rs13": "not an SNV in dbSNP", "rs14": "more than one GRCh38 position in dbSNP",
        "rs15": "on chromosome 3 in GRCh38", "rs16": "its GRCh37 base is not a dbSNP allele on either strand"}, reasons)
    p = subprocess.run([sys.executable, TOOL, "build-sites", "--template", tpl, "--dbsnp", "-", "--out", out + "2",
                        "--dropped", dropped + "2"], input=DBSNP, capture_output=True, text=True)
    check("build-sites on stdin gives the same table",
          p.returncode == 0 and [s["rsid"] for s in export_consumer.read_sites(out + "2")[1]] == [s["rsid"] for s in sites],
          p.stdout + p.stderr)


# --- 2-5. export ---------------------------------------------------------------
H = ("##fileformat=VCFv4.2\n##contig=<ID=chr20>\n##contig=<ID=chrX>\n##contig=<ID=chrY>\n##contig=<ID=chrM>\n"
     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")


def rec(c, p, ref, alt, gt, flt=".", end=None):
    info = f"END={end}" if end else "."
    return f"{c}\t{p}\t.\t{ref}\t{alt}\t30\t{flt}\t{info}\tGT:GQ\t{gt}:30\n"


GVCF_LINES = [
    rec("chr20", 90, "A", "<*>", "0/0", end=120),          # s1 at 100
    rec("chr20", 150, "C", "T,<*>", "0/1", "PASS"),        # s2
    rec("chr20", 151, "G", "<*>", "0/0", end=159),
    rec("chr20", 160, "G", "A,<*>", "1/1", "PASS"),        # s3
    rec("chr20", 170, "T", "C,<*>", "0/1", "LowQual"),     # s4 filtered
    rec("chr20", 175, "A", "<*>", "./.", end=185),         # s5 no-call block
    rec("chr20", 186, "A", "<*>", "0/0", end=199),         # gap 200-204: s6
    rec("chr20", 205, "A", "<*>", "0/0", end=209),
    rec("chr20", 210, "GTGA", "G,<*>", "0/1", "PASS"),     # s7 at 212 inside the deletion
    rec("chr20", 211, "T", "<*>", "0/0", end=213),         # DeepVariant's block inside the deleted bases
    rec("chr20", 214, "C", "<*>", "0/0", end=229),         # s8 at 220, strand -
    rec("chr20", 230, "G", "C,<*>", "0/1", "PASS"),        # s9: variant before the block at the same position
    rec("chr20", 230, "G", "<*>", "0/0", end=240),
    rec("chr20", 250, "ACG", "A,<*>", "0/0", "RefCall"),   # a reference call of a deletion...
    rec("chr20", 251, "C", "<*>", "0/0", end=260),         # ...and the block over s10 at 252: reference
    rec("chrX", 20000, "A", "G,<*>", "0/1", "PASS"),       # sx1 PAR1
    rec("chrX", 5000000, "C", "T,<*>", "1/1", "PASS"),     # sx2
    rec("chrX", 5000100, "C", "T,<*>", "0/1", "PASS"),     # sx3
    rec("chrX", 5000150, "G", "<*>", "0", end=5000300),    # sx4 haploid block
    rec("chrY", 3000000, "A", "C,<*>", "1", "PASS"),       # sy1
    rec("chrM", 50, "A", "<*>", "0/0", end=150),           # sm1 at 100
    rec("chrM", 200, "G", "A,<*>", "0/1", "PASS"),         # sm2 heteroplasmic
]
SITES = [
    ("rs101", "20", 1100, "A", 100, "A", "+"), ("rs102", "20", 1150, "C", 150, "C", "+"),
    ("rs103", "20", 1160, "G", 160, "G", "+"), ("rs104", "20", 1170, "T", 170, "T", "+"),
    ("rs105", "20", 1180, "A", 180, "A", "+"), ("rs106", "20", 1200, "C", 200, "C", "+"),
    ("rs107", "20", 1212, "G", 212, "G", "+"), ("rs108", "20", 1220, "A", 220, "T", "-"),
    ("rs109", "20", 1230, "G", 230, "G", "+"), ("rs110", "1", 1000, "A", 1000, "A", "+"),
    ("rs111", "20", 1252, "C", 252, "C", "+"),
    ("rs201", "X", 70000, "A", 20000, "A", "+"), ("rs202", "X", 4900000, "C", 5000000, "C", "+"),
    ("rs203", "X", 4900100, "C", 5000100, "C", "+"), ("rs204", "X", 4900200, "G", 5000200, "G", "+"),
    ("rs301", "Y", 2900000, "A", 3000000, "A", "+"), ("rs401", "MT", 100, "A", 100, "A", "+"),
    ("rs402", "MT", 200, "G", 200, "G", "+"),
]
# Keyed by the number after rs (the public-data scan reads an rsID followed by a genotype as a leak).
WANT_MALE = {101: "AA", 102: "CT", 103: "AA", 104: "--", 105: "--", 106: "--", 107: "--", 108: "AA", 109: "CG",
             110: "--", 111: "CC", 201: "AG", 202: "T", 203: "--", 204: "G", 301: "C", 401: "A", 402: "--"}


def export(work, name, gvcf_lines, sites, sex):
    d = os.path.join(work, name)
    gv, st = os.path.join(d, "in.g.vcf.gz"), os.path.join(d, "sites.tsv.gz")
    write(gv, H + "".join(gvcf_lines), gz=True)
    table(st, sites)
    rc, log = run("export", "--gvcf", gv, "--sites", st, "--sample", "S", "--sex", sex, "--par", PAR,
                  "--out-dir", os.path.join(d, "out"))
    return rc, log, os.path.join(d, "out")


def summary(out):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(out, "S_consumer_export.tsv")) if l.strip())


def test_export(work):
    rc, log, out = export(work, "male", GVCF_LINES, SITES, "male")
    check("export (male) exits 0", rc == 0, log)
    t23 = unzip(os.path.join(out, "S_23andme_v5.zip"))
    got = calls_23andme(t23)
    check("every site has the genotype the gVCF holds (male)", got == WANT_MALE,
          {k: (got.get(k), v) for k, v in WANT_MALE.items() if got.get(k) != v})
    s = summary(out)
    check("summary counts", s.get("sites") == "18" and s.get("called") == "11" and s.get("no_call") == "7"
          and s.get("no_call: filtered call") == "1" and s.get("no_call: no gVCF record") == "2"
          and s.get("no_call: no-call in the gVCF") == "1" and s.get("no_call: inside a deletion or a longer allele") == "1"
          and s.get("reference_mismatch") == "0", s)
    check("the zips hold only the .txt files (no .txt left beside them)",
          sorted(os.listdir(out)) == ["S_23andme_v5.zip", "S_ancestrydna_layout.zip", "S_consumer_export.tsv"],
          os.listdir(out))
    for layout, f in (("23andme", "S_23andme_v5.zip"), ("ancestrydna", "S_ancestrydna_layout.zip")):
        rc2, log2 = run("check", "--layout", layout, os.path.join(out, f))
        check(f"the strict reader passes the {layout} file", rc2 == 0, log2)
    anc = unzip(os.path.join(out, "S_ancestrydna_layout.zip"))
    rows = {c[0]: c[1:] for c in (l.split("\t") for l in anc.splitlines()) if c[0].startswith("rs")}
    check("AncestryDNA: codes 23, 24, 25 (X PAR on GRCh37), 26 and 0 0 for a no-call",
          rows["rs202"] == ["23", "4900000", "T", "T"] and rows["rs301"] == ["24", "2900000", "C", "C"]
          and rows["rs201"] == ["25", "70000", "A", "G"] and rows["rs401"] == ["26", "100", "A", "A"]
          and rows["rs104"] == ["20", "1170", "0", "0"], rows)
    order = [c[1] for c in (l.split("\t") for l in anc.splitlines()) if c[0].startswith("rs") and c[0] != "rsid"]
    check("AncestryDNA rows ordered by code (25 after 24)", order == sorted(order, key=int), order)
    check("23andMe header names the layout and the build", t23.startswith("# This data file was generated by")
          and "23andMe" in t23.splitlines()[0] and "build 37" in t23 and "# rsid\tchromosome\tposition\tgenotype\n" in t23)
    check("no AncestryDNA word in the 23andMe file (chip-to-vcf.sh would read it as AncestryDNA)",
          "ancestrydna" not in t23.lower())

    rc, log, out = export(work, "female", GVCF_LINES, SITES, "female")
    got = calls_23andme(unzip(os.path.join(out, "S_23andme_v5.zip")))
    check("female: chrY '--', chrX keeps the call's ploidy", rc == 0 and got[301] == "--" and got[202] == "TT"
          and got[203] == "CT" and got[204] == "G", (rc, log, {k: got.get(k) for k in (301, 202, 203, 204)}))

    # 5. a planted hole: drop the block over rs101
    holed = [l for l in GVCF_LINES if not l.startswith("chr20\t90\t")]
    rc, log, out = export(work, "hole", holed, SITES, "male")
    got = calls_23andme(unzip(os.path.join(out, "S_23andme_v5.zip")))
    check("planted hole: the site is '--'", rc == 0 and got.get(101) == "--", (rc, log, got.get(101)))
    check("planted hole: every other site unchanged", {k: v for k, v in got.items() if k != 101}
          == {k: v for k, v in WANT_MALE.items() if k != 101})

    # 6. the table shifted by one position (GRCh38) on a run of het SNVs
    seq = "ACGT" * 300
    lines, sites = [], []
    for k, p in enumerate(range(100, 600, 10)):
        alt = "T" if seq[p - 1] != "T" else "A"
        lines.append(rec("chr20", p, seq[p - 1], f"{alt},<*>", "0/1", "PASS"))
        lines.append(rec("chr20", p + 1, seq[p], "<*>", "0/0", end=p + 9))
        sites.append((f"rs{500 + k}", "20", p + 1000, seq[p - 1], p, seq[p - 1], "+"))
    rc, log, out = export(work, "run", lines, sites, "")
    check("het run: every site called", rc == 0 and summary(out).get("called") == "50", log)
    shifted = [(r[0], r[1], r[2], r[3], r[4] + 1, r[5], r[6]) for r in sites]
    rc, log, out = export(work, "shifted", lines, shifted, "")
    check("table shifted by one position: the export stops with the reference-base message",
          rc == 1 and "another reference base than the site table" in log, (rc, log))
    check("table shifted by one position: nothing written", not os.path.exists(out) or not os.listdir(out),
          os.listdir(out) if os.path.exists(out) else None)

    # 7. an unsorted gVCF
    rc, log, _ = export(work, "unsorted", [GVCF_LINES[3], GVCF_LINES[1]], SITES, "")
    check("an unsorted gVCF stops the export", rc == 1 and "not sorted" in log, (rc, log))


# --- 4. the strict reader, red controls -----------------------------------------
GOOD = ("# This data file was generated by a test\n#\n# rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t100\tAG\nrs2\t1\t200\t--\nrs3\tX\t50\tA\nrs4\tMT\t10\tG\n")


def test_reader(work):
    p = os.path.join(work, "r.txt")
    write(p, GOOD)
    check("the strict reader passes a good file", export_consumer.check_23andme(p) == [], export_consumer.check_23andme(p))
    bad = {
        "CRLF line ends": GOOD.replace("\n", "\r\n"),
        "no column header": GOOD.replace("# rsid\tchromosome\tposition\tgenotype\n", ""),
        "a bad genotype": GOOD.replace("AG", "AN"),
        "one letter on an autosome": GOOD.replace("rs2\t1\t200\t--", "rs2\t1\t200\tA"),
        "rows out of order": GOOD.replace("rs1\t1\t100", "rs1\t1\t300"),
        "a repeated rsID": GOOD.replace("rs2\t", "rs1\t"),
        "a cut last line": GOOD[:-1],
        "two letters on MT": GOOD.replace("rs4\tMT\t10\tG", "rs4\tMT\t10\tGG"),
    }
    for what, text in bad.items():
        write(p, text)
        check(f"the strict reader fails on {what}", export_consumer.check_23andme(p) != [])


def main():
    work = tempfile.mkdtemp()
    try:
        test_build(work)
        test_export(work)
        test_reader(work)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
