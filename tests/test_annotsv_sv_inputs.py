#!/usr/bin/env python3
"""bin/manta_inv_filter.py, bin/sv_svlen_fix.py and bin/cnvpytor_to_vcf.py on
made-up VCFs.

Checks:
  1. manta_inv_filter.py, on Manta's file and a convertInversion.py-style
     converted file:
       - a reciprocal pair (INV3 and INV5 junctions 100 bp apart) and a pair
         that shares one Manta EVENT keep their INV records, INV_JUNCTIONS=2;
       - a single INV3 junction 50 Mb long is no longer an INV record: its two
         breakends from Manta's file take its place, INV_JUNCTIONS=1, in
         position order; the control: the converted file holds that 50 Mb INV;
       - a single 100 kb junction stays an INV with INV_JUNCTIONS=1, a single
         400 bp one also gets INVDUP_CANDIDATE;
       - the deletion passes unchanged, and the summary line counts each kind;
  2. sv_svlen_fix.py: REF 137 bp / ALT 2 bp gets SVLEN -135 and END
     POS+134, with the old values kept (the control: the input says -136);
     a pure deletion keeps SVLEN and gets END POS+135 (POS-1+len(REF)-len(ALT)); a net
     insertion gets its net SVLEN and keeps END; symbolic, breakend and
     multi-ALT records pass unchanged; a second run changes nothing;
  3. cnvpytor_to_vcf.py: a Delly PASS deletion Manta has is not added, one
     Manta lacks is; Delly's LowQual and INV records are left out; a CNVpytor
     deletion Delly's added call covers is not added twice; CNVpytor calls
     with e-val >= 0.01 or q0 >= 0.5, and LOH, are left out; the CNVpytor
     deletion and duplication nobody else has are added, with their
     genotype; every record names its caller, the output is sorted; the
     control: Manta's file alone lacks the CNVpytor deletion.

Run: python3 tests/test_annotsv_sv_inputs.py
"""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    if ok:
        print("ok:   " + desc)
    else:
        FAILS += 1
        print("FAIL: " + desc + (": " + str(detail) if detail else ""))


def run(tool, *args, stdin=""):
    p = subprocess.run([sys.executable, os.path.join(BIN, tool)] + list(args), input=stdin,
                       capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def records(text):
    return [l.split("\t") for l in text.splitlines() if l and not l.startswith("#")]


def info(rec):
    d = {}
    for item in rec[7].split(";"):
        k, _, v = item.partition("=")
        d[k] = v if _ else True
    return d


HEAD = (
    "##fileformat=VCFv4.1\n"
    "##contig=<ID=chr1,length=248956422>\n"
    "##contig=<ID=chr20,length=64444167>\n"
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
)


def bnd(chrom, pos, rid, alt, mate, extra=""):
    return "{0}\t{1}\t{2}\tN\t{3}\t500\tPASS\tSVTYPE=BND;MATEID={4}{5}\tGT\t0/1\n".format(
        chrom, pos, rid, alt, mate, extra)


def inv(chrom, pos, rid, end, side, extra=""):
    return "{0}\t{1}\t{2}\tN\t<INV>\t500\tPASS\tEND={3};SVTYPE=INV;SVLEN={4}{5};INV{6}\tGT\t0/1\n".format(
        chrom, pos, rid, end, end - pos, extra, side)


DEL = "chr20\t40000000\tMantaDEL:9:0:0:0:0:0\tN\t<DEL>\t500\tPASS\tEND=40010000;SVTYPE=DEL;SVLEN=-10000\tGT\t0/1\n"

# Manta's own file: each junction a pair of breakends.
RAW = HEAD + "".join([
    # reciprocal inversion chr20:1,000,000-2,000,000: INV3 and INV5 junctions
    bnd("chr20", 1000000, "MantaBND:10:0:1:0:0:0:0", "N]chr20:2000000]", "MantaBND:10:0:1:0:0:0:1"),
    bnd("chr20", 1000101, "MantaBND:10:0:2:0:0:0:0", "[chr20:2000101[N", "MantaBND:10:0:2:0:0:0:1"),
    bnd("chr20", 2000000, "MantaBND:10:0:1:0:0:0:1", "N]chr20:1000000]", "MantaBND:10:0:1:0:0:0:0"),
    bnd("chr20", 2000101, "MantaBND:10:0:2:0:0:0:1", "[chr20:1000101[N", "MantaBND:10:0:2:0:0:0:0"),
    # one INV3 junction, 50 Mb
    bnd("chr20", 5000000, "MantaBND:20:0:1:0:0:0:0", "N]chr20:55000000]", "MantaBND:20:0:1:0:0:0:1"),
    bnd("chr20", 55000000, "MantaBND:20:0:1:0:0:0:1", "N]chr20:5000000]", "MantaBND:20:0:1:0:0:0:0"),
    # one INV5 junction, 400 bp
    bnd("chr20", 30000000, "MantaBND:30:0:1:0:0:0:0", "[chr20:30000400[N", "MantaBND:30:0:1:0:0:0:1"),
    bnd("chr20", 30000400, "MantaBND:30:0:1:0:0:0:1", "[chr20:30000000[N", "MantaBND:30:0:1:0:0:0:0"),
    # one INV3 junction, 100 kb
    bnd("chr1", 3000000, "MantaBND:40:0:1:0:0:0:0", "N]chr1:3100000]", "MantaBND:40:0:1:0:0:0:1"),
    bnd("chr1", 3100000, "MantaBND:40:0:1:0:0:0:1", "N]chr1:3000000]", "MantaBND:40:0:1:0:0:0:0"),
    DEL,
])
# What convertInversion.py writes for it (an INV5 record moves one base left).
EV = ";EVENT=MantaINV:50:0:1:0:0:0"
CONVERTED = HEAD + "".join([
    inv("chr1", 3000000, "MantaINV:40:0:1:0:0:0", 3100000, 3),
    # two junctions 3 kb apart that Manta put in one EVENT
    inv("chr1", 10000000, "MantaINV:50:0:1:0:0:0", 10500000, 3, EV),
    inv("chr1", 10003000, "MantaINV:50:0:2:0:0:0", 10503000, 5, EV),
    inv("chr20", 1000000, "MantaINV:10:0:1:0:0:0", 2000000, 3),
    inv("chr20", 1000100, "MantaINV:10:0:2:0:0:0", 2000100, 5),
    inv("chr20", 5000000, "MantaINV:20:0:1:0:0:0", 55000000, 3),
    inv("chr20", 29999999, "MantaINV:30:0:1:0:0:0", 30000399, 5),
    DEL,
])


def test_inv_filter(w):
    raw, conv = os.path.join(w, "raw.vcf.gz"), os.path.join(w, "conv.vcf")
    with gzip.open(raw, "wt") as f:
        f.write(RAW)
    with open(conv, "w") as f:
        f.write(CONVERTED)

    long_inv = [r for r in records(CONVERTED) if info(r).get("SVTYPE") == "INV" and int(info(r)["END"]) - int(r[1]) >= 1000000
                and r[2] == "MantaINV:20:0:1:0:0:0"]
    check("control: the converted file holds the 50 Mb single-junction INV", len(long_inv) == 1)

    rc, out, err = run("manta_inv_filter.py", raw, conv)
    check("manta_inv_filter.py exits 0", rc == 0, err)
    recs = records(out)
    by_id = {r[2]: r for r in recs}
    for rid in ("MantaINV:10:0:1:0:0:0", "MantaINV:10:0:2:0:0:0", "MantaINV:50:0:1:0:0:0", "MantaINV:50:0:2:0:0:0"):
        check("both junctions present: {0} keeps its INV, INV_JUNCTIONS=2".format(rid),
              rid in by_id and info(by_id[rid]).get("INV_JUNCTIONS") == "2", by_id.get(rid))
    check("the 50 Mb single junction is no INV record any more", "MantaINV:20:0:1:0:0:0" not in by_id)
    mb_invs = [r[2] for r in recs if info(r).get("SVTYPE") == "INV" and int(info(r)["END"]) - int(r[1]) > 1000000]
    check("no INV record longer than the largest paired inversion (1 Mb)", mb_invs == [], mb_invs)
    for rid in ("MantaBND:20:0:1:0:0:0:0", "MantaBND:20:0:1:0:0:0:1"):
        check("its breakend {0} is back, INV_JUNCTIONS=1".format(rid),
              rid in by_id and info(by_id[rid]).get("SVTYPE") == "BND" and info(by_id[rid]).get("INV_JUNCTIONS") == "1",
              by_id.get(rid))
    r = by_id.get("MantaINV:40:0:1:0:0:0")
    check("the 100 kb single junction stays an INV, INV_JUNCTIONS=1, no INVDUP_CANDIDATE",
          r is not None and info(r).get("INV_JUNCTIONS") == "1" and "INVDUP_CANDIDATE" not in info(r), r)
    r = by_id.get("MantaINV:30:0:1:0:0:0")
    check("the 400 bp single junction gets INVDUP_CANDIDATE",
          r is not None and info(r).get("INV_JUNCTIONS") == "1" and info(r).get("INVDUP_CANDIDATE") is True, r)
    check("the deletion passes unchanged", DEL.rstrip("\n").split("\t") in recs)
    order = [(r[0], int(r[1])) for r in recs]
    check("records sorted by the ##contig order, then POS",
          order == sorted(order, key=lambda x: (["chr1", "chr20"].index(x[0]), x[1])), order)
    check("the header declares INV_JUNCTIONS and INVDUP_CANDIDATE",
          "##INFO=<ID=INV_JUNCTIONS," in out and "##INFO=<ID=INVDUP_CANDIDATE," in out)
    check("the header lines come before #CHROM", out.index("##INFO=<ID=INV_JUNCTIONS,") < out.index("#CHROM"))
    check("the summary counts each kind",
          "4 INV records with both junctions, 2 single-junction INV records kept (1 shorter than 1000 bp "
          "marked INVDUP_CANDIDATE), 1 single-junction inversions of 1000000 bp or more" in err, err)

    rc, out, err = run("manta_inv_filter.py", raw, conv, "--max-single-bp", "100000000")
    check("with a 100 Mb limit the 50 Mb junction stays an INV (INV_JUNCTIONS=1)",
          any(r[2] == "MantaINV:20:0:1:0:0:0" and info(r).get("INV_JUNCTIONS") == "1" for r in records(out)), err)


SV = (
    "##fileformat=VCFv4.1\n##contig=<ID=chr20,length=64444167>\n"
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
)


def test_svlen(w):
    ref137, alt2 = "A" * 137, "AC"
    rows = [
        ("chr20", "1005400", "complex", ref137, alt2, "END=1005536;SVTYPE=DEL;SVLEN=-136"),
        ("chr20", "1100000", "pure", ref137, "A", "END=1100136;SVTYPE=DEL;SVLEN=-136"),
        ("chr20", "1200000", "ins", "ACG", "A" + "T" * 59, "END=1200002;SVTYPE=INS;SVLEN=59"),
        ("chr20", "1300000", "sym", "N", "<DEL>", "END=1310000;SVTYPE=DEL;SVLEN=-10000"),
        ("chr20", "1400000", "bnd", "N", "N]chr20:1500000]", "SVTYPE=BND"),
        ("chr20", "1600000", "multi", ref137, "A,AC", "END=1600136;SVTYPE=DEL;SVLEN=-136"),
    ]
    vcf = SV + "".join("\t".join(r[:5] + ("50", "PASS", r[5], "GT", "0/1")) + "\n" for r in rows)
    src = os.path.join(w, "sv.vcf")
    with open(src, "w") as f:
        f.write(vcf)
    by_id = {r[2]: r for r in records(vcf)}
    check("control: the input gives the REF 137 / ALT 2 record SVLEN -136", info(by_id["complex"])["SVLEN"] == "-136")

    rc, out, err = run("sv_svlen_fix.py", src)
    check("sv_svlen_fix.py exits 0", rc == 0, err)
    by_id = {r[2]: r for r in records(out)}
    c = info(by_id["complex"])
    check("REF 137 / ALT 2: SVLEN -135, END POS+134 (an interval of 135 bp for AnnotSV)",
          c.get("SVLEN") == "-135" and c.get("END") == "1005534", c)
    check("REF 137 / ALT 2: the old values kept", c.get("ORIG_SVLEN") == "-136" and c.get("ORIG_END") == "1005536", c)
    p = info(by_id["pure"])
    check("pure deletion: SVLEN stays -136, END becomes POS+135", p.get("SVLEN") == "-136" and p.get("END") == "1100135"
          and "ORIG_SVLEN" not in p, p)
    i = info(by_id["ins"])
    check("net insertion: SVLEN 57, END kept", i.get("SVLEN") == "57" and i.get("END") == "1200002", i)
    for rid in ("sym", "bnd", "multi"):
        check("{0} passes unchanged".format(rid), by_id[rid] == [x for x in records(vcf) if x[2] == rid][0])
    check("the header declares ORIG_SVLEN and ORIG_END", "##INFO=<ID=ORIG_SVLEN," in out and "##INFO=<ID=ORIG_END," in out)
    check("the summary counts the rewritten records", "SV length check: 3 sequence-resolved records" in err, err)
    again = os.path.join(w, "sv2.vcf")
    with open(again, "w") as f:
        f.write(out)
    rc, out2, err2 = run("sv_svlen_fix.py", again)
    check("a second run changes nothing", rc == 0 and out2 == out, err2)


BASE = (
    "##fileformat=VCFv4.1\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chr20,length=64444167>\n"
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    "chr20\t1000000\tMantaDEL:1\tN\t<DEL>\t500\tPASS\tEND=1100000;SVTYPE=DEL;SVLEN=-100000\tGT\t0/1\n"
    "chr20\t3000000\tMantaINV:2\tN\t<INV>\t500\tPASS\tEND=4000000;SVTYPE=INV;SVLEN=1000000;INV3\tGT\t0/1\n"
)
DELLY = (
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
    "chr20\t1000500\tDEL00000001\tN\t<DEL>\t60\tPASS\tIMPRECISE;SVTYPE=DEL;END=1099000;CIPOS=-50,50;CIEND=-50,50\tGT:GQ\t0/1:40\n"
    "chr20\t2000000\tDEL00000002\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;END=2050000\tGT:GQ\t1/1:40\n"
    "chr20\t2500000\tDEL00000003\tN\t<DEL>\t10\tLowQual\tIMPRECISE;SVTYPE=DEL;END=2600000\tGT\t0/1\n"
    "chr20\t3000000\tINV00000004\tN\t<INV>\t60\tPASS\tIMPRECISE;SVTYPE=INV;END=4000000;CT=3to3\tGT\t0/1\n"
)


def cnv(pos, rid, alt, end, p1, q0, gt):
    t = alt[1:4]
    return ("chr{0}\t{1}\t{2}\tN\t{3}\t.\tPASS\tEND={4};IMPRECISE;SVLEN={5};SVTYPE={6};pytorRD=0.5;pytorP1={7};"
            "pytorP2=0;pytorP3=0;pytorP4=0;pytorQ0={8};pytorPN=0;pytorDG=1000;pytorCL=rd_mean_shift\tGT:CN\t{9}:1\n"
            ).format(pos[0], pos[1], rid, alt, end, end - pos[1], t, p1, q0, gt)


CNVPYTOR = "##fileformat=VCFv4.2\n##source=CNVpytor\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" + "".join([
    cnv(("20", 2001000), "CNVpytor_del0", "<DEL>", 2049000, "1e-06", "0.01", "1/1"),   # Delly's added call has it
    cnv(("20", 3200000), "CNVpytor_del1", "<DEL>", 3700000, "1e-08", "0.02", "0/1"),   # inside Manta's INV: added
    cnv(("20", 5000000), "CNVpytor_del2", "<DEL>", 5100000, "0.5", "0.0", "0/1"),      # e-val too high
    cnv(("20", 6000000), "CNVpytor_del3", "<DEL>", 6100000, "1e-08", "0.9", "0/1"),    # q0 too high
    cnv(("20", 7000000), "CNVpytor_loh4", "<LOH>", 7100000, "1e-08", "0.0", "0/1"),    # LOH
    cnv(("1", 500000), "CNVpytor_dup5", "<DUP>", 600000, "1e-05", "0.0", "0/1"),       # added
])


def test_extra_callers(w):
    base, delly, cnvp = os.path.join(w, "base.vcf"), os.path.join(w, "delly.vcf.gz"), os.path.join(w, "cnv.vcf.gz")
    with open(base, "w") as f:
        f.write(BASE)
    with gzip.open(delly, "wt") as f:
        f.write(DELLY)
    with gzip.open(cnvp, "wt") as f:
        f.write(CNVPYTOR)
    check("control: Manta's file alone lacks the CNVpytor deletion", "CNVpytor_del1" not in BASE)

    rc, out, err = run("cnvpytor_to_vcf.py", "--base", base, "--delly", delly, "--cnvpytor", cnvp)
    check("cnvpytor_to_vcf.py exits 0", rc == 0, err)
    recs = records(out)
    ids = [r[2] for r in recs]
    check("the added records are exactly Delly's DEL00000002 and CNVpytor's del1 and dup5",
          ids == ["CNVpytor_dup5", "MantaDEL:1", "DEL00000002", "MantaINV:2", "CNVpytor_del1"], ids)
    by_id = {r[2]: r for r in recs}
    check("Manta's records name their caller", all(info(by_id[i]).get("SVCALLER") == "manta" for i in ("MantaDEL:1", "MantaINV:2")))
    d = by_id.get("DEL00000002")
    check("Delly's call: SVCALLER=delly, <DEL>, END kept, its genotype",
          d is not None and info(d).get("SVCALLER") == "delly" and d[4] == "<DEL>" and info(d).get("END") == "2050000"
          and d[8:] == ["GT", "1/1"], d)
    c = by_id.get("CNVpytor_del1")
    check("CNVpytor's deletion: SVCALLER=cnvpytor, END, no SVLEN (its positive one), IMPRECISE, its genotype",
          c is not None and info(c).get("SVCALLER") == "cnvpytor" and info(c).get("END") == "3700000"
          and "SVLEN" not in info(c) and info(c).get("IMPRECISE") is True and c[8:] == ["GT", "0/1"], c)
    c = by_id.get("CNVpytor_dup5")
    check("CNVpytor's duplication is a <DUP>", c is not None and c[4] == "<DUP>" and info(c).get("SVTYPE") == "DUP", c)
    check("the header declares SVCALLER", "##INFO=<ID=SVCALLER," in out)
    check("the summary counts what each caller gave",
          "delly: 2 PASS DEL/DUP calls, 1 added" in err and "cnvpytor: 3 DEL/DUP calls with e-val < 0.01 and q0 < 0.5, 2 added" in err,
          err)

    rc, out, err = run("cnvpytor_to_vcf.py", "--base", "-", "--base-name", "sv_vcf", stdin=BASE)
    check("no other caller: stdin in, the base out with SVCALLER=sv_vcf",
          rc == 0 and "none given" in err and [info(r).get("SVCALLER") for r in records(out)] == ["sv_vcf"] * 2, err)


def main():
    w = tempfile.mkdtemp()
    try:
        test_inv_filter(w)
        test_svlen(w)
        test_extra_callers(w)
    finally:
        shutil.rmtree(w)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else "{0} FAILED".format(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
