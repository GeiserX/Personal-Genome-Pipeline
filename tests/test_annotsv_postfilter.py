#!/usr/bin/env python3
"""bin/annotsv_postfilter.py on a made-up class 4 and 5 table.

Checks:
  1. exon edges, on a made-up AnnotSV genes file: a deletion from the last
     base of exon 1 to the base before exon 3 of one transcript is flagged
     and named; one that starts 10 bp inside the intron is not; on a minus
     strand transcript the exons are numbered from its own start; a
     duplication is not a deletion ("."); without the genes file the column
     says "not checked";
  2. the population frequency: B_loss_AFmax for a deletion, B_gain_AFmax
     for a duplication, AF 0.05 is common, an empty one is not;
  3. the caller from INFO SVCALLER;
  4. the support from step 22's consensus: a deletion in it has support 2
     and its two callers by the VCF's sample names; a deletion only Manta
     made has support 1; a duplication of the consensus does not count for
     a deletion at the same place; without a consensus file the support is
     unknown (None), not 1;
  5. red control: with an empty consensus (header only) the two-caller
     deletion has support 1, so the support-2 check above can fail;
  6. the command line writes the table, and --consensus adds the support.

Run: python3 tests/test_annotsv_postfilter.py
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import annotsv_postfilter as ap  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(("ok:   " if ok else "FAIL: ") + desc + ("" if ok else f": {detail}"))
    if not ok:
        FAILS += 1


# chrom txStart txEnd strand gene tx cdsStart cdsEnd exonStarts exonEnds
GENES = ("20\t1000000\t1010000\t+\tFIXG1\tNM_900001\t1000100\t1009900\t1000000,1006000,1009800,\t1000200,1006100,1010000,\n"
         "20\t3000000\t3005100\t-\tMING1\tNM_900005\t3000000\t3005100\t3000000,3005000,\t3000100,3005100,\n")

HEAD = ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "ID", "INFO", "Annotation_mode",
        "Gene_name", "B_gain_source", "B_gain_AFmax", "B_loss_source", "B_loss_AFmax", "ACMG_class",
        "AnnotSV_ranking_score"]


def row(i, start, end, t, cls, caller, gain_af="", loss_af="", loss_src=""):
    v = [f"20_{start}_{end}_{t}_1", "20", str(start), str(end), str(end - start if t != "DEL" else start - end), t,
         f"SV{i}", f"END={end};SVTYPE={t};SVCALLER={caller}", "full", f"GENE{i}", "", gain_af, loss_src, loss_af,
         str(cls), "0.95"]
    return "\t".join(v) + "\n"


CLASS45 = "\t".join(HEAD) + "\n" + "".join([
    row(1, 1000200, 1009800, "DEL", 4, "manta", loss_af="0.05", loss_src="gnomAD-SV_v4.1_DEL_fixture"),
    row(2, 10100999, 10102999, "DEL", 5, "manta"),
    row(3, 20000000, 20450000, "DEL", 4, "manta"),
    row(4, 30000000, 30500000, "DUP", 4, "delly", gain_af=""),
    row(5, 1000190, 1009800, "DEL", 4, "cnvpytor"),
    row(6, 3000100, 3005000, "DEL", 4, "manta"),
])

CONS_HEAD = ("##fileformat=VCFv4.1\n##source=SURVIVOR\n"
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmanta\tdelly\n")
CONSENSUS = CONS_HEAD + (
    "chr20\t10101000\tMantaDEL:1\tN\t<DEL>\t.\tPASS\tSUPP=2;SUPP_VEC=11;SVLEN=-2000;SVTYPE=DEL;SVMETHOD=SURVIVOR1.0.7;"
    "CHR2=chr20;END=10103000;STRANDS=+-\tGT\t0/1\t0/1\n"
    "chr20\t20000010\tMantaDUP:2\tN\t<DUP>\t.\tPASS\tSUPP=2;SUPP_VEC=11;SVLEN=450000;SVTYPE=DUP;END=20450010\tGT\t0/1\t0/1\n")


def main():
    work = tempfile.mkdtemp()
    try:
        genes = os.path.join(work, "genes.RefSeq.sorted.bed")
        with open(genes, "w") as f:
            f.write(GENES)
        c45 = os.path.join(work, "S_sv_class4_5.tsv")
        with open(c45, "w") as f:
            f.write(CLASS45)
        rows = ap.read_tsv(c45)
        tx = ap.read_transcripts(genes, {"20"})
        out = {r["Gene_name"]: r for r in ap.check_rows(rows, tx)}

        # 1. exon edges
        check("an intron-exact deletion (exon 1 end to exon 3 start) is on exon edges",
              out["GENE1"]["Exon_edge"] == "yes", out["GENE1"])
        check("it names the transcript and the exons",
              out["GENE1"]["Exon_edge_detail"] == "NM_900001 (FIXG1): exon 1 and exon 3 edges", out["GENE1"])
        check("a deletion starting 10 bp inside the intron is not", out["GENE5"]["Exon_edge"] == "no", out["GENE5"])
        check("a deletion in no transcript is not", out["GENE3"]["Exon_edge"] == "no", out["GENE3"])
        check("minus strand: exons numbered from the transcript's start",
              out["GENE6"]["Exon_edge_detail"] == "NM_900005 (MING1): exon 1 and exon 2 edges", out["GENE6"])
        check("a duplication is not checked as a deletion", out["GENE4"]["Exon_edge"] == ".", out["GENE4"])
        no_genes = {r["Gene_name"]: r for r in ap.check_rows(rows, None)}
        check("without the genes file: not checked", no_genes["GENE1"]["Exon_edge"] == "not checked", no_genes["GENE1"])

        # 2. population frequency
        check("the deletion's B_loss_AFmax and its source",
              (out["GENE1"]["Population_AF"], out["GENE1"]["Population_source"]) == ("0.05", "gnomAD-SV_v4.1_DEL_fixture"),
              out["GENE1"])
        check("AF 0.05 is common", out["GENE1"]["Common"] == "yes")
        check("no benign record: no AF, not common",
              (out["GENE2"]["Population_AF"], out["GENE2"]["Common"]) == ("", "no"), out["GENE2"])
        check("a duplication reads B_gain_AFmax (empty here)", out["GENE4"]["Population_AF"] == "", out["GENE4"])

        # 3. caller
        check("the caller from INFO SVCALLER", (out["GENE2"]["Caller"], out["GENE4"]["Caller"]) == ("manta", "delly"))

        # 4. consensus support
        cons = os.path.join(work, "S_sv_consensus.vcf")
        with open(cons, "w") as f:
            f.write(CONSENSUS)
        checked = ap.check_rows(rows, tx)
        sup = {r["Gene_name"]: s for r, s in zip(checked, ap.consensus_support(checked, cons))}
        check("a deletion in the consensus: support 2, manta and delly",
              sup["GENE2"] == {"Support": 2, "Support_callers": "manta,delly"}, sup["GENE2"])
        check("a deletion only Manta made: support 1", sup["GENE3"] == {"Support": 1, "Support_callers": "manta"},
              sup["GENE3"])
        check("without a consensus file the support is unknown, not 1",
              all(s["Support"] is None for s in ap.consensus_support(checked, os.path.join(work, "none.vcf"))))

        # 5. red control: an empty consensus
        empty = os.path.join(work, "empty.vcf")
        with open(empty, "w") as f:
            f.write(CONS_HEAD)
        sup_e = ap.consensus_support(checked, empty)
        check("red control: with an empty consensus the two-caller deletion has support 1",
              sup_e[1]["Support"] == 1, sup_e[1])
        check("red control: and no row has support 2", all(s["Support"] == 1 for s in sup_e), sup_e)

        # 6. the command line
        dst = os.path.join(work, "checked.tsv")
        p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "annotsv_postfilter.py"), "--class45", c45,
                            "--genes", genes, "--consensus", cons, "--out", dst], capture_output=True, text=True)
        check("the command exits 0", p.returncode == 0, p.stderr)
        with open(dst) as f:
            got = list(csv.DictReader(f, delimiter="\t"))
        check("one row per SV", len(got) == 6, len(got))
        check("the columns", list(got[0]) == ap.KEEP + ap.ADDED + ap.SUPPORT, list(got[0]) if got else None)
        check("the summary line", "deletions on exon edges: 2; common (benign-set AF >= 1%): 1" in p.stderr, p.stderr)
        g2 = next(r for r in got if r["Gene_name"] == "GENE2")
        check("the written support", (g2["Support"], g2["Support_callers"]) == ("2", "manta,delly"), g2)
        p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "annotsv_postfilter.py"), "--class45", c45,
                            "--out", dst], capture_output=True, text=True)
        with open(dst) as f:
            head = f.readline().rstrip("\n").split("\t")
        check("without --consensus no support columns, without --genes not checked",
              p.returncode == 0 and "Support" not in head and "deletions on exon edges: not checked" in p.stderr,
              p.stderr)
        open(c45, "w").close()
        p = subprocess.run([sys.executable, os.path.join(REPO, "bin", "annotsv_postfilter.py"), "--class45", c45,
                            "--genes", genes, "--out", dst], capture_output=True, text=True)
        with open(dst) as f:
            lines = f.read().splitlines()
        check("an empty class 4 and 5 table gives a header and no row",
              p.returncode == 0 and len(lines) == 1, (p.stderr, lines))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
