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
  6. the command line writes the table, and --consensus adds the support;
  7. the rule table and the depth check, through the command line with
     assets/sv_rule_regions.tsv (each check reads the written table, so on a
     postfilter without them every check fails with what it got):
     an in-frame deletion inside DSPP's coding tandem repeat is class 3 and
     the same deletion reaching 254 bp past the repeat stays class 4; a
     deletion inside IGH is class 3, one 66% inside it too, and one 39%
     inside it stays; a
     deletion inside the MHC is class 3; a class 5 deletion on the
     pseudogene RHEBP1, which AnnotSV names RHEB (its alias), is not mapped
     to RHEB and is class 3, while one on RHEB itself stays RHEB and 5; two
     Manta deletions at one place, DHFFC 0.95 and 0.5: the first is class 3
     (no depth support), the second stays 4; a duplication with DHBFC 1.1
     is class 3; Delly RDCN 2 on a deletion is no support, RDCN 1 is; on
     chrX a duplication with RDCN 2 counts; a record with no depth value is
     "not checked" and keeps its class; a CNVpytor call has support; a
     Delly call with FORMAT GT only (as step 05 adds it) has support over
     1,000 bp and is "not checked" at 800 bp; the
     summary line counts the rows moved to class 3;
  8. red controls: without --rules the DSPP and IGH deletions stay class 4;
     and the rule table itself: every row has a known rule, IGH lies on
     chr14 between 105.5 and 107 Mb, the MHC is GRCh38's 28,510,120 to
     33,480,577, and there are at least 2,000 coding tandem repeats.

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


RULES_TSV = os.path.join(REPO, "assets", "sv_rule_regions.tsv")
HEAD7 = ["AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type", "Samples_ID", "ID", "INFO", "FORMAT",
         "S1", "Annotation_mode", "Gene_name", "B_loss_AFmax", "ACMG_class", "AnnotSV_ranking_score"]


def row7(i, chrom, start, end, t, cls, caller, gene, fmt="GT", sample="0/1"):
    length = end - start if t != "DEL" else start - end
    v = [f"{chrom}_{start}_{end}_{t}_1", chrom, str(start), str(end), str(length), t, "S1", i,
         f"END={end};SVTYPE={t};SVCALLER={caller}", fmt, sample, "full", gene, "", str(cls), "0.99"]
    return "\t".join(v) + "\n"


ROWS7 = "".join([
    row7("TR_IN", "4", 87614600, 87615500, "DEL", 4, "manta", "DSPP", "GT:DHFFC", "0/1:0.45"),
    row7("TR_OUT", "4", 87614600, 87616700, "DEL", 4, "manta", "DSPP", "GT:DHFFC", "0/1:0.45"),
    row7("IGH_IN", "14", 106000000, 106100000, "DEL", 4, "manta", "IGHV3-23", "GT:DHFFC", "0/1:0.5"),
    row7("IGH_EDGE", "14", 105400000, 105700000, "DEL", 4, "manta", "TMEM121", "GT:DHFFC", "0/1:0.5"),
    row7("IGH_MOST", "14", 105550000, 105650000, "DEL", 4, "manta", "IGHA2", "GT:DHFFC", "0/1:0.5"),
    row7("MHC_IN", "6", 31000000, 31100000, "DEL", 4, "manta", "HLA-C", "GT:DHFFC", "0/1:0.5"),
    row7("ALIAS", "10", 46634800, 46635600, "DEL", 5, "manta", "RHEB", "GT:DHFFC", "0/1:0.3"),
    row7("REAL", "7", 151470000, 151480000, "DEL", 5, "manta", "RHEB", "GT:DHFFC", "0/1:0.3"),
    row7("D_LOW", "20", 40000000, 40300000, "DEL", 4, "manta", "G7", "GT:DHFFC", "0/1:0.95"),
    row7("D_OK", "20", 40000000, 40300000, "DEL", 4, "manta", "G8", "GT:DHFFC", "0/1:0.5"),
    row7("DUP_LOW", "20", 41000000, 41300000, "DUP", 4, "manta", "G9", "GT:DHBFC", "0/1:1.1"),
    row7("DELLY2", "20", 42000000, 42300000, "DEL", 4, "delly", "G10", "GT:RDCN", "0/1:2"),
    row7("DELLY1", "20", 43000000, 43300000, "DEL", 4, "delly", "G11", "GT:RDCN", "0/1:1"),
    row7("XDUP", "X", 50000000, 50300000, "DUP", 4, "delly", "G12", "GT:RDCN", "0/1:2"),
    row7("NODEPTH", "20", 44000000, 44300000, "DEL", 4, "manta", "G13"),
    row7("CNV", "20", 45000000, 45300000, "DEL", 4, "cnvpytor", "G14"),
    row7("DELLY_GT", "20", 46000000, 46300000, "DEL", 4, "delly", "G15"),
    row7("DELLY_SHORT", "20", 47000000, 47000800, "DEL", 4, "delly", "G16"),
])
# RHEB as AnnotSV's genes file would name it: once on the pseudogene RHEBP1
# (chr10, the region of the rule table) and once on RHEB itself (chr7).
GENES7 = ("10\t46634910\t46635466\t+\tRHEB\tNM_990001\t46634910\t46635466\t46634910,\t46635466,\n"
          "7\t151466000\t151520000\t-\tRHEB\tNM_005614\t151466100\t151519900\t151466000,151519000,\t151467000,151520000,\n")


def run7(work, name, with_rules=True):
    """{ID: written row} of the command line on ROWS7, {} when it failed."""
    c45 = os.path.join(work, name + ".tsv")
    with open(c45, "w") as f:
        f.write("\t".join(HEAD7) + "\n" + ROWS7)
    genes = os.path.join(work, name + "_genes.bed")
    with open(genes, "w") as f:
        f.write(GENES7)
    dst = os.path.join(work, name + "_checked.tsv")
    cmd = [sys.executable, os.path.join(REPO, "bin", "annotsv_postfilter.py"), "--class45", c45, "--genes", genes,
           "--out", dst] + (["--rules", RULES_TSV] if with_rules else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0 or not os.path.isfile(dst):
        print(f"  (the command failed: {p.stderr.strip()[-300:]})")
        return {}, p.stderr
    with open(dst) as f:
        return {r["ID"]: r for r in csv.DictReader(f, delimiter="\t")}, p.stderr


def rules_section(work):
    got, err = run7(work, "rules")

    def col(i, c):
        return got.get(i, {}).get(c)

    # the region rules
    check("an in-frame deletion inside DSPP's coding tandem repeat is class 3, not likely pathogenic",
          (col("TR_IN", "ACMG_class"), col("TR_IN", "AnnotSV_class"), col("TR_IN", "Rule")) == ("3", "4", "tandem_repeat"),
          got.get("TR_IN"))
    check("the same deletion reaching 254 bp past the repeat stays class 4",
          (col("TR_OUT", "ACMG_class"), col("TR_OUT", "Rule")) == ("4", ""), got.get("TR_OUT"))
    check("a deletion inside IGH is class 3 (V(D)J)",
          (col("IGH_IN", "ACMG_class"), col("IGH_IN", "Rule")) == ("3", "ig_tr_locus"), got.get("IGH_IN"))
    check("a deletion 39% inside IGH stays class 4", col("IGH_EDGE", "ACMG_class") == "4", got.get("IGH_EDGE"))
    check("a deletion 66% inside IGH is class 3 (gnomAD-SV's rule: over half)", col("IGH_MOST", "ACMG_class") == "3",
          got.get("IGH_MOST"))
    check("a deletion inside the MHC is class 3",
          (col("MHC_IN", "ACMG_class"), col("MHC_IN", "Rule")) == ("3", "mhc"), got.get("MHC_IN"))
    check("a class 5 deletion on the pseudogene RHEBP1 is not mapped to RHEB",
          (col("ALIAS", "Gene_name"), col("ALIAS", "ACMG_class"), col("ALIAS", "Rule")) == ("RHEBP1", "3", "gene_alias"),
          got.get("ALIAS"))
    check("a class 5 deletion on RHEB itself stays RHEB and class 5",
          (col("REAL", "Gene_name"), col("REAL", "ACMG_class"), col("REAL", "Rule")) == ("RHEB", "5", ""), got.get("REAL"))

    # the depth check
    check("a Manta deletion with DHFFC 0.95 has no depth support and is class 3",
          (col("D_LOW", "Depth_support"), col("D_LOW", "ACMG_class"), col("D_LOW", "Rule")) == ("no", "3", "no_depth_support"),
          got.get("D_LOW"))
    check("the same deletion with DHFFC 0.5 stays class 4",
          (col("D_OK", "Depth_support"), col("D_OK", "ACMG_class")) == ("yes", "4"), got.get("D_OK"))
    check("a duplication with DHBFC 1.1 is class 3",
          (col("DUP_LOW", "Depth_support"), col("DUP_LOW", "ACMG_class")) == ("no", "3"), got.get("DUP_LOW"))
    check("Delly RDCN 2 on a deletion: no support; RDCN 1: support",
          (col("DELLY2", "Depth_support"), col("DELLY1", "Depth_support")) == ("no", "yes"),
          (got.get("DELLY2"), got.get("DELLY1")))
    check("a chrX duplication with RDCN 2 counts (the sex is not known here)", col("XDUP", "Depth_support") == "yes",
          got.get("XDUP"))
    check("no depth value: not checked, class kept",
          (col("NODEPTH", "Depth_support"), col("NODEPTH", "ACMG_class")) == ("not checked", "4"), got.get("NODEPTH"))
    check("a CNVpytor call has depth support", col("CNV", "Depth_support") == "yes", got.get("CNV"))
    check("a Delly call step 05 added (FORMAT GT only) over 1,000 bp passed its RDCN gate there: support",
          (col("DELLY_GT", "Depth_support"), col("DELLY_GT", "ACMG_class")) == ("yes", "4"), got.get("DELLY_GT"))
    check("the same call of 800 bp had no RDCN gate: not checked, class kept",
          (col("DELLY_SHORT", "Depth_support"), col("DELLY_SHORT", "ACMG_class")) == ("not checked", "4"),
          got.get("DELLY_SHORT"))
    moved = sorted(i for i, r in got.items() if r["ACMG_class"] != r["AnnotSV_class"])
    want = sorted(["TR_IN", "IGH_IN", "IGH_MOST", "MHC_IN", "ALIAS", "D_LOW", "DUP_LOW", "DELLY2"])
    check("exactly the eight rows above are moved to class 3", moved == want, moved)
    check("the summary line counts them", "moved to class 3: 8; depth support no: 3, not checked: 2" in err, err)

    # 8. red controls
    plain, _ = run7(work, "norules", with_rules=False)
    check("red control: without --rules the DSPP and IGH deletions stay class 4",
          (plain.get("TR_IN", {}).get("ACMG_class"), plain.get("IGH_IN", {}).get("ACMG_class")) == ("4", "4"),
          (plain.get("TR_IN"), plain.get("IGH_IN")))
    with open(RULES_TSV) as f:
        table = [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]
    rules = {r[3] for r in table}
    check("the rule table: every row has 7 columns and a known rule",
          all(len(r) == 7 for r in table) and rules == {"tandem_repeat", "ig_tr_locus", "mhc", "gene_alias"}, rules)
    igh = [r for r in table if r[3] == "ig_tr_locus" and r[5] == "IGH"]
    check("the rule table: IGH on chr14 between 105.5 and 107 Mb",
          len(igh) == 1 and igh[0][0] == "chr14" and 105_500_000 < int(igh[0][1]) < int(igh[0][2]) < 107_000_000, igh)
    check("the rule table: the MHC is GRCh38's 28,510,120-33,480,577",
          [r[:3] for r in table if r[3] == "mhc"] == [["chr6", "28510119", "33480577"]])
    check("the rule table: at least 2,000 coding tandem repeats", sum(r[3] == "tandem_repeat" for r in table) >= 2000)


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

        # 7 and 8. the rule table and the depth check
        rules_section(work)
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
