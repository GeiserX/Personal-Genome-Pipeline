#!/usr/bin/env python3
"""bin/paraphase_report.py: step 35b's read-group check and its two tables.

Checks:
  1. check-bam: a BAM whose read groups all say PL:PACBIO passes (exit 0); one
     with PL:ILLUMINA, one with no @RG line, one with a PacBio and an Oxford
     Nanopore read group, and a file that is not a BAM stop (exit 3, or 1 for
     the file that cannot be read), each saying why;
  2. table, on a JSON shaped like Paraphase 4.1.0's: SMN1, SMN2 and the SMN
     exon 7-8 deletion from the smn1 region, a null SMN2 as 'no call'; PMS2,
     STRC and NCF1 from gene_cn; HBA1 and HBA2 counted from the haplotype
     names (a two-copy haplotype twice), the genotype with its number of
     alpha-globin genes, the 4.2 kb deletion named and HBA2 one copy short; a
     '--' allele said to be a large deletion that cannot be told apart; the
     RCCX copies with the CYP21A2 alleles; a region without gene_cn as the
     whole family's copy number;
  3. the variant table: one row per record a haplotype carries, with the
     haplotypes and their reads (ALT reads / depth); records no haplotype
     carries are left out;
  4. SMN1 from step 35's Parascopy table: from the row holding SMN1 exon 7,
     or from its homologous region on a row in SMN2; a quality under 20 or a
     filter other than PASS is not a call; the comparison says agree,
     disagree, or nothing when one side has no call.

Run: python3 tests/test_paraphase_report.py
"""
import gzip
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "paraphase_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))
import paraphase_report as pr  # noqa: E402

FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:600]}")
    if not ok:
        FAILS += 1


def bam(path, header):
    """A BAM with this header and no reference or read (one gzip member, as BGZF)."""
    text = header.encode()
    with gzip.open(path, "wb") as f:
        f.write(b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", 0))


def run(*args):
    p = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def rows(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


PARAPHASE_JSON = {
    "smn1": {"total_cn": 4, "genome_depth": 31.5, "region_depth": {"median": 120, "percentile80": 130},
             "region_specific_info": {"smn1_cn": 2, "smn2_cn": None, "smn_del78_cn": 0}},
    "hba": {"total_cn": 3, "genome_depth": 31.5, "region_depth": {"median": 90},
            "final_haplotypes": {"1121": "hba_hba1hap1", "1211": "hba_hba1hap2", "2111": "hba_hba2hap1",
                                 "2221": "hba_4p2delhap1", "1112": "hba_homology_hap1"},
            "two_copy_haplotypes": [],
            "region_specific_info": {"genotype": "-a/aa",
                                     "sv_called": {"hba_4p2delhap1": "4p2del,chr16:169742-174000"}}},
    "rccx": {"total_cn": 4, "region_specific_info": {"annotated_alleles": ["WT", "WT"], "phasing_success": True}},
    "pms2": {"total_cn": 4, "region_specific_info": {"gene_cn": 2}},
    "strc": {"total_cn": 4, "region_specific_info": {"gene_cn": 2, "intergenic_depth": 30}},
    "ncf1": {"total_cn": 6, "region_specific_info": {"gene_cn": 2, "gene_reads": 20, "pseudo_reads": 40}},
    "GBA": {"total_cn": 4, "genes_in_region": ["GBA1"], "region_specific_info": {"gene_cn": 2}},
    "OTOA": {"total_cn": 3, "genes_in_region": ["OTOA"], "region_specific_info": {}},
    "TNXB": {"total_cn": None, "genes_in_region": ["TNXB"], "region_specific_info": {}},
}

VCF = ("##fileformat=VCFv4.2\n"
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpms2_pms2hap1\tpms2_pms2hap2\tpms2_pms2clhap1\n"
       "chr7\t5970000\t.\tC\tT\t.\tPASS\tHPBOUND=1-2\tGT:DP:AD\t1:20:0,19\t0:22:22,0\t1:18:1,17\n"
       "chr7\t5980000\t.\tG\tA\t.\tPASS\tHPBOUND=1-2\tGT:DP:AD\t0:20:20,0\t0:22:22,0\t.:3:2,1\n")

PARASCOPY_HEAD = "chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n"


def main():
    work = tempfile.mkdtemp()
    try:
        # 1. check-bam
        b = os.path.join(work, "x.bam")
        bam(b, "@HD\tVN:1.6\n@RG\tID:a\tPL:PACBIO\tSM:S\n@RG\tID:b\tPL:PacBio\tSM:S\n")
        rc, out = run("check-bam", b)
        check("check-bam: two PL:PACBIO read groups pass", rc == 0 and "2 read group(s), all PL:PACBIO" in out, out)
        bam(b, "@HD\tVN:1.6\n@RG\tID:a\tPL:ILLUMINA\tSM:S\n")
        rc, out = run("check-bam", b)
        check("check-bam: PL:ILLUMINA stops with exit 3 and the reason (licence, step 35)",
              rc == 3 and "PL ILLUMINA, not PACBIO" in out and "licence" in out and "step 35" in out, (rc, out))
        bam(b, "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n")
        rc, out = run("check-bam", b)
        check("check-bam: no @RG line stops with exit 3", rc == 3 and "has no @RG line" in out, (rc, out))
        bam(b, "@RG\tID:a\tPL:PACBIO\n@RG\tID:b\tPL:ONT\n@RG\tID:c\n")
        rc, out = run("check-bam", b)
        check("check-bam: a PacBio and an ONT read group (and one with no PL) stop",
              rc == 3 and "PL ONT, none, not PACBIO" in out, (rc, out))
        with open(b, "wb") as f:
            f.write(b"not a bam at all")
        rc, out = run("check-bam", b)
        check("check-bam: a file that is not a BAM stops with exit 1", rc == 1 and "ERROR" in out, (rc, out))

        # 2 and 3. table
        j, vd = os.path.join(work, "S.paraphase.json"), os.path.join(work, "S_paraphase_vcfs")
        os.makedirs(vd)
        with open(j, "w") as f:
            json.dump(PARAPHASE_JSON, f)
        with open(os.path.join(vd, "S_pms2.vcf"), "w") as f:
            f.write(VCF)
        cn, var = os.path.join(work, "cn.tsv"), os.path.join(work, "var.tsv")
        rc, out = run("table", "--sample", "S", "--json", j, "--vcf-dir", vd, "--out-cn", cn, "--out-variants", var)
        check("table: exit 0", rc == 0, out)
        r = {(x["region"], x["gene"]): x for x in rows(cn)}
        g = lambda reg, gene, k="copy_number": r.get((reg, gene), {}).get(k)  # noqa: E731
        check("SMN1 2, SMN2 a no-call, the exon 7-8 deletion 0",
              (g("smn1", "SMN1"), g("smn1", "SMN2"), g("smn1", "SMN2", "status"), g("smn1", "SMN exon 7-8 deletion"))
              == ("2", ".", "no call", "0"), r)
        check("SMN2's no-call says why", "depth does not settle one copy or two" in (g("smn1", "SMN2", "detail") or ""))
        check("PMS2, STRC, NCF1 from gene_cn (2 each)",
              (g("pms2", "PMS2"), g("strc", "STRC"), g("ncf1", "NCF1")) == ("2", "2", "2"), r)
        check("HBA1 2, HBA2 1 (one HBA2 copy lost)", (g("hba", "HBA1"), g("hba", "HBA2")) == ("2", "1"), r)
        hd = g("hba", "HBA2", "detail") or ""
        check("HBA: genotype -a/aa with 3 alpha-globin genes, the 4.2 kb deletion named",
              "genotype -a/aa (3 alpha-globin genes)" in hd and "4.2 kb deletion (-a4.2): HBA2 lost" in hd, hd)
        check("the depth columns are kept", g("smn1", "SMN1", "region_depth") == "120"
              and g("smn1", "SMN1", "genome_depth") == "31.5", r.get(("smn1", "SMN1")))
        rd = g("rccx", "CYP21A2 (RCCX)", "detail") or ""
        check("RCCX: 4 copies, the CYP21A2 alleles", g("rccx", "CYP21A2 (RCCX)") == "4" and "WT / WT" in rd, rd)
        check("GBA1 from gene_cn, OTOA as the whole family",
              g("GBA", "GBA1") == "2" and "not of its pseudogenes" in (g("GBA", "GBA1", "detail") or "")
              and g("OTOA", "OTOA") == "3" and "whole gene family" in (g("OTOA", "OTOA", "detail") or ""), r)
        check("another region with no copy number is a no-call row",
              g("TNXB", "TNXB") == "." and g("TNXB", "TNXB", "status") == "no call", r)
        order = [x["region"] for x in rows(cn)]
        check("the named regions come first", order.index("TNXB") > max(order.index(k) for k in ("hba", "smn1", "OTOA")),
              order)
        v = rows(var)
        check("variants: the record two haplotypes carry, with their reads; not the one none carries",
              len(v) == 1 and v[0]["pos"] == "5970000" and v[0]["haplotypes"] == "pms2_pms2hap1,pms2_pms2clhap1"
              and v[0]["support"] == "pms2_pms2hap1:19/20,pms2_pms2clhap1:17/18" and v[0]["region"] == "pms2", v)

        data = json.loads(json.dumps(PARAPHASE_JSON))
        data["hba"]["region_specific_info"] = {"genotype": "--/aa", "sv_called": {}}
        data["hba"]["final_haplotypes"] = {"1": "hba_hba1hap1", "2": "hba_hba2hap1"}
        t = {(x[0], x[1]): x for x in pr.cn_table(data)}
        check("HBA '--/aa': one HBA1 and one HBA2, and a large deletion that cannot be told apart",
              "large deletion such as --SEA" in t[("hba", "HBA2")][4] and "(2 alpha-globin genes)" in t[("hba", "HBA2")][4]
              and t[("hba", "HBA1")][2] == "1" and t[("hba", "HBA2")][2] == "1", t.get(("hba", "HBA2")))
        data["hba"]["region_specific_info"] = {"genotype": "aa/aa", "sv_called": {}}
        data["hba"]["two_copy_haplotypes"] = ["hba_hba1hap1", "hba_hba2hap1"]
        t = {(x[0], x[1]): x for x in pr.cn_table(data)}
        check("HBA 'aa/aa' with identical copies: each two-copy haplotype counted twice, no deletion found",
              t[("hba", "HBA1")][2] == "2" and t[("hba", "HBA2")][2] == "2"
              and "no 3.7 kb or 4.2 kb deletion or duplication found" in t[("hba", "HBA2")][4], t.get(("hba", "HBA2")))

        # 4. SMN1 from Parascopy
        p = os.path.join(work, "smn.tsv")
        with open(p, "w") as f:
            f.write(PARASCOPY_HEAD + "chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t35,30\t"
                    "chr5:70049000-70077000:-\n")
        s1 = pr.parascopy_smn1(p)
        check("Parascopy SMN1 from the row in SMN1 (first psCN)", s1["cn"] == 2 and "quality 35" in s1["text"], s1)
        with open(p, "w") as f:
            f.write(PARASCOPY_HEAD + "chr5\t70049000\t70077000\tSMN1\tPASS\t4\t40\tPASS\t3,1\t30,25\t"
                    "chr5:70925000-70953000:-\n")
        s2 = pr.parascopy_smn1(p)
        check("Parascopy SMN1 from the homologous region of a row in SMN2 (second psCN)", s2["cn"] == 1, s2)
        with open(p, "w") as f:
            f.write(PARASCOPY_HEAD + "chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tPASS\t2,2\t12,30\t*\n")
        s3 = pr.parascopy_smn1(p)
        check("quality 12 is not a call", s3["cn"] is None and "not resolved" in s3["text"], s3)
        with open(p, "w") as f:
            f.write(PARASCOPY_HEAD + "chr5\t70925000\t70953000\tSMN1\tPASS\t4\t40\tFewReliable\t2,2\t35,30\t*\n")
        check("a filter other than PASS is not a call", pr.parascopy_smn1(p)["cn"] is None)
        check("compare: agree, disagree, no comparison",
              (pr.compare_smn1(2, s1), pr.compare_smn1(2, s2), pr.compare_smn1(None, s1), pr.compare_smn1(2, s3))
              == ("agree", "disagree", None, None))
    finally:
        shutil.rmtree(work)
    print("\nRESULT:", "ALL PASS" if FAILS == 0 else f"{FAILS} FAILED")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
