#!/usr/bin/env python3
"""Two hits in one gene across callers (bin/collect_summary.py two_hits): a
called deletion plus a small variant, on synthetic sample folders. Nothing
here describes a person: the calls are made up at real gene positions.

Checks:
  1. MUTYH: an AnnotSV exon deletion plus a heterozygous ClinVar P/LP hit
     elsewhere in the gene gives a "possibly biallelic" line, and the ACMG SF
     tier now reports the hit (MUTYH needs two variants: the deletion is one);
  2. the same two with read-backed phase on different copies (GT with "|",
     same PS) give "biallelic";
  3. the same two in cis (same copy, same PS) give no two-hit line, and the SF
     tier keeps the one hit apart;
  4. TAR syndrome: step 35's 1q21.1 deletion plus the RBM8A low-expression
     allele (rs139428292, assets/carrier_genes.tsv) gives the TAR line;
  5. a hit called 1/1 inside a called deletion is shown as hemizygous, in the
     carrier list and in the SF tier;
  6. SMN1: one SMN1 copy (Parascopy) plus an SMN1 hit, and HBA: three alpha
     genes plus an HBA2 hit, each give a "possibly biallelic" line;
  7. without any deletion caller's output the hit gives no line, and the
     summary names the callers it could not read;
  8. both reports print the two-hit line and the hemizygous genotype.
"""
import gzip
import os
import shutil
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import collect_summary  # noqa: E402
import render_report  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + name + ("" if ok else " -- " + str(detail)[:1500]))
    if not ok:
        FAILS.append(name)


def put(path, text, gz=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with (gzip.open(path, "wt") if gz else open(path, "w")) as f:
        f.write(text)
    return path


VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
ANNOTSV_COLS = ("SV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tCaller\tExon_edge\t"
                "Exon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\tFORMAT\tS\n")
# MUTYH spans chr1:45329163-45340893 (assets/carrier_genes.tsv).
SNV = ("chr1", 45332000, "C", "T")
DEL = ("chr1", 45335000, 45337000)
INSIDE = ("chr1", 45336000, "G", "A")
# RBM8A c.-21G>A (rs139428292), on the plus strand C>T.
TAR_ALLELE = ("chr1", 145927447, "C", "T")


def hit(v, gene, gt, ps=None):
    chrom, pos, ref, alt = v
    fmt, val = ("GT:GQ:DP:PS", f"{gt}:50:30:{ps}") if ps else ("GT:GQ:DP", f"{gt}:50:30")
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tGENEINFO={gene}:1;CLNSIG=Pathogenic;"
            f"CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\t{fmt}\t{val}\n")


def sv(span, genes, gt="0/1", ps=None):
    chrom, a, b = span
    fmt, val = ("GT:PS", f"{gt}:{ps}") if ps else ("GT", gt)
    return f"{chrom}\t{a}\t{b}\t-{b - a}\tDEL\t{genes}\t4\tmanta\tyes\t.\t\t\tno\t{fmt}\t{val}\n"


SMN_COLS = "chrom\tstart\tend\tagCN\tpsCN\tpsCN_qual\tpsCN_filter\tagCN_qual\tagCN_filter\thomologous_regions\n"
INTERVAL_COLS = "locus\tname\tchrom\tstart\tend\tstatus\tcopies\texpected\tchange\testimate\treason\n"


def step35(d, smn1="2", intervals="", hba_genes=None):
    """Step 35's tables: Parascopy's SMN table (SMN1 exon 7 at chr5:70951946,
    SMN2's at chr5:70076526), the interval table and the HBA count."""
    put(f"{d}/paralogs/S_smn_copy_number.tsv", SMN_COLS
        + f"chr5\t70951000\t70952000\t4\t{smn1},2\t40,40\tPASS\t40\tPASS\tchr5:70076000-70077000\n")
    if intervals:
        put(f"{d}/paralogs/S_cn_intervals.tsv", INTERVAL_COLS + intervals)
    if hba_genes is not None:
        put(f"{d}/paralogs/S_hba_copy_number.tsv",
            f"metric\tvalue\nstatus\tcalled\nalpha_globin_genes\t{hba_genes}\nestimate\t{hba_genes}.02\n"
            "message\tone gene lost on one chromosome\n")


def sample(work, name, hits, annotsv=None):
    d = os.path.join(work, name, "S")
    put(f"{d}/clinvar/S_clinvar_hits.vcf", VCF_HEAD + "".join(hits))
    if annotsv is not None:
        put(f"{d}/annotsv/S_sv_class4_5_checked.tsv", ANNOTSV_COLS + annotsv)
    return d


def carrier(d):
    return collect_summary.sec_carrier(d, "S")[1]


def two(d, gene):
    return [x for x in (carrier(d).get("two_hits") or {}).get("genes", []) if x["gene"] == gene]


def sf(d):
    return collect_summary.acmg_sf_tier(d, "S", None)


def main():
    work = tempfile.mkdtemp(prefix="two_hits_")
    try:
        # 1. a MUTYH exon deletion plus a heterozygous hit: possibly biallelic
        d = sample(work, "mutyh", [hit(SNV, "MUTYH", "0/1")], sv(DEL, "MUTYH"))
        t = two(d, "MUTYH")
        check("MUTYH deletion + het hit: one possibly biallelic line",
              len(t) == 1 and t[0]["verdict"] == "possibly biallelic", t)
        check("MUTYH deletion + het hit: the line names the deletion and the hit, phase unknown",
              t and "DEL chr1:45335000-45337000" in t[0]["deletion"] and t[0]["variant"] == "chr1:45332000 C>T"
              and "possibly in trans" in t[0]["why"], t)
        s = sf(d)
        rep = [x for x in s["clinvar_hits"] if x["gene"] == "MUTYH"]
        check("MUTYH deletion + het hit: the SF tier reports the hit, counted with the deletion",
              len(rep) == 1 and "possibly biallelic with DEL" in (rep[0].get("two_hit") or ""), rep)
        check("MUTYH deletion + het hit: nothing left apart",
              not [x for x in s["not_reported"] if x["gene"] == "MUTYH"], s["not_reported"])

        # 2. read-backed phase on different copies: biallelic
        d = sample(work, "trans", [hit(SNV, "MUTYH", "0|1", ps=45330000)], sv(DEL, "MUTYH", "1|0", ps=45330000))
        t = two(d, "MUTYH")
        check("phased in trans: biallelic", len(t) == 1 and t[0]["verdict"] == "biallelic"
              and "in trans" in t[0]["why"], t)
        # a different phase set links nothing: possibly
        d = sample(work, "otherps", [hit(SNV, "MUTYH", "0|1", ps=45330000)], sv(DEL, "MUTYH", "1|0", ps=45000000))
        t = two(d, "MUTYH")
        check("phased in two different blocks: possibly biallelic", t and t[0]["verdict"] == "possibly biallelic", t)

        # 3. in cis: no line, and the SF tier keeps the one hit apart
        d = sample(work, "cis", [hit(SNV, "MUTYH", "1|0", ps=45330000)], sv(DEL, "MUTYH", "1|0", ps=45330000))
        check("phased in cis: no two-hit line", two(d, "MUTYH") == [], two(d, "MUTYH"))
        s = sf(d)
        check("phased in cis: the SF tier lists the one hit apart, as before",
              not [x for x in s["clinvar_hits"] if x["gene"] == "MUTYH"]
              and [x for x in s["not_reported"] if x["gene"] == "MUTYH"], s)

        # 4. TAR: the 1q21.1 deletion of step 35 plus the RBM8A low-expression allele
        tar_row = ("1q21.1_TAR_syndrome_proximal_BP2-BP3\t1q21.1_TAR_syndrome_proximal_BP2-BP3_t1\tchr1\t145686999\t"
                   "146048495\tcalled\t1\t2\tdeletion\t1.04\t.\n")
        d = sample(work, "tar", [hit(TAR_ALLELE, "RBM8A", "1/1")])
        step35(d, intervals=tar_row)
        c = carrier(d)
        t = [x for x in c["two_hits"]["genes"] if x["gene"] == "RBM8A"]
        check("TAR: the 1q21.1 deletion + the RBM8A allele give one RBM8A line",
              len(t) == 1 and t[0]["verdict"] == "biallelic", c.get("two_hits"))
        line = [x for x in render_report.carrier_lines(c) if "RBM8A" in x and "both copies" in x]
        check("TAR: the line names the condition, the deletion and the allele",
              len(line) == 1 and "Thrombocytopenia-absent radius syndrome" in line[0]
              and "1q21.1_TAR_syndrome" in line[0] and "rs139428292" in line[0], line)
        check("TAR: the allele alone stays set aside, never a carrier finding",
              [x for x in c["set_aside"] if x["gene"] == "RBM8A"] and not
              [x for x in c["carriers"] if x["gene"] == "RBM8A"], c)
        # the allele without the deletion: no line
        d = sample(work, "tar_alone", [hit(TAR_ALLELE, "RBM8A", "0/1")])
        step35(d, intervals=tar_row.replace("\t1\t2\tdeletion\t1.04", "\t2\t2\tno change\t2.01"))
        check("TAR allele without a deletion: no line", two(d, "RBM8A") == [], two(d, "RBM8A"))

        # 5. a hit called 1/1 inside the deletion: hemizygous
        d = sample(work, "hemi", [hit(INSIDE, "MUTYH", "1/1")], sv(DEL, "MUTYH"))
        c = carrier(d)
        item = [x for x in c["carriers"] if x["gene"] == "MUTYH"]
        check("1/1 inside a deletion: the carrier list says hemizygous",
              len(item) == 1 and item[0]["genotype"] == "hemizygous" and "inside DEL" in item[0].get("hemizygous", ""),
              item)
        t = two(d, "MUTYH")
        check("1/1 inside a deletion: biallelic, the variant on the one copy left",
              t and t[0]["verdict"] == "biallelic" and t[0]["genotype"] == "hemizygous", t)
        rep = [x for x in sf(d)["clinvar_hits"] if x["gene"] == "MUTYH"]
        check("1/1 inside a deletion: the SF tier says hemizygous", rep and rep[0]["genotype"] == "hemizygous", rep)
        d = sample(work, "hom_out", [hit(SNV, "MUTYH", "1/1")], sv(DEL, "MUTYH"))
        item = [x for x in carrier(d)["carriers"] if x["gene"] == "MUTYH"]
        check("1/1 outside the deletion stays hom", item and item[0]["genotype"] == "hom", item)

        # 6. SMN1 and HBA
        d = sample(work, "smn", [hit(("chr5", 70951950, "G", "A"), "SMN1", "0/1")])
        step35(d, smn1="1")
        t = two(d, "SMN1")
        check("one SMN1 copy + an SMN1 hit: a line", t and "one SMN1 copy" in t[0]["deletion"], t)
        d = sample(work, "smn2", [hit(("chr5", 70951950, "G", "A"), "SMN1", "0/1")])
        step35(d, smn1="2")
        check("two SMN1 copies + an SMN1 hit: no line", two(d, "SMN1") == [], two(d, "SMN1"))
        d = sample(work, "hba", [hit(("chr16", 173600, "C", "T"), "HBA2", "0/1")])
        step35(d, hba_genes="3")
        t = two(d, "HBA2")
        check("three alpha genes + an HBA2 hit: possibly biallelic",
              t and t[0]["verdict"] == "possibly biallelic" and "3 of 4 alpha-globin genes" in t[0]["deletion"], t)
        d = sample(work, "hba4", [hit(("chr16", 173600, "C", "T"), "HBA2", "0/1")])
        step35(d, hba_genes="4")
        check("four alpha genes + an HBA2 hit: no line", two(d, "HBA2") == [], two(d, "HBA2"))

        # 7. no deletion caller: no line, and the gap is named
        d = sample(work, "none", [hit(SNV, "MUTYH", "0/1")])
        t = carrier(d)["two_hits"]
        check("no deletion caller: no line", t["genes"] == [], t)
        check("no deletion caller: the callers not read are named",
              any("AnnotSV" in x for x in t["deletions_not_read"])
              and any("step 35" in x for x in t["deletions_not_read"]), t)
        check("no deletion caller: the report line says so",
              any("not read: AnnotSV" in x for x in render_report.carrier_lines(carrier(d))),
              render_report.carrier_lines(carrier(d)))

        # 8. both reports
        d = sample(work, "report", [hit(INSIDE, "MUTYH", "1/1"), hit(("chr1", 45330500, "A", "G"), "MUTYH", "0/1")],
                   sv(DEL, "MUTYH"))
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        check("text report: the two-hit line", "Biallelic: MUTYH" in txt, txt[-3000:])
        check("HTML report: the two-hit line", "Biallelic: MUTYH" in html)
        check("text report: the hemizygous genotype", "45336000 G&gt;A" not in txt and "hemizygous (called homozygous "
              "inside DEL" in txt, txt[-3000:])
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print(f"\n{len(FAILS)} failed" if FAILS else "\nall passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
