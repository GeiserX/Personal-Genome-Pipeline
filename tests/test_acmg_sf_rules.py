#!/usr/bin/env python3
"""The ACMG SF v3.3 reporting rules (assets/acmg_sf.tsv, bin/collect_summary.py
acmg_sf_tier) and the one-person carrier list (assets/carrier_genes.tsv), on
synthetic sample folders.

Checks:
  1. the catalog: exactly the 84 genes of SF v3.3 Table 1 (Lee et al., Genet
     Med 2025;27(8):101454), one row each; the X-linked rows are GLA, OTC and
     ABCD1 with Table 1's rule each (GLA and OTC: hemi, het and homozygous;
     ABCD1: hemi, homozygous or two het); the recessive genes need two
     variants; HFE counts p.C282Y homozygotes only;
  2. MUTYH: one heterozygous ClinVar hit is not a secondary finding (it is
     listed apart with the reason), two are, a homozygote is;
  3. HFE: p.C282Y heterozygous is not an SF, homozygous is, a homozygous
     p.His63Asp is not;
  4. X-linked: a female heterozygote (0/1) in GLA and in OTC is reported, in
     ABCD1 it is not; a hemizygous ABCD1 call (GT 1) is, and so are two ABCD1
     heterozygous variants;
  5. structural: an AnnotSV class 5 exon deletion in BRCA1 is in the SF
     section, a class 5 deletion of the 3' end of EPCAM counts for MSH2, and a
     class 5 deletion in MUTYH (no SV rule) is not;
  6. severity: BTD p.Asp444His homozygous is reported as partial, not
     profound; a BTD homozygote without it carries no partial note;
  7. the carrier list of one sample, with no partner: a heterozygous CFTR hit
     is a carrier line, a HEXA pseudodeficiency allele (c.739C>T) is set
     aside and makes no carrier line, another HEXA hit does; the
     couple screen sets the same allele aside;
  8. both reports print the SF verdicts, the reasons and the carrier list;
  9. without the assets folder the ClinVar hits and the clinical counts stay,
     and both lists say why they were not made.
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


# ACMG SF v3.3 Table 1, gene by gene.
SF_V33 = """APC RET BRCA1 BRCA2 PALB2 SDHD SDHAF2 SDHC SDHB MAX TMEM127 BMPR1A SMAD4 TP53 MLH1 MSH2 MSH6 PMS2
MEN1 MUTYH NF2 STK11 PTEN RB1 TSC1 TSC2 VHL WT1 FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYH11 PKP2 DSP DSC2 TMEM43
DSG2 RYR2 CASQ2 TRDN TNNT2 LMNA FLNC TTN BAG3 DES RBM20 TNNC1 PLN COL3A1 LDLR APOB PCSK9 MYH7 MYBPC3
TNNI3 TPM1 MYL3 ACTC1 PRKAG2 MYL2 KCNQ1 KCNH2 SCN5A CALM1 CALM2 CALM3 BTD CYP27A1 GLA HFE OTC GAA ABCD1
ACVRL1 ENG TTR RYR1 CACNA1S HNF1A RPE65 ATP7B""".split()
# Table 1 "Variants to report" of the X-linked genes, and the rule each must be.
XL_RULES = {"GLA": ("All hemi, het, homozygous P and LP", "het_female_too"),
            "OTC": ("All hemi, het, homozygous P and LP", "het_female_too"),
            "ABCD1": ("All hemi, homozygous or two het. P and/or LP", "hemizygous")}
AR_TWO = {"MUTYH", "CASQ2", "TRDN", "BTD", "CYP27A1", "GAA", "RPE65", "ATP7B"}

VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def hit(chrom, pos, ref, alt, gene, gt, sig="Pathogenic"):
    return (f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tGENEINFO={gene}:1;CLNSIG={sig};"
            f"CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts\tGT:GQ:DP\t{gt}:50:30\n")


def sample(work, name, hits, annotsv=None):
    d = os.path.join(work, name, "S")
    put(f"{d}/clinvar/S_clinvar_hits.vcf", VCF_HEAD + "".join(hits))
    if annotsv is not None:
        put(f"{d}/annotsv/S_sv_class4_5_checked.tsv", ANNOTSV_COLS + annotsv)
    return d


ANNOTSV_COLS = ("SV_chrom\tSV_start\tSV_end\tSV_length\tSV_type\tGene_name\tACMG_class\tCaller\tExon_edge\t"
                "Exon_edge_detail\tPopulation_AF\tPopulation_source\tCommon\n")


def sv(chrom, a, b, genes, cls="5", kind="DEL"):
    return f"{chrom}\t{a}\t{b}\t-{b - a}\t{kind}\t{genes}\t{cls}\tmanta\tno\t.\t\t\tno\n"


def sf_of(d):
    return collect_summary.collect("S", d)["sections"]["clinical"]["data"].get("acmg_sf") or {}


def sf_of_clinvar(d):
    """The SF tier of a folder with ClinVar hits only (no clinical VCF): the
    Nextflow layout gives sec_clinical the clinical VCF, so the tier is read
    straight from the hits here."""
    return collect_summary.acmg_sf_tier(d, "S", None)


def reported(sf, gene):
    return [x for x in sf.get("clinvar_hits", []) if x["gene"] == gene]


def apart(sf, gene):
    return [x for x in sf.get("not_reported", []) if x["gene"] == gene]


def main():
    work = tempfile.mkdtemp(prefix="acmg_sf_")
    try:
        # 1. the catalog
        p = os.path.join(REPO, "assets", "acmg_sf.tsv")
        try:
            rows = collect_summary.asset_rows("acmg_sf.tsv")
            rules = collect_summary.acmg_sf_rules()
        except Exception as e:  # the pre-fix state: no catalog, no loader
            rows, rules = [], {}
            check("catalog: assets/acmg_sf.tsv loads", False, f"{type(e).__name__}: {e}")
        genes = [r["gene"] for r in rows]
        check("catalog: exactly 84 rows, one per gene", len(genes) == 84 and len(set(genes)) == 84, len(genes))
        check("catalog: the genes are SF v3.3 Table 1's", set(genes) == set(SF_V33) and len(SF_V33) == 84,
              sorted(set(genes) ^ set(SF_V33)))
        check("catalog: every row has an HGNC id", all(r.get("hgnc_id", "").startswith("HGNC:") for r in rows))
        xl = {g for g, r in rules.items() if r["inheritance"] == "XL"}
        check("catalog: the X-linked rows are GLA, OTC and ABCD1", xl == set(XL_RULES), sorted(xl))
        for g, (text, rule) in XL_RULES.items():
            r = rules.get(g, {})
            check(f"catalog: {g} is '{text}', rule {rule}",
                  (r.get("table_text"), r.get("report_if")) == (text, rule), r)
        two = {g for g, r in rules.items() if r["report_if"] == "biallelic"}
        check("catalog: the recessive two-variant genes", two == AR_TWO, sorted(two ^ AR_TWO))
        hfe = rules.get("HFE", {})
        check("catalog: HFE counts p.C282Y (chr6:26092913 G>A) homozygotes only",
              hfe.get("report_if") == "hom_variant" and list(hfe.get("variants", {})) == [("chr6", 26092913, "G", "A")],
              hfe)
        check("catalog: BRCA1, BRCA2, MLH1, MSH2 and LDLR take class 4/5 SVs",
              {g for g, r in rules.items() if r["sv_cn"]} == {"BRCA1", "BRCA2", "MLH1", "MSH2", "LDLR"},
              sorted(g for g, r in rules.items() if r["sv_cn"]))
        # A dropped gene makes the count check fail: the check can go red.
        dropped = os.path.join(work, "dropped")
        os.makedirs(dropped)
        if os.path.isfile(p):
            with open(p) as f, open(os.path.join(dropped, "acmg_sf.tsv"), "w") as g:
                g.write("".join(x for x in f if not x.startswith("TTN\t")))
            keep = collect_summary.ASSETS_DIR
            collect_summary.ASSETS_DIR = dropped
            try:
                n = len(collect_summary.acmg_sf_rules())
            finally:
                collect_summary.ASSETS_DIR = keep
            check("catalog control: a copy without TTN has 83 genes, not 84", n == 83, n)

        # 2. MUTYH: chr1 MUTYH p.Gly396Asp and p.Tyr179Cys (GRCh38)
        g396 = hit("chr1", 45332803, "C", "T", "MUTYH", "0/1")
        y179 = hit("chr1", 45334483, "T", "C", "MUTYH", "0/1")
        sf = sf_of_clinvar(sample(work, "mutyh_het", [g396]))
        check("MUTYH: one heterozygous variant is not a secondary finding", reported(sf, "MUTYH") == [], sf)
        check("MUTYH: it is listed apart with the reason",
              [x["reason"] for x in apart(sf, "MUTYH")] == ["one heterozygous variant; MUTYH is reported with "
                                                            "P and LP (2 variants)"], apart(sf, "MUTYH"))
        sf = sf_of_clinvar(sample(work, "mutyh_two", [g396, y179]))
        check("MUTYH: two heterozygous variants are a secondary finding (phase unknown)",
              len(reported(sf, "MUTYH")) == 2 and not apart(sf, "MUTYH"), sf)
        sf = sf_of_clinvar(sample(work, "mutyh_hom", [hit("chr1", 45332803, "C", "T", "MUTYH", "1/1")]))
        check("MUTYH: a homozygote is a secondary finding", len(reported(sf, "MUTYH")) == 1, sf)

        # 3. HFE
        c282y_het = hit("chr6", 26092913, "G", "A", "HFE", "0/1")
        sf = sf_of_clinvar(sample(work, "hfe_het", [c282y_het]))
        check("HFE: p.C282Y heterozygous is not an SF", reported(sf, "HFE") == [] and len(apart(sf, "HFE")) == 1, sf)
        sf = sf_of_clinvar(sample(work, "hfe_hom", [hit("chr6", 26092913, "G", "A", "HFE", "1/1")]))
        check("HFE: p.C282Y homozygous is an SF", len(reported(sf, "HFE")) == 1, sf)
        sf = sf_of_clinvar(sample(work, "hfe_h63d", [hit("chr6", 26090951, "C", "G", "HFE", "1/1")]))
        check("HFE: p.His63Asp homozygous is not an SF", reported(sf, "HFE") == [], sf)

        # 4. X-linked, a female heterozygote (diploid 0/1) and a hemizygous call (GT 1)
        sf = sf_of_clinvar(sample(work, "xl_female", [
            hit("chrX", 101398000, "C", "T", "GLA", "0/1"), hit("chrX", 38400000, "G", "A", "OTC", "0/1"),
            hit("chrX", 153725000, "C", "T", "ABCD1", "0/1")]))
        check("X-linked: a female heterozygote in GLA is reported (Table 1: hemi, het, homozygous)",
              len(reported(sf, "GLA")) == 1, sf)
        check("X-linked: a female heterozygote in OTC is reported (Table 1: hemi, het, homozygous)",
              len(reported(sf, "OTC")) == 1, sf)
        check("X-linked: a female heterozygote in ABCD1 is not (Table 1: hemi, homozygous or two het)",
              reported(sf, "ABCD1") == [] and len(apart(sf, "ABCD1")) == 1, sf)
        sf = sf_of_clinvar(sample(work, "xl_male", [hit("chrX", 153725000, "C", "T", "ABCD1", "1")]))
        check("X-linked: a hemizygous ABCD1 variant is reported", len(reported(sf, "ABCD1")) == 1, sf)
        sf = sf_of_clinvar(sample(work, "xl_two", [hit("chrX", 153725000, "C", "T", "ABCD1", "0/1"),
                                                    hit("chrX", 153726000, "G", "A", "ABCD1", "0/1")]))
        check("X-linked: two heterozygous ABCD1 variants are reported", len(reported(sf, "ABCD1")) == 2, sf)

        # 5. structural: BRCA1 exon deletion, EPCAM 3' end, MUTYH (no SV rule)
        d = sample(work, "sv", [], sv("17", 43070000, 43080000, "BRCA1") + sv("2", 47380000, 47390000, "EPCAM")
                   + sv("1", 45330000, 45336000, "MUTYH"))
        sf = sf_of_clinvar(d)
        st = {(x["gene"], x.get("via")) for x in sf.get("structural", [])}
        check("structural: a class 5 BRCA1 exon deletion is in the SF section", ("BRCA1", None) in st, sf)
        check("structural: a class 5 deletion of EPCAM's 3' end counts for MSH2",
              any(g == "MSH2" and v and "EPCAM" in v for g, v in st), st)
        check("structural: a class 5 deletion in MUTYH (no SV rule) is not", not any(g == "MUTYH" for g, _ in st), st)

        # 6. severity: BTD p.Asp444His (chr3:15645186 G>C)
        sf = sf_of_clinvar(sample(work, "btd_d444h", [hit("chr3", 15645186, "G", "C", "BTD", "1/1")]))
        btd = reported(sf, "BTD")
        check("BTD: D444H/D444H is reported as partial, not profound",
              len(btd) == 1 and (btd[0].get("severity") or "").startswith("partial, not profound"), sf)
        sf = sf_of_clinvar(sample(work, "btd_profound", [hit("chr3", 15647000, "G", "A", "BTD", "1/1")]))
        btd = reported(sf, "BTD")
        check("BTD: a homozygote without D444H carries no partial note",
              len(btd) == 1 and not btd[0].get("severity"), sf)

        # 7. the carrier list of one sample, no partner; HEXA c.739C>T pseudodeficiency
        d = sample(work, "carrier", [hit("chr7", 117559590, "ATCT", "A", "CFTR", "0/1"),
                                     hit("chr15", 72350584, "G", "A", "HEXA", "0/1"),
                                     hit("chr15", 72346579, "G", "GGATA", "HEXA", "0/1")])
        summ = collect_summary.collect("S", d)
        car = summ["sections"].get("carrier", {}).get("data") or {}
        lines = {(x["gene"], x["variant"]) for x in car.get("carriers", [])}
        check("carrier list: one sample, no partner, a CFTR carrier line", ("CFTR", "chr7:117559590 ATCT>A") in lines, car)
        check("carrier list: the HEXA pseudodeficiency allele makes no carrier line",
              ("HEXA", "chr15:72350584 G>A") not in lines
              and [x["kind"] for x in car.get("set_aside", [])] == ["pseudodeficiency"], car)
        check("carrier list: the other HEXA hit does", ("HEXA", "chr15:72346579 G>GGATA") in lines, car)
        import couple_report
        genes = couple_report.read_genes(os.path.join(REPO, "assets", "carrier_genes.tsv"))
        aside = []
        try:
            f, _ = couple_report.clinvar_findings(f"{d}/clinvar/S_clinvar_hits.vcf", genes, "female", {}, aside)
        except TypeError as e:  # the pre-fix signature: no set-aside list
            f = [{"gene": str(e)}]
        check("couple screen: the same pseudodeficiency allele is no carrier finding",
              sorted(x["gene"] for x in f) == ["CFTR", "HEXA"] and len(aside) == 1 and "pseudodeficiency" in aside[0],
              (f, aside))

        # 8. both reports
        d = sample(work, "report", [g396, c282y_het, hit("chr3", 15645186, "G", "C", "BTD", "1/1"),
                                    hit("chr7", 117559590, "ATCT", "A", "CFTR", "0/1")],
                   sv("17", 43070000, 43080000, "BRCA1"))
        put(f"{d}/clinical/S_clinical.vcf.gz", "##fileformat=VCFv4.2\n" + VCF_HEAD.split("\n", 1)[1], gz=True)
        summ = collect_summary.collect("S", d)
        txt, html = render_report.text_report(summ), render_report.html_report(summ)
        for what, needle in (("MUTYH not reported, with the reason", "MUTYH chr1:45332803 C>T het (ClinVar): not "
                              "reported, one heterozygous variant"),
                             ("HFE not reported", "HFE chr6:26092913 G>A het (ClinVar): not reported"),
                             ("BTD partial", "partial, not profound"),
                             ("the BRCA1 deletion", "DEL chr17:43070000-43080000"),
                             ("the carrier list", "Carrier list (")):
            check(f"text report: {what}", needle in txt, [x for x in txt.splitlines() if "ACMG" in x or "Carrier" in x])
            check(f"HTML report: {what}", needle in html.replace("&gt;", ">"))

        # 9. without the assets folder
        keep = collect_summary.ASSETS_DIR
        collect_summary.ASSETS_DIR = os.path.join(work, "no-assets")
        collect_summary._ASSETS.clear()
        try:
            summ = collect_summary.collect("S", d)
        finally:
            collect_summary.ASSETS_DIR = keep
            collect_summary._ASSETS.clear()
        cv, cl = summ["sections"]["clinvar"], summ["sections"]["clinical"]
        check("no assets: the ClinVar hits stay", cv["state"] == "ok" and cv["data"]["count"] == 4, cv["state"])
        check("no assets: both lists say why they were not made",
              "acmg_sf.tsv" in (cl["data"].get("acmg_sf") or {}).get("unavailable", "")
              and "carrier_genes.tsv" in summ["sections"]["carrier"]["data"].get("unavailable", ""),
              (cl, summ["sections"]["carrier"]))
        txt = render_report.text_report(summ)
        check("no assets: the text report says so", "Secondary findings: not assessed" in txt
              and "Carrier list: not made" in txt)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    print("\nRESULT:", "ALL PASS" if not FAILS else "%d FAILED" % len(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
