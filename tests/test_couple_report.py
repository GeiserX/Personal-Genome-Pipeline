#!/usr/bin/env python3
"""bin/couple_report.py (step 39 and the Nextflow COUPLE_CARRIER process) on
synthetic step outputs at the real coordinates of assets/carrier_genes.tsv.

Checks:
  1. the carrier list: every ACMG tier 3 gene (97 autosomal recessive, 16
     X-linked) is in it, the ACOG genes are marked, every row has a GRCh38
     span on chr1 to chrX, and a chrX gene outside the PARs is X-linked;
  2. one variant in the same recessive gene in each partner is one shared
     gene (1 in 4 per pregnancy), a variant in another gene is one partner
     only, and a hit in a gene outside the list (APOE) is left out;
  3. a Delly deletion over an exon of a gene where the other partner has a
     SNV counts as shared; a deletion inside an intron, one that failed
     Delly's filter and one without the exon list do not;
     a Delly deletion longer than 1,000 bp counts only when its RDCN (read
     depth) is below a normal region's copies: a 50 Mb PASS het deletion
     with RDCN 2 and a 2 kb one without an estimate are set aside, the same
     50 Mb deletion with RDCN 1 counts, a 600 bp precise one counts on split
     reads alone, and on a male's chrX RDCN 1 is set aside and RDCN 0 counts;
  4. a homozygous hit gives 1 in 2, two het hits in one partner say the
     chance can be higher;
  5. X-linked: a female's FMR1 premutation (ExpansionHunter) gives sons 1 in
     2 with the premutation note; a male's DMD hit makes every daughter a
     carrier; FXN expansions in both partners are shared;
  6. SMN1: one copy in each partner (Parascopy) is shared; an unresolved call
     is listed as not assessed;
  7. a partner without Delly, Parascopy or ExpansionHunter is named in the
     residual-risk lines; the partner tables are the couple table;
  8. red controls: an emptied gene list stops the screen (exit 2) and writes
     no table; a pedigree without partners exits 1; a couple where one
     partner has no ClinVar screen is not screened;
  9. the same inputs give the same bytes;
 10. genes: the list is rebuilt from a small Orphanet and GENCODE extract,
     keeping the ACMG rows, adding an autosomal and an X-linked Orphanet gene,
     a PAR gene as autosomal, leaving out chrY, and stopping on a panel gene
     GENCODE does not have.

Run: python3 tests/test_couple_report.py
"""
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "couple_report.py")
GENES = os.path.join(REPO, "assets", "carrier_genes.tsv")
PAR = os.path.join(REPO, "assets", "par_grch38.bed")
sys.path.insert(0, os.path.join(REPO, "bin"))
import couple_report as cr  # noqa: E402

FAILS = []

ACMG = """HBB XPC TYR CYP21A2 PAH CFTR TNXB HEXA GJB2 DHCR7 ATP7B ASPA ACADM NPHS1 PMM2 FKTN SLC26A4 ERCC2
DYNC2H1 CEP290 GBE1 GAA CHRNE G6PC1 OCA2 COL7A1 ABCC8 ALDOB FANCC GRIP1 BCKDHB ANO10 NAGA SMPD1 USH2A MMUT CPT2
AHI1 DHDDS SLC19A3 GALT CYP11A1 TF MMACHC ABCA3 GBA1 MCOLN1 GNPTAB AGA PCDH15 FAH AIRE BBS2 CYP27A1 CCDC88C FMO3
TMEM216 CNGB3 MCPH1 SLC37A4 PRF1 SCO2 AGXT ACADVL ASL EVC2 ARSA MVK PKHD1 BTD ALPL BBS1 CLCN1 CYP27B1 POLG MCCC2
MLC1 ACAT1 CC2D2A SLC26A2 CBS LRP2 IDUA FKRP RNASEH2B RARS2 HBA1 HBA2 SMN1 HPS1 HPS3 ELP1 FXN DLD NEB CLRN1
BLM""".split()
ACMG_XL = "ABCD1 AFF2 ARX DMD F8 F9 FMR1 GLA L1CAM MID1 NR0B1 OTC PLP1 RPGR RS1 SLC6A8".split()
ACOG = "CFTR SMN1 HBB HBA1 HBA2 FMR1 HEXA ASPA ELP1 BLM ABCC8 FANCC GBA1 G6PC1 TMEM216 BCKDHB MCOLN1 SMPD1 PCDH15 CLRN1".split()


def check(desc, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + desc + ("" if ok else "  " + str(detail)))
    if not ok:
        FAILS.append(desc)


def run(*args):
    p = subprocess.run([sys.executable, SCRIPT] + list(args), capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def put(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    return path


GENE_ROWS = cr.read_genes(GENES)


def mid(gene, off=0):
    g = GENE_ROWS[gene]
    return g["chrom"], (g["start"] + g["end"]) // 2 + off


HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def hit(gene, gt="0/1", off=0, gid="1", rev="criteria_provided,_multiple_submitters,_no_conflicts", at=None):
    c, p = at or mid(gene, off)
    return (f"{c}\t{p}\t9000{off % 97}\tA\tG\t.\tPASS\tGENEINFO={gene}:{gid};CLNSIG=Pathogenic;"
            f"CLNREVSTAT={rev}\tGT\t{gt}\n")


def deletion(chrom, start, end, gt="0/1", filt="PASS"):
    return f"{chrom}\t{start}\tDEL0\tN\t<DEL>\t.\t{filt}\tSVTYPE=DEL;END={end}\tGT:FT\t{gt}:PASS\n"


def deldepth(chrom, start, end, rdcn, gt="0/1", precise="IMPRECISE"):
    """A Delly 2.x deletion with its read-depth copy number (FORMAT RDCN)."""
    return (f"{chrom}\t{start}\tDEL1\tN\t<DEL>\t.\tPASS\t{precise};SVTYPE=DEL;END={end}\t"
            f"GT:FT:RDCN\t{gt}:PASS:{rdcn}\n")


def eh(locus, alleles):
    """An ExpansionHunter record; alleles: the repeat count of each allele (invented)."""
    gt = "/".join("1" for _ in alleles)
    sizes = "/".join(str(a) for a in alleles)
    return f"chrX\t1\t.\tC\t<STR>\t.\tPASS\tVARID={locus};REPID={locus}\tGT:REPCN\t{gt}:{sizes}\n"


def smn(cn, qual="45", filt="PASS"):
    return ("chrom\tstart\tend\tlocus\tagCN_filter\tagCN\tagCN_qual\tpsCN_filter\tpsCN\tpsCN_qual\thomologous_regions\n"
            f"chr5\t70951000\t70953000\tSMN1\tPASS\t{cn + 2}\t40\t{filt}\t{cn},2\t{qual},40\tchr5:70076000-70078000\n")


def rows(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f]


def main():
    # 1. the carrier list
    acmg = {g for g, r in GENE_ROWS.items() if "ACMG" in r["lists"]}
    check("all 113 ACMG tier 3 genes are in the list", acmg == set(ACMG) | set(ACMG_XL),
          sorted(acmg ^ (set(ACMG) | set(ACMG_XL))))
    check("97 of them autosomal recessive, 16 X-linked",
          sum(GENE_ROWS[g]["inheritance"] == "AR" for g in ACMG) == 97
          and all(GENE_ROWS[g]["inheritance"] == "XL" for g in ACMG_XL))
    check("the ACOG genes are marked", {g for g, r in GENE_ROWS.items() if "ACOG" in r["lists"]} == set(ACOG))
    check("every row has a span on chr1 to chrX",
          all(r["chrom"] in cr.PRIMARY[:-1] and 0 < r["start"] <= r["end"] for r in GENE_ROWS.values()))
    par = cr.family_report.read_par(PAR)
    check("a chrX gene outside the PARs is X-linked, inside a PAR autosomal",
          all((r["inheritance"] == "XL") == (r["chrom"] == "chrX" and not cr.family_report.in_par(par, "chrX", r["start"]))
              for r in GENE_ROWS.values()))
    check("more than the panel: the Orphanet genes are there", len(GENE_ROWS) > 3000, len(GENE_ROWS))

    with tempfile.TemporaryDirectory() as t:
        ped = put(f"{t}/pedigree.tsv", "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\tB\nB\tmale\t.\t.\tA\n")
        cx, cy = mid("CYP21A2")
        exons = put(f"{t}/exons.bed", f"{cx}\t{cy - 200}\t{cy + 200}\n{cx}\t{cy + 5000}\t{cy + 5100}\n")
        tx, ty = mid("TNXB", 20000)
        a_hits = put(f"{t}/A.vcf", HEAD + hit("NAGA") + hit("APOE", at=("chr19", 44908684)))
        b_hits = put(f"{t}/B.vcf", HEAD + hit("NAGA", off=7) + hit("UGT1A1") + hit("CYP21A2") + hit("TNXB", off=20000))
        a_sv = put(f"{t}/A_sv.vcf", HEAD + deletion(cx, cy - 1000, cy + 1000)
                   + deletion(tx, ty - 300, ty + 300)               # intron of TNXB: no exon in the list
                   + deletion(cx, cy - 1000, cy + 1000, filt="LowQual"))
        base = ["screen", "--pedigree", ped, "--genes", GENES, "--par-bed", PAR,
                "--clinvar", f"A={a_hits}", "--clinvar", f"B={b_hits}"]
        out = f"{t}/out"
        rc, log = run(*base, "--delly", f"A={a_sv}", "--exons", exons, "--out-dir", out)
        check("screen exits 0", rc == 0, log)
        tab = rows(f"{out}/couple_A_B.tsv")
        res = {r["gene"]: r for r in tab}
        # 2.
        check("NAGA, one variant in each partner: shared, 1 in 4 per pregnancy",
              res.get("NAGA", {}).get("result") == "shared" and res["NAGA"]["offspring_risk"] == "1 in 4 per pregnancy", res.get("NAGA"))
        check("UGT1A1, a variant in one partner only: one_partner", res.get("UGT1A1", {}).get("result") == "one_partner", res.get("UGT1A1"))
        check("APOE (not in the list) is left out", "APOE" not in res and "APOE" not in cr.read_genes(GENES))
        # 3.
        check("a deletion over a CYP21A2 exon and the other partner's SNV: shared",
              res.get("CYP21A2", {}).get("result") == "shared" and "deletion" in res["CYP21A2"]["findings_a"]
              and "1 exon of the gene (Delly)" in res["CYP21A2"]["findings_a"], res.get("CYP21A2"))
        check("the filtered copy of the deletion is not counted", res["CYP21A2"]["findings_a"].count("deletion") == 1,
              res["CYP21A2"]["findings_a"])
        check("a deletion inside a TNXB intron does not count: TNXB is one partner's only",
              res.get("TNXB", {}).get("result") == "one_partner" and res["TNXB"]["findings_a"] == ".", res.get("TNXB"))
        check("shared genes come first", [r["result"] for r in tab][:2] == ["shared", "shared"], [r["gene"] for r in tab])
        check("the partner tables are the couple table",
              open(f"{out}/A_couple_carrier.tsv").read() == open(f"{out}/couple_A_B.tsv").read()
              == open(f"{out}/B_couple_carrier.tsv").read())
        md = open(f"{out}/couple_A_B.md").read()
        check("the report names NAGA as shared and points to a genetic counsellor",
              "**NAGA**" in md.split("## X-linked")[0] and "genetic counsellor" in md)
        # 7.
        check("B's missing Delly calls, both partners' Parascopy and ExpansionHunter are named",
              "Deletions of B: no Delly calls" in md and "SMN1 copy number of A: no Parascopy" in md
              and "FMR1 and FXN repeats of B: no ExpansionHunter" in md, md)
        check("stdout says the counts", "A and B: 2 shared genes" in log, log)
        # 9.
        out2 = f"{t}/out2"
        run(*base, "--delly", f"A={a_sv}", "--exons", exons, "--out-dir", out2)
        check("the same inputs give the same bytes",
              all(open(f"{out}/{n}").read() == open(f"{out2}/{n}").read() for n in ("couple_A_B.tsv", "couple_A_B.md")))
        # 3. without the exon list
        rc, log = run(*base, "--delly", f"A={a_sv}", "--out-dir", f"{t}/noexons")
        r3 = {r["gene"]: r for r in rows(f"{t}/noexons/couple_A_B.tsv")}
        check("without the exon list the deletion is not used, and the report says so",
              r3["CYP21A2"]["result"] == "one_partner"
              and "without the GENCODE exon list" in open(f"{t}/noexons/couple_A_B.md").read(), r3.get("CYP21A2"))

        # 4. homozygous, and two hets in one partner
        a_hom = put(f"{t}/A_hom.vcf", HEAD + hit("NAGA", gt="1/1"))
        b_two = put(f"{t}/B_two.vcf", HEAD + hit("NAGA", off=7) + hit("NAGA", off=9))
        run("screen", "--pedigree", ped, "--genes", GENES, "--clinvar", f"A={a_hom}", "--clinvar", f"B={b_two}",
            "--out-dir", f"{t}/hom")
        r4 = rows(f"{t}/hom/couple_A_B.tsv")[0]
        check("a homozygous partner: 1 in 2 per pregnancy", r4["offspring_risk"] == "1 in 2 per pregnancy", r4)
        check("two het changes of B: the chance can be higher", "B has two changes in the gene" in r4["explanation"], r4)

        # 5. X-linked and repeats
        a_eh = put(f"{t}/A_eh.vcf", HEAD.replace("\tS\n", "\tA\n") + eh("FMR1", [30, 85]) + eh("FXN", [9, 120]))
        b_eh = put(f"{t}/B_eh.vcf", HEAD + eh("FMR1", [29]) + eh("FXN", [12, 700]))
        b_dmd = put(f"{t}/B_dmd.vcf", HEAD + hit("DMD", gt="1"))
        a_none = put(f"{t}/A_none.vcf", HEAD)
        run("screen", "--pedigree", ped, "--genes", GENES, "--par-bed", PAR, "--clinvar", f"A={a_none}",
            "--clinvar", f"B={b_dmd}", "--eh", f"A={a_eh}", "--eh", f"B={b_eh}", "--out-dir", f"{t}/xl")
        r5 = {r["gene"]: r for r in rows(f"{t}/xl/couple_A_B.tsv")}
        check("FMR1 alleles of 30 and 85 in the female partner: X-linked, sons 1 in 2, premutation note",
              r5.get("FMR1", {}).get("result") == "x_linked" and r5["FMR1"]["offspring_risk"] == "A: sons 1 in 2"
              and "premutation (55 to 200)" in r5["FMR1"]["findings_a"]
              and "grow into a full mutation" in r5["FMR1"]["explanation"], r5.get("FMR1"))
        check("the male partner's FMR1 allele of 29 is no finding", r5.get("FMR1", {}).get("findings_b") == ".", r5.get("FMR1"))
        check("DMD in the male partner: every daughter a carrier",
              r5.get("DMD", {}).get("offspring_risk") == "B: daughters carriers" and "hemizygous" in r5["DMD"]["findings_b"],
              r5.get("DMD"))
        check("FXN expanded in both: shared", r5.get("FXN", {}).get("result") == "shared", r5.get("FXN"))

        # 6. SMN1
        sa, sb, su = put(f"{t}/A_smn.tsv", smn(1)), put(f"{t}/B_smn.tsv", smn(1)), put(f"{t}/U_smn.tsv", smn(1, qual="5"))
        run("screen", "--pedigree", ped, "--genes", GENES, "--clinvar", f"A={a_none}", "--clinvar", f"B={a_none}",
            "--smn", f"A={sa}", "--smn", f"B={sb}", "--out-dir", f"{t}/smn")
        r6 = {r["gene"]: r for r in rows(f"{t}/smn/couple_A_B.tsv")}
        check("SMN1 one copy in each partner: shared", r6.get("SMN1", {}).get("result") == "shared", r6)
        run("screen", "--pedigree", ped, "--genes", GENES, "--clinvar", f"A={a_none}", "--clinvar", f"B={a_none}",
            "--smn", f"A={su}", "--smn", f"B={sb}", "--out-dir", f"{t}/smnu")
        md6 = open(f"{t}/smnu/couple_A_B.md").read()
        check("an SMN1 call of quality 5 is not a finding and is named as not resolved",
              [r["findings_a"] for r in rows(f"{t}/smnu/couple_A_B.tsv") if r["gene"] == "SMN1"] == ["."]
              and "SMN1 copy number of A: SMN1 copy number 1 (quality 5, filter PASS), not resolved" in md6, md6)

        # 3b. read depth: a long Delly deletion counts only when RDCN shows the loss
        genes_mid = {g: mid(g) for g in ("GBA1", "PAH", "CFTR", "HBB", "DMD", "OTC")}
        dexons = put(f"{t}/depth_exons.bed", "".join(
            f"{c}\t{p - 50}\t{p + 50}\n" for c, p in sorted(genes_mid.values())))
        mb = 25000000
        (gc, gp), (pc, pp), (fc, fp) = genes_mid["GBA1"], genes_mid["PAH"], genes_mid["CFTR"]
        (hc, hp), (dc, dp), (oc, op) = genes_mid["HBB"], genes_mid["DMD"], genes_mid["OTC"]
        a_depth = put(f"{t}/A_depth.vcf", HEAD
                      + deldepth(gc, gp - mb, gp + mb, 2)                        # 50 Mb, depth shows 2 copies
                      + deldepth(pc, pp - mb, pp + mb, 1)                        # 50 Mb, depth shows 1 copy
                      + deldepth(fc, fp - 300, fp + 300, 2, precise="PRECISE")  # 600 bp, split reads alone
                      + deldepth(hc, hp - 1000, hp + 1000, "."))                # 2 kb, no depth estimate
        b_depth = put(f"{t}/B_depth.vcf", HEAD
                      + deldepth(dc, dp - 1000000, dp + 1000000, 1, gt="1")      # male chrX: 1 copy is normal
                      + deldepth(oc, op - 2500, op + 2500, 0, gt="1"))           # male chrX: 0 copies, a loss
        empty_hits = put(f"{t}/none_hits.vcf", HEAD)
        rc, log = run("screen", "--pedigree", ped, "--genes", GENES, "--par-bed", PAR,
                      "--clinvar", f"A={empty_hits}", "--clinvar", f"B={empty_hits}",
                      "--delly", f"A={a_depth}", "--delly", f"B={b_depth}", "--exons", dexons,
                      "--out-dir", f"{t}/depth")
        rd = {r["gene"]: r for r in rows(f"{t}/depth/couple_A_B.tsv")} if rc == 0 else {}
        check("read depth: the screen exits 0", rc == 0, log)
        check("a 50 Mb PASS het deletion with RDCN 2 over GBA1 is set aside", "GBA1" not in rd, rd.get("GBA1"))
        check("a 50 Mb PASS het deletion with RDCN 1 over PAH still counts",
              rd.get("PAH", {}).get("result") == "one_partner"
              and "(50000.0 kb) het, 1 exon of the gene (Delly)" in rd["PAH"]["findings_a"], rd.get("PAH"))
        check("a 600 bp precise deletion counts on split reads alone, whatever its RDCN",
              rd.get("CFTR", {}).get("result") == "one_partner" and "(600 bp) het" in rd["CFTR"]["findings_a"],
              rd.get("CFTR"))
        check("a 2 kb deletion without a depth estimate (RDCN .) is set aside", "HBB" not in rd, rd.get("HBB"))
        check("a male's chrX deletion with RDCN 1 (his normal copy) is set aside", "DMD" not in rd, rd.get("DMD"))
        check("a male's chrX deletion with RDCN 0 counts, hemizygous",
              rd.get("OTC", {}).get("result") == "x_linked" and "hemizygous" in rd["OTC"]["findings_b"], rd.get("OTC"))
        mdd = open(f"{t}/depth/couple_A_B.md").read() if rc == 0 else ""
        check("the report says how many long deletions read depth set aside",
              "| A | 0 / 0 | 2 over exons of listed genes; 2 longer than 1000 bp set aside, read depth does not confirm the loss |" in mdd
              and "| B | 0 / 0 | 1 over exons of listed genes; 1 longer than 1000 bp set aside, read depth does not confirm the loss |" in mdd,
              mdd)

        # 8. red controls
        empty = put(f"{t}/empty.tsv", "# emptied\n" + "\t".join(cr.GENE_COLUMNS) + "\n")
        rc, log = run("screen", "--pedigree", ped, "--genes", empty, "--clinvar", f"A={a_hits}",
                      "--clinvar", f"B={b_hits}", "--out-dir", f"{t}/red")
        check("red: an emptied gene list stops the screen (exit 2) and writes no table",
              rc == 2 and "holds no gene" in log and not os.path.exists(f"{t}/red/couple_A_B.tsv"), (rc, log))
        lone = put(f"{t}/lone.tsv", "sample\tsex\tfather\tmother\tpartner\nA\tfemale\t.\t.\t.\nB\tmale\t.\t.\t.\n")
        rc, log = run("screen", "--pedigree", lone, "--genes", GENES, "--clinvar", f"A={a_hits}", "--out-dir", f"{t}/lone")
        check("a pedigree without partners exits 1", rc == 1 and "no two people as partners" in log, (rc, log))
        rc, log = run("screen", "--pedigree", ped, "--genes", GENES, "--clinvar", f"A={a_hits}", "--out-dir", f"{t}/half")
        check("a couple where one partner has no ClinVar screen is not screened (exit 1)",
              rc == 1 and "B has no ClinVar screen" in log and not os.path.exists(f"{t}/half/couple_A_B.tsv"), (rc, log))
        rc, log = run("screen", "--pedigree", ped, "--genes", GENES, "--clinvar", f"C={a_hits}", "--out-dir", f"{t}/c")
        check("an input for a sample outside the pedigree stops it (exit 2)", rc == 2 and "'C' is not in the pedigree" in log, log)

        # 10. genes
        gtf = put(f"{t}/g.gtf", "##description: evidence-based annotation of the human genome (GRCh38), version 99\n" + "".join(
            f'{c}\tHAVANA\tgene\t{s}\t{e}\t.\t+\t.\tgene_id "X"; gene_name "{g}";\n'
            for g, c, s, e in (("NAGA", "chr22", 100, 900), ("AAA1", "chr3", 10, 20), ("XXL1", "chrX", 5000000, 5000100),
                               ("PARG1", "chrX", 20000, 30000), ("PARG1", "chrY", 20000, 30000), ("YYY1", "chrY", 9000000, 9000100),
                               ("C9orf72", "chr9", 1, 2))))
        p9 = put(f"{t}/p9.xml", "<JDBOR><DisorderList>" + "".join(
            f"<Disorder><OrphaCode>{c}</OrphaCode><TypeOfInheritanceList><TypeOfInheritance><Name>{m}</Name>"
            "</TypeOfInheritance></TypeOfInheritanceList></Disorder>"
            for c, m in (("1", "Autosomal recessive"), ("2", "X-linked recessive"), ("3", "Autosomal dominant"))) + "</DisorderList></JDBOR>")

        def dis(code, name, sym, kind="Disease-causing germline mutation(s) in"):
            return (f"<Disorder><OrphaCode>{code}</OrphaCode><Name>{name}</Name><DisorderGeneAssociationList>"
                    f"<DisorderGeneAssociation><Gene><Symbol>{sym}</Symbol></Gene><DisorderGeneAssociationType><Name>{kind}"
                    "</Name></DisorderGeneAssociationType><DisorderGeneAssociationStatus><Name>Assessed</Name>"
                    "</DisorderGeneAssociationStatus></DisorderGeneAssociation></DisorderGeneAssociationList></Disorder>")
        p6 = put(f"{t}/p6.xml", '<JDBOR date="2026-01-02 03:04:05"><DisorderList>'
                 + dis("1", "Recessive one", "AAA1") + dis("2", "X one", "XXL1") + dis("1", "Pseudo", "PARG1")
                 + dis("1", "On Y", "YYY1") + dis("3", "Dominant", "DOM1") + dis("1", "Candidate", "AAA1", "Candidate gene tested in")
                 + dis("1", "Upper case orf", "C9ORF72") + "</DisorderList></JDBOR>")
        seed = put(f"{t}/seed.tsv", "\t".join(cr.GENE_COLUMNS) + "\nNAGA\t.\t0\t0\tAR\tACMG\tSchindler disease\n"
                   "OLD1\tchr1\t1\t2\tAR\tOrphanet\tdropped\n")
        rc, log = run("genes", "--table", seed, "--orphanet-genes", p6, "--orphanet-inheritance", p9, "--gtf", gtf,
                      "--par-bed", PAR, "--out", f"{t}/new.tsv")
        new = cr.read_genes(f"{t}/new.tsv") if rc == 0 else {}
        check("genes: the ACMG row is kept with its GENCODE span", new.get("NAGA", {}).get("start") == 100
              and new["NAGA"]["condition"] == "Schindler disease", (rc, log, new.get("NAGA")))
        check("genes: an autosomal Orphanet gene is AR, a chrX one XL, a PAR one AR, chrY and dominant left out",
              new.get("AAA1", {}).get("inheritance") == "AR" and new.get("XXL1", {}).get("inheritance") == "XL"
              and new.get("PARG1", {}).get("inheritance") == "AR" and new.get("PARG1", {}).get("chrom") == "chrX"
              and "YYY1" not in new and "DOM1" not in new, sorted(new))
        check("genes: an old Orphanet row is made again, not kept; an upper-case orf symbol is matched",
              "OLD1" not in new and "C9orf72" in new and new["AAA1"]["condition"] == "Recessive one", sorted(new))
        check("genes: the header names the releases",
              "release 99" in open(f"{t}/new.tsv").read() and "2026-01-02" in open(f"{t}/new.tsv").read())
        bad = put(f"{t}/bad.tsv", "\t".join(cr.GENE_COLUMNS) + "\nNOSUCH\t.\t0\t0\tAR\tACMG\tx\n")
        rc, log = run("genes", "--table", bad, "--orphanet-genes", p6, "--orphanet-inheritance", p9, "--gtf", gtf,
                      "--par-bed", PAR, "--out", f"{t}/bad_out.tsv")
        check("genes: a panel gene GENCODE does not have stops it", rc == 2 and "NOSUCH" in log, (rc, log))

    print("\nRESULT:", "ALL PASS" if not FAILS else "%d FAILED" % len(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
