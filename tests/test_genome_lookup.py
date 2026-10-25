#!/usr/bin/env python3
"""bin/genome_lookup.py (scripts/pgp-lookup.sh) on a synthetic gVCF, annotated VCF and summary.

Checks:
  1. a query is read as an rsID, a position (chr7:1,234 / 7:1234 / chrMT:5)
     or a gene symbol; anything else is refused;
  2. resolve: an rsID through the consumer site table (with its GRCh38 base),
     through assets/traits.tsv, and, when neither has it, through the ID and
     Existing_variation of the annotated VCF (rs1001 does not match rs10011);
     a gene through GENCODE's gene lines (case-insensitive as a fallback);
     an unknown rsID or gene exits 3;
  3. report, a site: a position inside a reference block is ref/ref (with the
     table's base), a position no record covers is a no call (never the
     reference), a no-call block is a no call, a variant gives its genotype
     with depth and quality and the annotated VCF's gene, consequence,
     ClinVar and gnomAD, a position inside a deletion names the deletion, and
     a filtered call says so;
  4. report, a gene: every annotated variant of the gene (not one of another
     gene at the same span), the share of bases with a call, and the reports
     that mention the gene (a section that did not run is not read); without
     an annotated VCF the gVCF's own alt calls are listed;
  5. the command line resolves and reports end to end;
  6. a region given to Vcf keeps only the records that overlap it (a
     reference block by its END, a deletion by its span), so a whole file
     read without an index stays in bounded memory.

Run: python3 tests/test_genome_lookup.py
"""
import gzip
import io
import json
import os
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "bin"))
import genome_lookup as gl  # noqa: E402

SCRIPT = os.path.join(REPO, "bin", "genome_lookup.py")
TRAITS = os.path.join(REPO, "assets", "traits.tsv")
FAILS = 0


def check(desc, ok, detail=""):
    global FAILS
    print(f"[{'PASS' if ok else 'FAIL'}] {desc}{'' if ok else ' -- ' + str(detail)[:900]}")
    if not ok:
        FAILS += 1


HEAD = "\n".join([
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr7,length=1000000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH1",
])

GVCF = HEAD + "\n" + "\n".join([
    "chr7\t100\t.\tA\t<*>\t0\t.\tEND=200\tGT:GQ:MIN_DP\t0/0:50:31",
    "chr7\t250\t.\tC\tT,<*>\t42.5\tPASS\t.\tGT:GQ:DP\t0/1:40:33",
    "chr7\t251\t.\tG\t<*>\t0\t.\tEND=299\tGT:GQ:MIN_DP\t0/0:45:30",
    # 300-399: no record at all
    "chr7\t400\t.\tT\t<*>\t0\t.\tEND=410\tGT:GQ:MIN_DP\t./.:0:0",
    "chr7\t500\t.\tACGT\tA,<*>\t38\tPASS\t.\tGT:GQ:DP\t0/1:35:29",
    "chr7\t600\t.\tG\tA,<*>\t3.1\tRefCall\t.\tGT:GQ:DP\t0/1:5:4",
    "chr7\t650\t.\tT\tG,<*>\t50\tPASS\t.\tGT:GQ:DP\t1/1:50:40",
]) + "\n"

CSQ = ("Allele|Consequence|IMPACT|SYMBOL|Gene|HGVSc|HGVSp|Existing_variation|CANONICAL|CLIN_SIG|"
       "gnomADg_AF|ClinVar_CLNSIG|ClinVar_CLNREVSTAT")
ANN_HEAD = HEAD.replace(
    "#CHROM", f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
              f'Format: {CSQ}">\n#CHROM')


def csq(*entries):
    return "CSQ=" + ",".join("|".join(e) for e in entries)


E250_OTHER = ["T", "upstream_gene_variant", "MODIFIER", "GENEZ", "ENSG2", "", "", "rs1001", "YES", "", "0.12", "", ""]
E250_MAIN = ["T", "missense_variant", "MODERATE", "GENEA", "ENSG1", "ENST1:c.10C>T", "ENSP1:p.Arg4Cys",
             "rs1001&COSV1", "YES", "pathogenic", "0.0003", "Pathogenic", "criteria_provided&_multiple_submitters"]
ANN = ANN_HEAD + "\n" + "\n".join([
    f"chr7\t250\t.\tC\tT\t42.5\tPASS\t{csq(E250_OTHER, E250_MAIN)}\tGT:GQ:DP\t0/1:40:33",
    "chr7\t500\t.\tACGT\tA\t38\tPASS\t" + csq(["-", "inframe_deletion", "MODERATE", "GENEA", "ENSG1", "", "",
                                              "rs10011", "YES", "", "0.2", "", ""]) + "\tGT:GQ:DP\t0/1:35:29",
    "chr7\t650\t.\tT\tG\t50\tPASS\t" + csq(["G", "intron_variant", "MODIFIER", "GENEZ", "ENSG2", "", "", "", "YES",
                                           "", "0.5", "", ""]) + "\tGT:GQ:DP\t1/1:50:40",
]) + "\n"

GTF = "\n".join([
    "##description: synthetic",
    'chr7\tHAVANA\tgene\t90\t700\t.\t+\t.\tgene_id "ENSG1"; gene_type "protein_coding"; gene_name "GENEA"; level 2;',
    'chr7\tHAVANA\tgene\t640\t660\t.\t-\t.\tgene_id "ENSG2"; gene_type "lncRNA"; gene_name "GENEZ"; level 2;',
]) + "\n"

SITES = "rsid\tchrom\tpos37\tref37\tpos38\tref38\tstrand\nrs2002\t7\t140\tA\t150\tA\t+\nrs3003\tMT\t73\tA\t73\tA\t+\n"

SUMMARY = {"sections": {
    "pharmcat": {"title": "PharmCAT", "step": "07", "state": "ok",
                 "data": {"genes": [{"gene": "GENEA", "diplotype": "*1/*2", "phenotype": "Intermediate"},
                                    {"gene": "OTHER", "diplotype": "*1/*1"}]}},
    "traits": {"title": "Traits panel", "step": "42", "state": "ok",
               "data": {"rows": [{"rsid": "rs2002", "gene": "GENEB", "genotype": "A/A", "reading": "synthetic"}]}},
    "clinvar": {"title": "ClinVar screen", "step": "06", "state": "missing",
                "data": {"rows": [{"gene": "GENEA", "note": "MUST-NOT-PRINT"}]}},
}}


def run_report(**kw):
    args = ["report"]
    for k, v in kw.items():
        args += [f"--{k}", str(v)]
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = gl.main(args)
    return rc, buf.getvalue()


def main():
    tmp = tempfile.mkdtemp()
    p = lambda n: os.path.join(tmp, n)  # noqa: E731
    with gzip.open(p("g.vcf.gz"), "wt") as f:
        f.write(GVCF)
    with open(p("ann.vcf"), "w") as f:
        f.write(ANN)
    with open(p("genes.gtf"), "w") as f:
        f.write(GTF)
    with gzip.open(p("sites.tsv.gz"), "wt") as f:
        f.write(SITES)
    with open(p("summary.json"), "w") as f:
        json.dump(SUMMARY, f)

    # 1. queries
    check("rsID query", gl.parse_query("RS4988235") == ("rsid", "rs4988235"))
    check("chr:pos query with commas", gl.parse_query("chr7:1,234") == ("site", ("7", 1234)))
    check("pos query without chr", gl.parse_query("7:1234") == ("site", ("7", 1234)))
    check("chrMT is chrM", gl.parse_query("chrMT:5") == ("site", ("MT", 5)) and gl.contig_name("MT") == "chrM")
    check("gene query", gl.parse_query("CFTR") == ("gene", "CFTR"))
    for bad in ("chr7:0", "a b", "rs12;x"):
        try:
            gl.parse_query(bad)
            check(f"refuses {bad!r}", False, "accepted")
        except ValueError:
            check(f"refuses {bad!r}", True)

    def resolve(*args):
        out = io.StringIO()
        with redirect_stdout(out):
            rc = gl.main(["resolve", *args])
        return rc, [line.split("\t") for line in out.getvalue().splitlines()]

    # 2. resolve
    rc, rows = resolve("--query", "rs2002", "--sites", p("sites.tsv.gz"), "--annotated", p("ann.vcf"))
    check("rsID from the site table, with its GRCh38 base", rc == 0 and rows and rows[0][:5] == ["site", "chr7", "150", "150", "A"], rows)
    rc, rows = resolve("--query", "rs4988235", "--table", TRAITS)
    check("rsID from assets/traits.tsv", rc == 0 and rows and rows[0][:5] == ["site", "chr2", "135851076", "135851076", "G"], rows)
    rc, rows = resolve("--query", "rs1001", "--sites", p("sites.tsv.gz"), "--annotated", p("ann.vcf"))
    check("rsID from the annotated VCF's Existing_variation, once (rs10011 is another rsID)",
          rc == 0 and [r[:3] for r in rows] == [["site", "chr7", "250"]], rows)
    rc, rows = resolve("--query", "rs999", "--sites", p("sites.tsv.gz"), "--annotated", p("ann.vcf"))
    check("an unknown rsID exits 3", rc == 3 and not rows, (rc, rows))
    rc, rows = resolve("--query", "GENEA", "--genes", p("genes.gtf"))
    check("gene from GENCODE's gene lines", rc == 0 and rows and rows[0][:4] == ["gene", "chr7", "90", "700"], rows)
    rc, rows = resolve("--query", "genea", "--genes", p("genes.gtf"))
    check("gene name case-insensitive as a fallback", rc == 0 and rows and rows[0][2] == "90", rows)
    rc, rows = resolve("--query", "NOPE1", "--genes", p("genes.gtf"))
    check("an unknown gene exits 3", rc == 3, rc)
    rc, rows = resolve("--query", "GENEA")
    check("a gene without the gene lines exits 3", rc == 3, rc)

    base = dict(gvcf=p("g.vcf.gz"), annotated=p("ann.vcf"), summary=p("summary.json"), chrom="chr7")

    # 3. sites
    rc, out = run_report(query="rs2002", kind="site", start=150, end=150, ref="A", **base)
    check("a position in a reference block is ref/ref with the table's base",
          rc == 0 and "genotype: A/A (homozygous reference, inside the reference block chr7:100-200)" in out
          and "lowest depth in the block 31" in out and "genotype quality 50" in out, out)
    check("a reference site has no annotated variant, said so",
          "no variant of this sample here in the annotated VCF" in out, out)
    check("an rsID the reports mention is listed (traits row)",
          "Traits panel (step 42):" in out and "rsid=rs2002" in out, out)
    rc, out = run_report(query="chr7:150", kind="site", start=150, end=150, **base)
    check("a position in a reference block without a known base is ref/ref",
          "genotype: ref/ref (homozygous reference" in out, out)
    rc, out = run_report(query="chr7:350", kind="site", start=350, end=350, **base)
    check("a position no record covers is a no call, never the reference",
          "genotype: no call (no gVCF record covers this position" in out and "reference" not in out.split("genotype:")[1].split("\n")[0], out)
    rc, out = run_report(query="chr7:405", kind="site", start=405, end=405, **base)
    check("a no-call block is a no call", "genotype: no call (reference block chr7:400-410 has no genotype)" in out, out)
    rc, out = run_report(query="rs1001", kind="site", start=250, end=250, **base)
    check("a variant gives its genotype, depth and quality",
          "genotype: C/T (heterozygous; REF C, ALT T)" in out and "depth 33" in out and "QUAL 42.5" in out, out)
    check("a variant gets the gene and consequence of its most severe entry (missense in GENEA, not GENEZ upstream)",
          "gene: GENEA" in out and "consequence: missense_variant" in out and "gene: GENEZ" not in out, out)
    check("ClinVar significance and review status, and the gnomAD frequency, from the annotated VCF",
          "ClinVar (pipeline's ClinVar file): Pathogenic" in out
          and "ClinVar review status: criteria provided, multiple submitters" in out
          and "gnomAD genomes AF: 0.0003" in out and "clinical significance (VEP cache): pathogenic" in out, out)
    check("the reports of the variant's genes are read (PharmCAT names GENEA)",
          "PharmCAT (step 07):" in out and "diplotype=*1/*2" in out and "MUST-NOT-PRINT" not in out, out)
    rc, out = run_report(query="chr7:502", kind="site", start=502, end=502, **base)
    check("a position inside a deletion names the deletion",
          "lies inside the variant called at chr7:500 ACGT>A, genotype ACGT/A (heterozygous)" in out, out)
    rc, out = run_report(query="chr7:600", kind="site", start=600, end=600, **base)
    check("a filtered call says it did not pass", "filter: RefCall" in out, out)

    # 4. genes
    rc, out = run_report(query="GENEA", kind="gene", start=90, end=700, source="GENCODE 50 (GENEA)", **base)
    check("gene: its variants are listed with consequence, ClinVar and gnomAD",
          "chr7:250 C>T C/T (heterozygous; COSV1,rs1001" in out and "consequence: missense_variant" in out
          and "HGVSp: ENSP1:p.Arg4Cys" in out and "ClinVar (pipeline's ClinVar file): Pathogenic" in out
          and "gnomAD genomes AF: 0.0003" in out and "chr7:500 ACGT>A" in out, out)
    check("gene: a variant of another gene in the span is not listed", "chr7:650" not in out, out)
    # 90-99 no record, 100-200 call, 250 call, 251-299 call, 300-399 none, 400-410 no-call,
    # 500-503 call, 600 call, 650 call: 101 + 1 + 49 + 4 + 1 + 1 = 157 of 611
    check("gene: the share of bases with a call", "157 of 611 bases have a genotype call (25.7%)" in out, out)
    check("gene: the reports that mention it", "PharmCAT (step 07):" in out and "phenotype=Intermediate" in out
          and "MUST-NOT-PRINT" not in out and "OTHER" not in out, out)
    nob = dict(base, annotated=None)
    nob.pop("annotated")
    rc, out = run_report(query="GENEA", kind="gene", start=90, end=700, **nob)
    check("gene without an annotated VCF: the gVCF's alt calls are listed, a RefCall too, blocks not",
          "chr7:250 C>T C/T" in out and "chr7:500 ACGT>A" in out and "chr7:650 T>G G/G (homozygous" in out
          and "filter RefCall" in out and "chr7:100" not in out, out)
    nos = dict(nob)
    nos.pop("summary")
    rc, out = run_report(query="GENEA", kind="gene", start=90, end=700, **nos)
    check("without summary.json the answer says how to get one", "no summary.json for this sample" in out, out)

    # 6. a region keeps only the records that overlap it, so an unindexed
    # whole genome read whole stays in bounded memory
    def kept(path, chrom, start, end):
        return [r.pos for r in gl.Vcf(path, chrom, start, end).records]
    check("region: a reference block is kept by its END", kept(p("g.vcf.gz"), "chr7", 150, 150) == [100],
          kept(p("g.vcf.gz"), "chr7", 150, 150))
    check("region: a deletion is kept by the bases it spans", kept(p("g.vcf.gz"), "7", 503, 503) == [500],
          kept(p("g.vcf.gz"), "7", 503, 503))
    check("region: a gap keeps nothing", kept(p("g.vcf.gz"), "chr7", 300, 399) == [],
          kept(p("g.vcf.gz"), "chr7", 300, 399))
    check("region: another contig keeps nothing", kept(p("g.vcf.gz"), "chr8", 1, 1000) == [],
          kept(p("g.vcf.gz"), "chr8", 1, 1000))
    check("region: no region keeps every record", len(gl.Vcf(p("g.vcf.gz")).records) == 7)

    # 5. command line
    r = subprocess.run([sys.executable, SCRIPT, "resolve", "--query", "rs2002", "--sites", p("sites.tsv.gz")],
                       capture_output=True, text=True)
    check("command line: resolve", r.returncode == 0 and r.stdout.startswith("site\tchr7\t150\t150\tA\t"), r)
    r = subprocess.run([sys.executable, SCRIPT, "report", "--query", "chr7:250", "--kind", "site", "--chrom", "chr7",
                        "--start", "250", "--end", "250", "--gvcf", p("g.vcf.gz"), "--annotated", p("ann.vcf")],
                       capture_output=True, text=True)
    check("command line: report", r.returncode == 0 and "genotype: C/T" in r.stdout and "no summary.json" in r.stdout, r)
    r = subprocess.run([sys.executable, SCRIPT, "resolve", "--query", "a b"], capture_output=True, text=True)
    check("command line: a bad query exits 2", r.returncode == 2 and "ERROR" in r.stderr, r)

    print(f"\n{'FAILED: ' + str(FAILS) if FAILS else 'All checks passed.'}")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
