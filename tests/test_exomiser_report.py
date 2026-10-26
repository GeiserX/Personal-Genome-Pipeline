#!/usr/bin/env python3
"""bin/exomiser_report.py: the symptom-driven ranking of step 48 (Exomiser).

Every input is synthetic: a few-term OBO file, one- and four-sample VCF
headers, a pedigree, and Exomiser variants TSVs in the column layout of
Exomiser 15's TsvVariantResultsWriter. The cases:

  1. prepare: the terms are checked against the ontology (an alternative id is
     read as its term, a repeat is dropped), the labels are written, and the
     analysis holds the VCF, the terms, the exome preset and the output options
  2. prepare stops on a term that is malformed, unknown, obsolete (naming its
     replacement), on no term at all, on a VCF of several samples without
     --pedigree, and on --affected without --pedigree
  3. the family file: the proband affected, --affected and --unaffected as
     named, everyone else unknown, parents only when sequenced; a stray name,
     a proband named unaffected, and a VCF sample the pedigree lacks all stop
  4. report: the top ranks only, the contributing variants only, both
     variants of a compound heterozygote, the scores and ACMG class, the
     words with the note; a header-only file gives no rank and says so; a
     missing column and a short line stop it
  5. the summary section and both reports show the ranking (collect_summary.py
     and render_report.py); the module's -stub table has the same columns and
     reads as no gene

Run: python3 tests/test_exomiser_report.py
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "exomiser_report.py")
sys.path.insert(0, os.path.join(REPO, "bin"))

FAILS = []

OBO = """format-version: 1.2

[Term]
id: HP:0001156
name: Brachydactyly
alt_id: HP:0009803

[Term]
id: HP:0001363
name: Craniosynostosis

[Term]
id: HP:0011304
name: Broad thumb

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0001250
name: obsolete Seizure
is_obsolete: true
replaced_by: HP:0012345

[Typedef]
id: part_of
name: part of
"""

HEADER = ["#RANK", "ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "MOI", "P-VALUE", "EXOMISER_GENE_COMBINED_SCORE",
          "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_VARIANT_SCORE", "EXOMISER_VARIANT_SCORE",
          "CONTRIBUTING_VARIANT", "WHITELIST_VARIANT", "VCF_ID", "RS_ID", "CONTIG", "START", "END", "REF", "ALT",
          "CHANGE_LENGTH", "QUAL", "FILTER", "GENOTYPE", "FUNCTIONAL_CLASS", "HGVS", "EXOMISER_ACMG_CLASSIFICATION",
          "EXOMISER_ACMG_EVIDENCE", "EXOMISER_ACMG_DISEASE_ID", "EXOMISER_ACMG_DISEASE_NAME", "CLINVAR_VARIATION_ID",
          "CLINVAR_PRIMARY_INTERPRETATION", "CLINVAR_STAR_RATING", "GENE_CONSTRAINT_LOEUF",
          "GENE_CONSTRAINT_LOEUF_LOWER", "GENE_CONSTRAINT_LOEUF_UPPER", "MAX_FREQ_SOURCE", "MAX_FREQ", "ALL_FREQ",
          "MAX_PATH_SOURCE", "MAX_PATH", "ALL_PATH"]


def variant(rank, gene, moi, contig, pos, ref, alt, gt, effect, hgvs, acmg, evidence="", disease_id="", disease="",
            clinvar="NOT_PROVIDED", stars="0", contributing="1", score="0.9000", freq_src="", freq=""):
    v = {k: "" for k in HEADER}
    v.update({"#RANK": str(rank), "ID": "%s-%s-%s-%s_%s" % (contig, pos, ref, alt, moi), "GENE_SYMBOL": gene,
              "ENTREZ_GENE_ID": "1", "MOI": moi, "P-VALUE": "0.0010", "EXOMISER_GENE_COMBINED_SCORE": score,
              "EXOMISER_GENE_PHENO_SCORE": "0.8000", "EXOMISER_GENE_VARIANT_SCORE": "0.9500",
              "EXOMISER_VARIANT_SCORE": "0.9500", "CONTRIBUTING_VARIANT": contributing, "WHITELIST_VARIANT": "0",
              "CONTIG": contig, "START": str(pos), "END": str(pos), "REF": ref, "ALT": alt, "CHANGE_LENGTH": "0",
              "QUAL": "50.0000", "FILTER": "PASS", "GENOTYPE": gt, "FUNCTIONAL_CLASS": effect, "HGVS": hgvs,
              "EXOMISER_ACMG_CLASSIFICATION": acmg, "EXOMISER_ACMG_EVIDENCE": evidence,
              "EXOMISER_ACMG_DISEASE_ID": disease_id, "EXOMISER_ACMG_DISEASE_NAME": disease,
              "CLINVAR_PRIMARY_INTERPRETATION": clinvar, "CLINVAR_STAR_RATING": stars,
              "MAX_FREQ_SOURCE": freq_src, "MAX_FREQ": freq})
    return v


ROWS = [
    variant(1, "GENEA", "AD", "chr10", 123256215, "T", "G", "0/1", "missense_variant",
            "GENEA:ENST1:c.1688A>C:p.(Glu563Ala)", "LIKELY_PATHOGENIC", "PM2,PP3_Strong,PP4", "OMIM:100001",
            "Synthetic syndrome one", "LIKELY_PATHOGENIC", "1", score="0.9981"),
    variant(2, "GENEB", "AR",
            "chr6", 1000, "G", "A", "0/1", "splice_donor_variant", "GENEB:ENST2:c.10+1G>A:p.?",
            "UNCERTAIN_SIGNIFICANCE", "PVS1_Strong", freq_src="GNOMAD_E_SAS", freq="0.0032"),
    variant(2, "GENEB", "AR", "chr6", 2000, "C", "T", "1|0", "stop_gained", "GENEB:ENST2:c.300C>T:p.(Gln100*)",
            "PATHOGENIC", "PVS1,PM2"),
    variant(2, "GENEB", "AR",
            "chr6", 3000, "A", "C", "0/1", "missense_variant", "GENEB:ENST2:c.400A>C:p.(Lys134Thr)",
            "NOT_AVAILABLE", contributing="0"),
    variant(3, "GENEC", "XR", "chrX", 5000, "A", "G", "1", "missense_variant", "GENEC:ENST3:c.5A>G:p.(Asp2Gly)",
            "UNCERTAIN_SIGNIFICANCE"),
    variant(4, "GENED", "AD", "chr1", 7000, "A", "T", "0/1", "missense_variant", "GENED:ENST4:c.7A>T:p.(Asn3Ile)",
            "VUS"),
]


def check(name, ok, detail=""):
    print(("[PASS] " if ok else "[FAIL] ") + name + ("" if ok else ": " + detail))
    if not ok:
        FAILS.append(name)


def run(*args):
    return subprocess.run([sys.executable, SCRIPT] + list(args), capture_output=True, text=True)


def write(path, text):
    with open(path, "w") as f:
        f.write(text)
    return path


def vcf(path, samples):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples) + "\n")
        f.write("\t".join(["chr1", "100", ".", "A", "G", "50", "PASS", ".", "GT"] + ["0/1"] * len(samples)) + "\n")
    return path


def tsv(path, rows, header=HEADER):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r.get(k, "") for k in header) + "\n")
    return path


def read_rows(path):
    with open(path) as f:
        lines = f.read().rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def prepare(d, hpo, samples=("S1",), extra=(), obo=True):
    v = vcf(os.path.join(d, "in.vcf.gz"), list(samples))
    args = ["prepare", "--hpo", hpo, "--sample", "S1", "--vcf", v, "--assembly", "hg38", "--out-dir", "results",
            "--analysis", os.path.join(d, "a.yml"), "--terms", os.path.join(d, "terms.tsv")]
    if obo:
        args += ["--obo", write(os.path.join(d, "hp.obo"), OBO)]
    return run(*(args + list(extra)))


def case_prepare():
    with tempfile.TemporaryDirectory() as d:
        r = prepare(d, "HP:0001156, HP:0001363,HP:0009803 HP:0011304")
        check("prepare: valid terms pass", r.returncode == 0, r.stderr)
        a = json.load(open(os.path.join(d, "a.yml")))
        an = a["analysis"]
        check("prepare: the alternative id is its term and the repeat is dropped",
              an["hpoIds"] == ["HP:0001156", "HP:0001363", "HP:0011304"], str(an["hpoIds"]))
        check("prepare: the VCF and the assembly", an["vcf"].endswith("in.vcf.gz") and an["genomeAssembly"] == "hg38")
        steps = [list(x)[0] for x in an["steps"]]
        check("prepare: the exome preset's steps, prioritisers last",
              steps[-2:] == ["omimPrioritiser", "hiPhivePrioritiser"] and "variantEffectFilter" in steps
              and "frequencyFilter" in steps, str(steps))
        check("prepare: no pedigree for one sample", "ped" not in an and "proband" not in an)
        oo = a["outputOptions"]
        check("prepare: output files named after the sample, TSV and HTML",
              oo["outputFileName"] == "S1_exomiser" and oo["outputDirectory"] == "results"
              and set(oo["outputFormats"]) == {"TSV_VARIANT", "TSV_GENE", "HTML"}, str(oo))
        terms = read_rows(os.path.join(d, "terms.tsv"))
        check("prepare: the terms table has the labels",
              [(t["hpo_id"], t["label"]) for t in terms] == [("HP:0001156", "Brachydactyly"),
                                                            ("HP:0001363", "Craniosynostosis"),
                                                            ("HP:0011304", "Broad thumb")], str(terms))
        r = prepare(d, "HP:0001156", obo=False)
        check("prepare: without the ontology a well-formed term passes, with no label",
              r.returncode == 0 and read_rows(os.path.join(d, "terms.tsv"))[0]["label"] == "", r.stderr)


def case_prepare_stops():
    with tempfile.TemporaryDirectory() as d:
        for hpo, want, name in [("HP:12345", "not an HPO id", "a malformed term"),
                                ("Seizure", "not an HPO id", "a name instead of an id"),
                                ("HP:0009999", "not a term of the HPO release", "an unknown term"),
                                ("HP:0001250", "obsolete", "an obsolete term"),
                                (" , ", "no HPO term given", "no term")]:
            r = prepare(d, hpo)
            check("prepare stops on %s" % name, r.returncode == 1 and want in r.stderr, r.stderr)
        r = prepare(d, "HP:0001250")
        check("prepare names the replacement of an obsolete term", "use HP:0012345" in r.stderr, r.stderr)
        r = prepare(d, "HP:0001156", samples=("S1", "S2"))
        check("prepare stops on a VCF of two samples without --pedigree",
              r.returncode == 1 and "has 2 samples" in r.stderr, r.stderr)
        r = prepare(d, "HP:0001156", extra=["--affected", "S2"])
        check("prepare stops on --affected without --pedigree", r.returncode == 1 and "need --pedigree" in r.stderr,
              r.stderr)


PEDIGREE = """sample\tsex\tfather\tmother\tpartner
S1\tfemale\tDAD\tMUM\t.
SIB\tmale\tDAD\tMUM\t.
DAD\tmale\t.\t.\tMUM
MUM\tfemale\t.\t.\tDAD
"""


def case_family():
    with tempfile.TemporaryDirectory() as d:
        ped = write(os.path.join(d, "pedigree.tsv"), PEDIGREE)
        out = os.path.join(d, "fam.ped")
        fam = ["--pedigree", ped, "--ped", out]
        r = prepare(d, "HP:0001156", samples=("S1", "SIB", "DAD"), extra=fam + ["--unaffected", "DAD",
                                                                                 "--affected", "SIB"])
        check("family: prepare passes", r.returncode == 0, r.stderr)
        lines = [ln.split("\t") for ln in open(out).read().splitlines()]
        got = {x[1]: x for x in lines}
        check("family: one PED line per member of the VCF", sorted(got) == ["DAD", "S1", "SIB"], str(lines))
        check("family: the proband is affected, with its sequenced father and its unsequenced mother as 0",
              got["S1"] == ["family", "S1", "DAD", "0", "2", "2"], str(got["S1"]))
        check("family: a member named affected", got["SIB"][5] == "2" and got["SIB"][4] == "1", str(got["SIB"]))
        check("family: a member named unaffected", got["DAD"][5] == "1", str(got["DAD"]))
        an = json.load(open(os.path.join(d, "a.yml")))["analysis"]
        check("family: the analysis names the PED file and the proband", an.get("ped") == out and an.get("proband") == "S1")
        r = prepare(d, "HP:0001156", samples=("S1", "SIB", "DAD"), extra=fam)
        got = {x[1]: x for x in (ln.split("\t") for ln in open(out).read().splitlines())}
        check("family: a member not named is unknown (0)", got["SIB"][5] == "0" and got["DAD"][5] == "0", str(got))
        for extra, want, name in [(["--affected", "MUM"], "not a member of the family VCF", "a name not in the VCF"),
                                  (["--unaffected", "S1"], "cannot be unaffected", "the proband named unaffected"),
                                  (["--affected", "SIB", "--unaffected", "SIB"], "both affected and unaffected",
                                   "a name in both lists")]:
            r = prepare(d, "HP:0001156", samples=("S1", "SIB", "DAD"), extra=fam + extra)
            check("family stops on %s" % name, r.returncode == 1 and want in r.stderr, r.stderr)
        r = prepare(d, "HP:0001156", samples=("S1", "STRANGER"), extra=fam)
        check("family stops on a VCF sample the pedigree does not name",
              r.returncode == 1 and "STRANGER" in r.stderr, r.stderr)


def report(d, rows, top="3", header=HEADER, raw=None):
    v = os.path.join(d, "v.tsv")
    if raw is not None:
        write(v, raw)
    else:
        tsv(v, rows, header)
    t = write(os.path.join(d, "terms.tsv"), "hpo_id\tlabel\tgiven_as\nHP:0001156\tBrachydactyly\tHP:0001156\n")
    r = run("report", "--sample", "S1", "--variants", v, "--terms", t, "--out", os.path.join(d, "top.tsv"),
            "--text", os.path.join(d, "words.txt"), "--top", top)
    return r


def case_report():
    with tempfile.TemporaryDirectory() as d:
        r = report(d, ROWS)
        check("report passes", r.returncode == 0, r.stderr)
        top = read_rows(os.path.join(d, "top.tsv"))
        check("report: the top 3 ranks only, contributing variants only",
              [(x["rank"], x["gene"], x["variant"]) for x in top] == [
                  ("1", "GENEA", "chr10:123256215 T>G"), ("2", "GENEB", "chr6:1000 G>A"),
                  ("2", "GENEB", "chr6:2000 C>T"), ("3", "GENEC", "chrX:5000 A>G")], str(top))
        a = top[0]
        check("report: scores, inheritance and ACMG of the first rank",
              (a["exomiser_score"], a["phenotype_score"], a["variant_score"], a["moi"], a["acmg"], a["acmg_evidence"])
              == ("0.998", "0.800", "0.950", "AD", "LIKELY_PATHOGENIC", "PM2,PP3_Strong,PP4"), str(a))
        check("report: disease and ClinVar", a["disease"] == "Synthetic syndrome one (OMIM:100001)"
              and a["clinvar"] == "LIKELY_PATHOGENIC (1 star)", str(a))
        check("report: a frequency is named with its source, an unseen variant says so",
              top[1]["max_frequency"] == "GNOMAD_E_SAS 0.0032%" and a["max_frequency"] == "not seen", str(top[1]))
        words = open(os.path.join(d, "words.txt")).read()
        check("words: the terms searched", "Symptoms searched: Brachydactyly (HP:0001156)." in words, words)
        check("words: a rank in a sentence",
              "1. GENEA (autosomal dominant): Exomiser score 0.998 (symptom match 0.800, variant 0.950)." in words
              and "chr10:123256215 T>G heterozygous, missense_variant" in words
              and "ACMG likely pathogenic (PM2,PP3_Strong,PP4)" in words, words)
        check("words: both variants of the compound heterozygote, and a hemizygous one",
              "chr6:1000 G>A heterozygous" in words and "chr6:2000 C>T heterozygous" in words
              and "chrX:5000 A>G hemizygous" in words and "X-linked recessive" in words, words)
        check("words: the note says it is not a diagnosis", "not a diagnosis" in words and "GENED" not in words, words)
        r = report(d, [])
        check("report: a header-only file gives no rank and says so",
              r.returncode == 0 and read_rows(os.path.join(d, "top.tsv")) == []
              and "No gene passed" in open(os.path.join(d, "words.txt")).read(), r.stderr)
        r = report(d, ROWS, header=[h for h in HEADER if h != "EXOMISER_ACMG_CLASSIFICATION"])
        check("report stops on a missing column", r.returncode == 1 and "EXOMISER_ACMG_CLASSIFICATION" in r.stderr,
              r.stderr)
        raw = "\t".join(HEADER) + "\n" + "\t".join(["1", "x", "GENEA"]) + "\n"
        r = report(d, None, raw=raw)
        check("report stops on a short line", r.returncode == 1 and "one field per column" in r.stderr, r.stderr)
        r = report(d, None, raw="")
        check("report stops on an empty file", r.returncode == 1 and "not an Exomiser variants TSV" in r.stderr,
              r.stderr)


def case_reports():
    import collect_summary as cs
    import render_report as rr
    with tempfile.TemporaryDirectory() as d:
        sd = os.path.join(d, "S1", "exomiser")
        os.makedirs(sd)
        r = report(sd, ROWS)
        os.rename(os.path.join(sd, "top.tsv"), os.path.join(sd, "S1_exomiser_top.tsv"))
        os.rename(os.path.join(sd, "terms.tsv"), os.path.join(sd, "S1_exomiser_terms.tsv"))
        p, data = cs.sec_exomiser(os.path.join(d, "S1"), "S1")
        check("summary: the section reads the table", p is not None and r.returncode == 0)
        check("summary: ranks with their variants",
              [(g["rank"], g["gene"], len(g["variants"])) for g in data["genes"]]
              == [("1", "GENEA", 1), ("2", "GENEB", 2), ("3", "GENEC", 1)], str(data["genes"]))
        check("summary: the terms", data["terms"] == [{"id": "HP:0001156", "label": "Brachydactyly"}], str(data["terms"]))
        check("summary: the section is listed under step 48",
              ("exomiser", "Symptom-driven ranking (Exomiser)", "48", cs.sec_exomiser) in cs.SECTIONS
              and "48" in cs.RUN_ALL_STEPS)
        text = "\n".join(rr.exomiser_text(data, "S1"))
        check("text report: rank, gene and variant",
              "1. GENEA (AD) score 0.998" in text and "chr6:2000 C>T" in text and "not a diagnosis" in text, text)
        html = "\n".join(rr.exomiser_html(data, "S1"))
        check("HTML report: a table row per variant, the terms and the note",
              html.count("<tr>") == 1 + 4 and "Brachydactyly (HP:0001156)" in html and "not a diagnosis" in html,
              html)
        bad = os.path.join(sd, "S1_exomiser_top.tsv")
        write(bad, "rank\tgene\n1\tGENEA\n")
        try:
            cs.sec_exomiser(os.path.join(d, "S1"), "S1")
            check("summary: a table without its columns is unreadable", False, "no error")
        except cs.Unreadable:
            check("summary: a table without its columns is unreadable", True)
        # The -stub table of EXOMISER_REPORT reads as "no gene passed", not unreadable.
        import re
        import exomiser_report as er
        nf = open(os.path.join(REPO, "modules", "local", "exomiser", "main.nf")).read()
        m = re.search(r"printf '([^']*)' > \$\{meta\.id\}_exomiser_top\.tsv", nf.split("stub:")[-1])
        head = m.group(1).replace("\\\\t", "\t").replace("\\\\n", "\n") if m else ""
        check("stub: the top table carries the report's columns", head == "\t".join(er.OUT_HEADER) + "\n", repr(head))
        write(bad, head)
        try:
            _, data = cs.sec_exomiser(os.path.join(d, "S1"), "S1")
            check("stub: the summary reads it as no gene", data.get("genes") == [], str(data))
        except cs.Unreadable as e:
            check("stub: the summary reads it as no gene", False, str(e))


if __name__ == "__main__":
    case_prepare()
    case_prepare_stops()
    case_family()
    case_report()
    case_reports()
    if FAILS:
        print("FAILED: %d check(s)" % len(FAILS))
        sys.exit(1)
    print("All checks passed.")
