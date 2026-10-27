#!/usr/bin/env python3
"""exomiser_report.py: the symptom-driven ranking of step 48 (Exomiser).

  exomiser_report.py prepare --hpo HP:0001156,HP:0001363 --sample S --vcf S.vcf.gz
                             --assembly hg38 --out-dir DIR --analysis S_analysis.yml
                             --terms S_exomiser_terms.tsv [--obo hp.obo]
                             [--pedigree pedigree.tsv --ped S_exomiser.ped
                              --affected a,b --unaffected c,d]
  exomiser_report.py report  --sample S --variants S_exomiser.variants.tsv
                             --terms S_exomiser_terms.tsv
                             --out S_exomiser_top.tsv --text S_exomiser.txt [--top 10]

prepare checks the HPO terms (HP: and seven digits; with --obo, each must be a
current term of the ontology Exomiser's phenotype data holds: an alternative id
is read as its term, an obsolete one stops the step with its replacement) and
writes the Exomiser analysis file: the VCF, the terms, Exomiser's exome preset
(coding and splice variants, gnomAD frequencies, the hiPHIVE and OMIM
prioritisers) and the output options (TSV_VARIANT, TSV_GENE and HTML into
--out-dir, files named <sample>_exomiser.*). The file is JSON, which YAML reads.

With --pedigree (step 38's pedigree.tsv) the VCF is the family's joint call:
prepare writes --ped, a PED file of the members the VCF holds, the parents a
member has in the VCF, the sex, and the affected status Exomiser filters on:
the sample (the proband) is affected, the members named in --affected are too,
those in --unaffected are not, and every other member is unknown (0), which
Exomiser does not filter on. A name that is not a member of the VCF stops it.

report reads Exomiser's variants TSV and writes the top --top ranks (a gene
under one mode of inheritance is one rank) with each variant that contributes
to the gene's score: one row per variant, columns OUT_HEADER. The text file
says the same in words, with what the ranking is and is not.

Used by scripts/48-exomiser.sh and the EXOMISER_PREPARE and EXOMISER_REPORT
processes (bin/ is on the task PATH); tests/test_exomiser_report.py checks it.
Standard library only.
"""
import argparse
import csv
import gzip
import json
import os
import re
import sys

TERM = re.compile(r"^HP:[0-9]{7}$")
TERMS_HEADER = ["hpo_id", "label", "given_as"]
OUT_HEADER = ["rank", "gene", "moi", "exomiser_score", "phenotype_score", "variant_score", "p_value",
              "variant", "genotype", "consequence", "hgvs", "acmg", "acmg_evidence", "disease",
              "clinvar", "max_frequency"]
# Columns of Exomiser's variants TSV (TsvVariantResultsWriter, 15.x) the report reads.
NEED = ["#RANK", "GENE_SYMBOL", "MOI", "P-VALUE", "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE",
        "EXOMISER_GENE_VARIANT_SCORE", "CONTRIBUTING_VARIANT", "CONTIG", "START", "REF", "ALT", "GENOTYPE",
        "FUNCTIONAL_CLASS", "HGVS", "EXOMISER_ACMG_CLASSIFICATION", "EXOMISER_ACMG_EVIDENCE",
        "EXOMISER_ACMG_DISEASE_ID", "EXOMISER_ACMG_DISEASE_NAME", "CLINVAR_PRIMARY_INTERPRETATION",
        "CLINVAR_STAR_RATING", "MAX_FREQ_SOURCE", "MAX_FREQ"]
MOI_WORDS = {"AD": "autosomal dominant", "AR": "autosomal recessive", "XD": "X-linked dominant",
             "XR": "X-linked recessive", "MT": "mitochondrial", "ANY": "any mode of inheritance"}
# Exomiser's exome preset (examples/preset-exome-analysis.yml of Exomiser 15.1.0).
INHERITANCE = {"AUTOSOMAL_DOMINANT": 0.1, "AUTOSOMAL_RECESSIVE_HOM_ALT": 0.1, "AUTOSOMAL_RECESSIVE_COMP_HET": 2.0,
               "X_DOMINANT": 0.1, "X_RECESSIVE_HOM_ALT": 0.1, "X_RECESSIVE_COMP_HET": 2.0, "MITOCHONDRIAL": 0.2}
FREQUENCY_SOURCES = ["UK10K", "GNOMAD_E_AFR", "GNOMAD_E_AMR", "GNOMAD_E_EAS", "GNOMAD_E_NFE", "GNOMAD_E_SAS",
                     "GNOMAD_G_AFR", "GNOMAD_G_AMR", "GNOMAD_G_EAS", "GNOMAD_G_NFE", "GNOMAD_G_SAS"]
PATHOGENICITY_SOURCES = ["REVEL", "MVP", "ALPHA_MISSENSE", "SPLICE_AI"]
REMOVED_EFFECTS = ["FIVE_PRIME_UTR_EXON_VARIANT", "FIVE_PRIME_UTR_INTRON_VARIANT", "THREE_PRIME_UTR_EXON_VARIANT",
                   "THREE_PRIME_UTR_INTRON_VARIANT", "NON_CODING_TRANSCRIPT_EXON_VARIANT",
                   "NON_CODING_TRANSCRIPT_INTRON_VARIANT", "CODING_TRANSCRIPT_INTRON_VARIANT",
                   "UPSTREAM_GENE_VARIANT", "DOWNSTREAM_GENE_VARIANT", "INTERGENIC_VARIANT",
                   "REGULATORY_REGION_VARIANT"]
# Genes written to Exomiser's own files: the report shows --top of them.
OUTPUT_GENES = 50

NOTE = ("A ranking to bring to a clinical geneticist, not a diagnosis. A high score means a rare variant, "
        "predicted to be damaging, in a gene linked to symptoms like the ones given; it does not show that the "
        "variant causes them. A low score or an empty list does not rule a genetic cause out: this search reads "
        "coding and splice variants against the terms given, and not structural variants, repeat expansions, "
        "mitochondrial heteroplasmy or variants outside genes.")


class Fail(Exception):
    pass


def split_list(text):
    return [x for x in re.split(r"[,\s]+", text or "") if x]


# --- HPO terms ------------------------------------------------------------------

def read_obo(path):
    """{id: (name, obsolete, replaced_by or '', consider list)} and {alt_id: id} of an OBO file."""
    terms, alt, cur = {}, {}, None

    def close():
        if cur and cur.get("id", "").startswith("HP:"):
            terms[cur["id"]] = (cur.get("name", ""), cur.get("obsolete", False),
                                cur.get("replaced_by", ""), cur.get("consider", []))
            for a in cur.get("alt_id", []):
                alt[a] = cur["id"]

    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("["):
                close()
                cur = {} if line == "[Term]" else None
                continue
            if cur is None or ": " not in line:
                continue
            k, v = line.split(": ", 1)
            v = v.split(" ! ")[0].strip()
            if k in ("id", "name", "replaced_by"):
                cur[k] = v
            elif k == "is_obsolete":
                cur["obsolete"] = v == "true"
            elif k in ("alt_id", "consider"):
                cur.setdefault(k, []).append(v)
    close()
    return terms, alt


def check_terms(given, obo=None):
    """[(id, label, given_as)] of the given terms, in order, once each; Fail on
    a malformed, unknown or obsolete term."""
    if not given:
        raise Fail("no HPO term given: name the symptoms as HPO ids, e.g. HPO=HP:0001250,HP:0001263 "
                   "(find them at https://hpo.jax.org)")
    bad = [t for t in given if not TERM.match(t)]
    if bad:
        raise Fail("not an HPO id: %s. An HPO id is HP: and seven digits, e.g. HP:0001250 (Seizure)"
                   % ", ".join(bad))
    terms, alt = read_obo(obo) if obo else ({}, {})
    out, seen, problems = [], set(), []
    for t in given:
        tid = alt.get(t, t)
        if obo:
            if tid not in terms:
                problems.append("%s is not a term of the HPO release in %s" % (t, obo))
                continue
            name, obsolete, repl, consider = terms[tid]
            if obsolete:
                hint = (" use %s instead" % repl) if repl else (
                    (" consider %s" % ", ".join(consider)) if consider else "")
                problems.append("%s (%s) is obsolete in this HPO release;%s" % (t, name, hint or " pick a current term"))
                continue
        else:
            name = ""
        if tid in seen:
            continue
        seen.add(tid)
        out.append((tid, name, t))
    if problems:
        raise Fail("; ".join(problems) + ". Look the terms up at https://hpo.jax.org")
    return out


# --- the family file ------------------------------------------------------------

def vcf_samples(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    raise Fail("%s has no #CHROM header line" % path)


def read_pedigree(path):
    """Rows of step 38's pedigree TSV (sample, sex, father, mother, partner)."""
    rows, header = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if header is None:
                header = [c.strip().lower() for c in line.lstrip("#").split("\t")]
                if "sample" not in header:
                    raise Fail("%s: the first line must be a header with a 'sample' column" % path)
                continue
            if line.startswith("#"):
                continue
            vals = [v.strip() for v in line.split("\t")]
            vals += [""] * (len(header) - len(vals))
            r = dict(zip(header, vals))
            rows[r["sample"]] = r
    if header is None:
        raise Fail("%s is empty" % path)
    return rows


def write_ped(pedigree, members, proband, affected, unaffected, out):
    """The PED file Exomiser reads: one line per member of the VCF."""
    rows = read_pedigree(pedigree)
    missing = [m for m in members if m not in rows]
    if missing:
        raise Fail("the family VCF holds %s, which %s does not name" % (", ".join(missing), pedigree))
    if proband not in members:
        raise Fail("%s is not a sample of the family VCF (%s)" % (proband, ", ".join(members)))
    for label, names in (("--affected", affected), ("--unaffected", unaffected)):
        stray = [n for n in names if n not in members]
        if stray:
            raise Fail("%s names %s, not a member of the family VCF (%s)"
                       % (label, ", ".join(stray), ", ".join(members)))
    both = sorted(set(affected) & set(unaffected))
    if both:
        raise Fail("%s is named both affected and unaffected" % ", ".join(both))
    if proband in unaffected:
        raise Fail("%s is the proband, the person the symptoms are of: it cannot be unaffected" % proband)
    sex = {"male": "1", "female": "2"}
    with open(out, "w") as f:
        for m in members:
            r = rows[m]
            father = r.get("father", "") if r.get("father", "") in members else "0"
            mother = r.get("mother", "") if r.get("mother", "") in members else "0"
            status = "2" if m == proband or m in affected else "1" if m in unaffected else "0"
            f.write("\t".join(["family", m, father, mother, sex.get(r.get("sex", "").lower(), "0"), status]) + "\n")


def cmd_prepare(a):
    terms = check_terms(split_list(a.hpo), a.obo)
    analysis = {
        "genomeAssembly": a.assembly,
        "vcf": a.vcf,
        "hpoIds": [t[0] for t in terms],
        "inheritanceModes": INHERITANCE,
        "analysisMode": "PASS_ONLY",
        "frequencySources": FREQUENCY_SOURCES,
        "pathogenicitySources": PATHOGENICITY_SOURCES,
        "steps": [
            {"geneBlacklistFilter": {}},
            {"failedVariantFilter": {}},
            {"alleleBalanceFilter": {}},
            {"variantEffectFilter": {"remove": REMOVED_EFFECTS}},
            {"frequencyFilter": {"maxFrequency": 2.0}},
            {"pathogenicityFilter": {"keepNonPathogenic": True}},
            {"inheritanceFilter": {}},
            {"omimPrioritiser": {}},
            {"hiPhivePrioritiser": {}},
        ],
    }
    if a.pedigree:
        if not a.ped:
            raise Fail("--pedigree needs --ped, the PED file to write")
        write_ped(a.pedigree, vcf_samples(a.vcf), a.sample, split_list(a.affected), split_list(a.unaffected), a.ped)
        analysis["ped"] = a.ped
        analysis["proband"] = a.sample
    elif split_list(a.affected) or split_list(a.unaffected):
        raise Fail("--affected and --unaffected name family members: they need --pedigree and the family VCF")
    else:
        n = len(vcf_samples(a.vcf))
        if n != 1:
            raise Fail("%s has %d samples; one person's VCF has one (a family VCF goes with --pedigree)" % (a.vcf, n))
    doc = {"analysis": analysis,
           "outputOptions": {"outputContributingVariantsOnly": False, "numGenes": OUTPUT_GENES,
                             "outputDirectory": a.out_dir, "outputFileName": "%s_exomiser" % a.sample,
                             "outputFormats": ["TSV_VARIANT", "TSV_GENE", "HTML"]}}
    with open(a.analysis, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    with open(a.terms, "w") as f:
        f.write("\t".join(TERMS_HEADER) + "\n")
        for t in terms:
            f.write("\t".join(t) + "\n")
    print("HPO terms: " + ", ".join("%s (%s)" % (t[0], t[1]) if t[1] else t[0] for t in terms))


# --- the report -----------------------------------------------------------------

def read_terms(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or any(k not in rows[0] for k in TERMS_HEADER):
        raise Fail("%s has no term rows or lacks a column of %s" % (path, ", ".join(TERMS_HEADER)))
    return rows


def number(v):
    try:
        return "%.3f" % float(v)
    except (TypeError, ValueError):
        return "."


def genotype_word(gt):
    g = gt.replace("|", "/")
    if g in ("0/1", "1/0"):
        return "heterozygous"
    if g == "1/1":
        return "homozygous"
    if g == "1":
        return "hemizygous"
    return gt


def read_variants(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
        header = rows[0].keys() if rows else None
    if header is None:
        with open(path, encoding="utf-8", errors="replace") as f:
            first = f.readline().rstrip("\n").split("\t")
        if not first or first[0] != "#RANK":
            raise Fail("%s is not an Exomiser variants TSV (no #RANK header)" % path)
        header = first
    missing = [k for k in NEED if k not in header]
    if missing:
        raise Fail("%s lacks the column%s %s: is it Exomiser 15's TSV_VARIANT output?"
                   % (path, "s" if len(missing) > 1 else "", ", ".join(missing)))
    bad = next((i for i, r in enumerate(rows, 2) if None in r or any(r.get(k) is None for k in NEED)), None)
    if bad:
        raise Fail("%s line %d does not have one field per column" % (path, bad))
    return rows


def top_rows(rows, top):
    """The contributing variants of the first `top` ranks, in rank order."""
    out, ranks = [], []
    for r in rows:
        try:
            rank = int(r["#RANK"])
        except ValueError:
            raise Fail("rank '%s' is not a number" % r["#RANK"])
        if r["CONTRIBUTING_VARIANT"] != "1":
            continue
        if rank not in ranks:
            if len(ranks) >= top:
                continue
            ranks.append(rank)
        disease = r["EXOMISER_ACMG_DISEASE_NAME"]
        if r["EXOMISER_ACMG_DISEASE_ID"]:
            disease = "%s (%s)" % (disease, r["EXOMISER_ACMG_DISEASE_ID"]) if disease else r["EXOMISER_ACMG_DISEASE_ID"]
        clinvar = r["CLINVAR_PRIMARY_INTERPRETATION"]
        clinvar = "." if clinvar in ("", "NOT_PROVIDED") else "%s (%s star%s)" % (
            clinvar, r["CLINVAR_STAR_RATING"] or "0", "" if r["CLINVAR_STAR_RATING"] == "1" else "s")
        freq = "%s %s%%" % (r["MAX_FREQ_SOURCE"], r["MAX_FREQ"]) if r["MAX_FREQ"] else "not seen"
        out.append({"rank": str(rank), "gene": r["GENE_SYMBOL"], "moi": r["MOI"],
                    "exomiser_score": number(r["EXOMISER_GENE_COMBINED_SCORE"]),
                    "phenotype_score": number(r["EXOMISER_GENE_PHENO_SCORE"]),
                    "variant_score": number(r["EXOMISER_GENE_VARIANT_SCORE"]),
                    "p_value": r["P-VALUE"] or ".",
                    "variant": "%s:%s %s>%s" % (r["CONTIG"], r["START"], r["REF"], r["ALT"]),
                    "genotype": r["GENOTYPE"] or ".", "consequence": r["FUNCTIONAL_CLASS"] or ".",
                    "hgvs": r["HGVS"] or ".", "acmg": r["EXOMISER_ACMG_CLASSIFICATION"] or ".",
                    "acmg_evidence": r["EXOMISER_ACMG_EVIDENCE"] or ".", "disease": disease or ".",
                    "clinvar": clinvar, "max_frequency": freq})
    return out


def words(rows, terms):
    """The text file: what was searched, each rank in a sentence, then the note."""
    named = ", ".join("%s (%s)" % (t["label"], t["hpo_id"]) if t["label"] else t["hpo_id"] for t in terms)
    lines = ["Symptoms searched: %s." % named]
    if not rows:
        lines.append("No gene passed Exomiser's filters with a variant that fits a mode of inheritance.")
    by_rank = {}
    for r in rows:
        by_rank.setdefault(r["rank"], []).append(r)
    for rank, vs in by_rank.items():
        g = vs[0]
        head = "%s. %s (%s): Exomiser score %s (symptom match %s, variant %s)." % (
            rank, g["gene"], MOI_WORDS.get(g["moi"], g["moi"]), g["exomiser_score"], g["phenotype_score"],
            g["variant_score"])
        parts = []
        for v in vs:
            p = "%s %s, %s %s, ACMG %s" % (v["variant"], genotype_word(v["genotype"]), v["consequence"], v["hgvs"],
                                          v["acmg"].replace("_", " ").lower())
            if v["acmg_evidence"] != ".":
                p += " (%s)" % v["acmg_evidence"]
            if v["disease"] != ".":
                p += ", linked disease %s" % v["disease"]
            if v["clinvar"] != ".":
                p += ", ClinVar %s" % v["clinvar"].replace("_", " ").lower()
            p += ", highest population frequency %s" % v["max_frequency"]
            parts.append(p)
        lines.append(head + " " + "; ".join(parts) + ".")
    lines.append(NOTE)
    return lines


def cmd_report(a):
    terms = read_terms(a.terms)
    rows = top_rows(read_variants(a.variants), a.top)
    with open(a.out, "w") as f:
        f.write("\t".join(OUT_HEADER) + "\n")
        for r in rows:
            f.write("\t".join(r[k].replace("\t", " ") for k in OUT_HEADER) + "\n")
    text = words(rows, terms)
    with open(a.text, "w") as f:
        f.write("\n".join(text) + "\n")
    print("\n".join(text[:1 + min(len(rows), a.top)]))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("prepare")
    q.add_argument("--hpo", required=True)
    q.add_argument("--obo")
    q.add_argument("--sample", required=True)
    q.add_argument("--vcf", required=True)
    q.add_argument("--assembly", default="hg38", choices=["hg38", "hg19"])
    q.add_argument("--out-dir", required=True)
    q.add_argument("--analysis", required=True)
    q.add_argument("--terms", required=True)
    q.add_argument("--pedigree")
    q.add_argument("--ped")
    q.add_argument("--affected", default="")
    q.add_argument("--unaffected", default="")
    r = sub.add_parser("report")
    r.add_argument("--sample", required=True)
    r.add_argument("--variants", required=True)
    r.add_argument("--terms", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--text", required=True)
    r.add_argument("--top", type=int, default=10)
    a = p.parse_args(argv)
    try:
        (cmd_prepare if a.cmd == "prepare" else cmd_report)(a)
    except (Fail, OSError) as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
